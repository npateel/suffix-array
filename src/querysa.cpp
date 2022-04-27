#include <sdsl/csa_bitcompressed.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>

#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <stdio.h>
#include <argp.h>

const char *argp_program_version = "buildsa 1.0";
const char *argp_program_bug_address = "<npateel@terpmail.umd.edu>";

/* Program documentation. */
static char doc[] =
    "buildsa -- Builds a suffix array given a reference genome (FASTA) and "
    "output location";

/* A description of the arguments we accept. */
static char args_doc[] = "INDEX QUERYFILE QUERYMODE OUTPUT";

/* The options we understand. */
static struct argp_option options[] = {
    {0, 'b', "benchmarking_file", 0, "Path to benchmarking file"}, {0}};

/* Used by main to communicate with parse_opt. */
struct arguments {
  char *index;
  char *queries;
  char *query_mode;
  char *output;
  char *benchmarking_file;
};

/* Parse a single option. */
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = (struct arguments *)state->input;

  switch (key) {
    case 'b':
      arguments->benchmarking_file = arg;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 4) /* Too many arguments. */
        argp_usage(state);
      if (state->arg_num == 0) {
        arguments->index = arg;
      } else if (state->arg_num == 1) {
        arguments->queries = arg;
      } else if (state->arg_num == 2) {
        if (strcmp(arg, "naive") == 0 || strcmp(arg, "simpaccel") == 0) {
          arguments->query_mode = arg;
        } else {
          argp_usage(state);
          break;
        }
      } else if (state->arg_num == 3) {
        arguments->output = arg;
      }
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 2) /* Not enough arguments. */
        argp_usage(state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

int lcp(std::string const &s1, std::string const &s2) {
  if (s1.length() > s2.length()) {
    return lcp(s2, s1);
  }
  // s1 guaranteed smaller now
  int prefix = 0;
  for (uint i = 0; i < s1.length(); i++) {
    if (s1.at(i) != s2.at(i)) {
      return prefix;
    }
    prefix++;
  }
  return prefix;
}

/* Our argp parser. */
static struct argp argp = {options, parse_opt, args_doc, doc};

// not only does this do a comparison, it will return the lcp length of query
// with seq at (seqidx + minlcp)
int lcpcompare(const std::string &seq, int seqidx, const std::string &query,
               int minlcp) {
  std::string segment = seq.substr(seqidx, query.length() - minlcp);
  int lcp = minlcp;
  for (uint i = (uint)minlcp;
       i < query.length() && (uint)seqidx + (uint)i < seq.length(); i++) {
    if (query.at(i) == seq.at(seqidx + i)) {
      lcp++;
    } else {
      if (query.at(i) < seq.at(seqidx + i)) {
        // query string is before, return lcp -1
        return lcp * -1 - 1;
      } else {
        // query string is after, return lcp + 1
        return lcp + 1;
      }
    }
  }
  // query string is equal
  if (lcp == (int)query.length()) {
    return 0;
  } else {
    // query string comes after
    return 1;
  }
}

void lcpsearch(int startidx, int endidx, std::string const &seq,
               sdsl::csa_wt<> &csa, std::string const &name,
               std::string const &query,
               std::unordered_map<std::string, std::pair<int, int>> &results,
               std::unordered_map<std::string, double> &times) {
  bool found = false;
  int smallest, largest;
  int start = startidx;
  int end = endidx;
  int startlcp = lcp(query, seq.substr(csa[start], query.length()));
  int endlcp = lcp(query, seq.substr(csa[end - 1], query.length()));
  auto starttime = std::chrono::steady_clock::now();
  int minlcp;
  while (start <= end) {
    minlcp = std::min(startlcp, endlcp);
    int mid = (start + end) / 2;
    int compare = lcpcompare(seq, csa[mid], query, minlcp);
    // std::cout << start << "," << end << ", q: " << query << ", startlcp: " <<
    // startlcp << ", endlcp: " << endlcp << ", comapre: "<< compare  <<
    // std::endl;
    //  query < mid
    if (compare < 0) {
      end = mid - 1;
      endlcp = compare * -1 - 1;

    } else if (compare > 0) {
      start = mid + 1;
      startlcp = compare - 1;

    } else {
      // im too lazy to do a speedup here. Probably could help but \_(:/)_/
      // shouldn't be a majority of cases. Assuming that after
      found = true;
      if (mid <= startidx || query.compare(seq.substr(csa[mid - 1], query.length())) > 0) {
        // we know we have the first one
        smallest = mid;
        break;
      }
      end = mid - 1;
      endlcp = query.length();
    }
  }
  // do the same thing for largest index
  start = startidx;
  end = endidx;
  startlcp = lcp(query, seq.substr(csa[start], query.length()));
  endlcp = lcp(query, seq.substr(csa[end - 1], query.length()));
  while (start <= end) {
    minlcp = std::min(startlcp, endlcp);
    int mid = (start + end) / 2;
    int compare = lcpcompare(seq, csa[mid], query, minlcp);
    // std::cout << start << ",  " << end << ", q: " << query << ", startlcp: "
    // << startlcp << ", endlcp: " << endlcp << ", comapre: " << compare  <<
    // std::endl;
    //  query < mid
    if (compare < 0) {
      end = mid - 1;
      endlcp = compare * -1 - 1;

    } else if (compare > 0) {
      start = mid + 1;
      startlcp = compare - 1;
    } else {
      if (mid >= endidx - 1 ||
          query.compare(seq.substr(csa[mid + 1], query.length())) < 0) {
        // we know we have the last one
        largest = mid;
        break;
      }
      start = mid + 1;
      startlcp = query.length();
    }
  }
  auto endtime = std::chrono::steady_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::nanoseconds>(endtime - starttime)
          .count();
  if (found) {
    // std::cout << "smallest: " << smallest << " largest: " << largest <<
    // std::endl;
    results[name] = {smallest, largest};
  } else {
    results[name] = {-1, -2};
  }
  times[name] = duration;
}

void binsearch(int startidx, int endidx, std::string const &seq,
               sdsl::csa_wt<> &csa, std::string const &name,
               std::string const &query,
               std::unordered_map<std::string, std::pair<int, int>> &results,
               std::unordered_map<std::string, double> &times) {
  bool found = false;
  int smallest, largest;
  int start = startidx;
  int end = endidx;
  auto starttime = std::chrono::steady_clock::now();
  while (start <= end) {
    int mid = (start + end) / 2;
    std::string segment = seq.substr(csa[mid], query.length());
    // std::cout << start << "," << end << ", q: " << query << ", segment: " <<
    // segment << std::endl;
    int compare = query.compare(segment);
    // query < mid
    if (compare < 0) {
      end = mid - 1;

    } else if (compare > 0) {
      start = mid + 1;

    } else {
      found = true;
      if (mid == startidx ||
          query.compare(seq.substr(csa[mid - 1], query.length())) > 0) {
        // we know we have the first one
        smallest = mid;
        break;
      }
      end = mid - 1;
    }
  }
  // do the same thing for largest index
  start = smallest;
  end = endidx;
  while (start <= end) {
    int mid = (start + end) / 2;
    std::string segment = seq.substr(csa[mid], query.length());
    // std::cout << start << "," << end << ", q: " << query << ", segment: " <<
    // segment << std::endl;
    int compare = query.compare(segment);
    // query < mid
    if (compare < 0) {
      end = mid - 1;

    } else if (compare > 0) {
      start = mid + 1;

    } else {
      if (mid == endidx - 1 ||
          query.compare(seq.substr(csa[mid + 1], query.length())) < 0) {
        // we know we have the last one
        largest = mid;
        break;
      }
      start = mid + 1;
    }
  }
  auto endtime = std::chrono::steady_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::nanoseconds>(endtime - starttime)
          .count();
  if (found) {
    /// std::cout << "smallest: " << smallest << " largest: " << largest <<
    /// std::endl;
    results[name] = {smallest, largest};
  } else {
    results[name] = {-1, -2};
  }
  times[name] = duration;
}

void naive(std::unordered_map<std::string, std::pair<int, int>> &prefix_table,
           std::string &seq, sdsl::csa_wt<> &csa,
           std::unordered_map<std::string, std::string> &queries, int k,
           std::unordered_map<std::string, std::pair<int, int>> &results,
           std::unordered_map<std::string, double> &times) {
  // set starting and ending positions
  int start = 0;
  int end = csa.size() - 1;
  for (const auto &[name, query] : queries) {
    // now do binary search!
    binsearch(start, end, seq, csa, name, query, results, times);
  }
}

void naiveprefix(
    std::unordered_map<std::string, std::pair<int, int>> &prefix_table,
    std::string &seq, sdsl::csa_wt<> &csa,
    std::unordered_map<std::string, std::string> &queries, int k,
    std::unordered_map<std::string, std::pair<int, int>> &results,
    std::unordered_map<std::string, double> &times) {
  // set starting and ending positions
  for (const auto &[name, query] : queries) {
    std::string prefix = query.substr(0, k);
    int start, end;
    //if (prefix_table.find(prefix) != prefix_table.end()) {
      start = prefix_table[prefix].first;
      end = prefix_table[prefix].second;
    //} else {
      //start = 0;
      //end = csa.size() ;
    //}

    // now do binary search!
    binsearch(start, end, seq, csa, name, query, results, times);
  }
}

void lcpnoprefix(
    std::unordered_map<std::string, std::pair<int, int>> &prefix_table,
    std::string &seq, sdsl::csa_wt<> &csa,
    std::unordered_map<std::string, std::string> &queries, int k,
    std::unordered_map<std::string, std::pair<int, int>> &results,
    std::unordered_map<std::string, double> &times) {
  int start = 0;
  int end = csa.size() - 1;
  for (const auto &[name, query] : queries) {
    // now do binary search!
    lcpsearch(start, end, seq, csa, name, query, results, times);
  }
}
void lcpprefix(
    std::unordered_map<std::string, std::pair<int, int>> &prefix_table,
    std::string &seq, sdsl::csa_wt<> &csa,
    std::unordered_map<std::string, std::string> &queries, int k,
    std::unordered_map<std::string, std::pair<int, int>> &results,
    std::unordered_map<std::string, double> &times) {
  for (const auto &[name, query] : queries) {
    std::string prefix = query.substr(0, k);
    int start, end;
    //if (prefix_table.find(prefix) != prefix_table.end()) {
      start = prefix_table[prefix].first;
      end = prefix_table[prefix].second;
      if (end == 0) {
        end = csa.size() - 1;
      }
    //} else {
      //start = 0;
      //end = csa.size();
    //}
    //std::cout << "start " << start << " end " << end << std::endl;
    // now do binary search!
    lcpsearch(start, end, seq, csa, name, query, results, times);
  }
}

int main(int argc, char **argv) {
  struct arguments arguments;

  /* Parse our arguments; every option seen by parse_opt will
     be reflected in arguments. */
  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  // load index file

  std::unordered_map<std::string, std::pair<int, int>> prefix_table;
  std::string seq;
  sdsl::csa_wt<> csa;
  std::ifstream infile(arguments.index);
  {
    cereal::BinaryInputArchive iarchive(infile);
    iarchive(prefix_table, seq);
  }

  bool bench = (arguments.benchmarking_file != NULL);
  std::ofstream bfile(arguments.benchmarking_file, std::ofstream::app);

// for (const auto& [key, val] : prefix_table) {
//      std::cout << "k: " << key << ",  val: " << std::to_string(val.first) <<
//      "," << std::to_string(val.second) << std::endl;
//    }

  csa.load(infile);
  if (bench) bfile << arguments.index << "," << arguments.query_mode;


  // determine size of k
  int k = -1;
  if (prefix_table.size() > 0) {
    k = prefix_table.begin()->first.length();
  }

  //std::cout << "loading queries" << std::endl;
  // load in queries
  std::unordered_map<std::string, std::string> queries;
  std::vector<std::string> listofquerynames;
  std::string queryname;
  std::ifstream queryfile(arguments.queries);
  if (queryfile.is_open()) {
    std::string line;
    while (std::getline(queryfile, line)) {
      if (line.at(0) == '>') {
        queryname = line.substr(1);
        queries[queryname] = "";
        listofquerynames.push_back(queryname);
      } else {
        queries[queryname].append(line);
      }
    }
    queryfile.close();
  } else {
    exit(1);
  }

// for (const auto& [key, val] : queries) {
//      std::cout << "k: " << key << ",  val: " << val << std::endl;
//    }



  if (bench) bfile << "," << queries.begin()->second.length();
  //std::cout << "computing results" << std::endl;
  // get results

  std::unordered_map<std::string, std::pair<int, int>> results;
  std::unordered_map<std::string, double> times;

  if (strcmp(arguments.query_mode, "simpaccel") == 0) {
    if (k == -1) {
      lcpnoprefix(prefix_table, seq, csa, queries, k, results, times);

    } else {
      lcpprefix(prefix_table, seq, csa, queries, k, results, times);
    }
  } else {
    if (k == -1) {
      naive(prefix_table, seq, csa, queries, k, results, times);
    } else {
      naiveprefix(prefix_table, seq, csa, queries, k, results, times);
    }
  }

  //std::cout << "serializing results" << std::endl;
  // serialize results
  std::ofstream outputfile(arguments.output);
  for (std::string queryname : listofquerynames) {
    std::pair<int, int> positions = results[queryname];
    int numpositions = positions.second - positions.first + 1;
    outputfile << queryname << "\t" << numpositions;
    if (positions.first != -1) {
      for (int pos = positions.first; pos <= positions.second; pos++) {
        outputfile << "\t" << csa[pos];
      }
    }
    outputfile << std::endl;
  }
  outputfile.close();
  if (bench) {
    double sum = 0;
    for (const auto &[k, v] : times) {
      sum += v;
    }
    double avg = sum / times.size();
    bfile << "," << avg << "\n";
    bfile.close();
  }

  exit(0);
}
