#include <sdsl/csa_bitcompressed.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
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
static char args_doc[] = "REFERENCE OUTPUT";

/* The options we understand. */
static struct argp_option options[] = {
    {"preftab", 777, "k", 0, "Prefix table length"},
    {0, 'b', "benchmarking_file", 0, "Path to benchmarking file"},
    {0}};

/* Used by main to communicate with parse_opt. */
struct arguments {
  int preftab;
  char *reference_file;
  char *output_file;
  char *benchmarking_file;
};

/* Parse a single option. */
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = (struct arguments *)state->input;

  switch (key) {
    case 777: {
      std::string a(arg);
      arguments->preftab = std::stoi(a);
      break;
    }
    case 'b': {
      arguments->benchmarking_file = arg;
    }

    case ARGP_KEY_ARG:
      if (state->arg_num >= 2) /* Too many arguments. */
        argp_usage(state);
      if (state->arg_num == 0) {
        arguments->reference_file = arg;
      } else if (state->arg_num == 1) {
        arguments->output_file = arg;
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

/* Our argp parser. */
static struct argp argp = {options, parse_opt, args_doc, doc};

int main(int argc, char **argv) {
  struct arguments arguments;
  arguments.preftab = -1;

  /* Parse our arguments; every option seen by parse_opt will
     be reflected in arguments. */
  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  std::ifstream ref(arguments.reference_file);
  bool benchmarking = (arguments.benchmarking_file != NULL);
  std::ofstream bfile(arguments.benchmarking_file, std::ofstream::app);

  std::string seq;
  bool firstLine = true;
  if (ref.is_open()) {
    std::string line;
    while (std::getline(ref, line)) {
      // if (line.at(0) != '>') {
      if (!firstLine) {
        seq += line;
      }
      firstLine = false;
    }
    ref.close();

  } else {
    // file dont exist :(
    exit(1);
  }

  // work file
  std::string work = "work";
  std::ofstream workfile(work);
  workfile << seq.c_str();
  workfile.close();
  if (benchmarking) bfile << seq.length();
  if (benchmarking) bfile << "," << arguments.preftab;

  sdsl::cache_config cc(
      false);  // do not delete temp files after csa construction;

  sdsl::csa_wt<> csa;

  auto start = std::chrono::steady_clock::now();
  sdsl::construct(csa, work, 1);
  auto end = std::chrono::steady_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
          .count() /
      1.0e9;

  std::cout << "Suffix Construction Time for file " << arguments.reference_file
            << " was " << duration << std::endl;
  if (benchmarking) bfile << "," << duration;

  // prefix -> (startidx, endindex)
  std::unordered_map<std::string, std::pair<int, int>> prefix_table;
  if (arguments.preftab != -1) {
    // 0 is $
    start = std::chrono::steady_clock::now();
    std::string prefix = seq.substr(csa[1], arguments.preftab);
    if (prefix.length() == (uint)arguments.preftab) {
      prefix_table[prefix].first = 1;
    }

    std::string current;
    for (uint i = 2; i < csa.size(); i++) {
      current = seq.substr(csa[i], arguments.preftab);
      // std::cout << current << std::endl;
      //  found our first diff sequence
      if (current.compare(prefix) != 0) {
        if (prefix.length() == (uint)arguments.preftab) {
          prefix_table[prefix].second = i - 1;
        }
        if (current.length() == (uint)arguments.preftab) {
          prefix_table[current].first = i;
          prefix_table[current].second = i;
        }
      }
      prefix = current;
    }
    prefix_table[current].second = csa.size() - 1;
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1.0e9;
    std::cout << "Preftable Construction Time for file "
              << arguments.reference_file << " was " << duration << std::endl;
  }
  if (benchmarking) {
    bfile << ",";
    if (arguments.preftab != -1) {
      bfile << duration;
    }
  }

  //   for (uint i = 0; i < csa.size(); i++) {
  //     std::cout << std::to_string(csa[i]) + " , ";
  //   }
  //   std::cout << std::endl;

  // for (const auto& [key, val] : prefix_table) {
  //   std::cout << "k: " << key << ",  val: " << std::to_string(val.first) <<
  //   "," << std::to_string(val.second) << std::endl;
  // }
  // std::cout << std::endl;

  // archive everything
  std::ofstream outfile(arguments.output_file);
  // to preserver order of cereal data
  {
    cereal::BinaryOutputArchive oarchive(outfile);
    oarchive(prefix_table, seq);
  }
  // cereal went out of scope, so contents are flushed.
  csa.serialize(outfile);
  outfile.close();

  if (benchmarking)
    bfile << "," << std::filesystem::file_size(arguments.output_file)
          << std::endl;
}