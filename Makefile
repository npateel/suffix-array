
CC = g++

CFLAGS = -std=c++2a -g -Wall -Wno-deprecated-declarations -I$(HOME)/include -Iinclude -L$(HOME)/lib

LDFLAGS = -lsdsl -ldivsufsort -ldivsufsort64


.PHONY = all clean

QUERY = bin/querysa
BUILD = bin/buildsa

TARGETS = $(QUERY) $(BUILD)

all: $(TARGETS)

$(QUERY): build/querysa.o
	$(CC) $(CFLAGS) -o $(QUERY) build/querysa.o $(LDFLAGS)



$(BUILD): build/buildsa.o
	$(CC) $(CFLAGS) -o $(BUILD) build/buildsa.o $(LDFLAGS)

build/%.o: src/%.cpp
	$(CC) $(CFLAGS) $< -c -o $@ $(LDFLAGS)



clean:
	rm -rf build/*.o $(TARGET)
