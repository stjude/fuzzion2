CXX=g++
CXXFLAGS=-std=c++11 -O3

SRC_PREFIX=src/cpp
BUILD_PREFIX=build
BIN_PREFIX=$(BUILD_PREFIX)/bin
OBJ_PREFIX=$(BUILD_PREFIX)/obj

FUZZION2_NAME=fuzzion2
FUZZION2_BIN=$(BIN_PREFIX)/$(FUZZION2_NAME)
FUZZION2_SRC_BASENAMES=fuzzion2.cpp bamread.cpp bin.cpp fastq.cpp kmer.cpp \
	minimizer.cpp pairread.cpp pattern.cpp rank.cpp refgen.cpp ubam.cpp \
	util.cpp
FUZZION2_SRCS=$(FUZZION2_SRC_BASENAMES:%.cpp=$(SRC_PREFIX)/%.cpp)
FUZZION2_OBJS=$(FUZZION2_SRC_BASENAMES:%.cpp=$(OBJ_PREFIX)/%.o)
FUZZION2_LDLIBS=-lhts -lpthread

KMERANK_NAME=kmerank
KMERANK_BIN=$(BIN_PREFIX)/$(KMERANK_NAME)
KMERANK_SRC_BASENAMES=kmerank.cpp bin.cpp kmer.cpp rank.cpp refgen.cpp util.cpp
KMERANK_SRCS=$(KMERANK_SRC_BASENAMES:%.cpp=$(SRC_PREFIX)/%.cpp)
KMERANK_OBJS=$(KMERANK_SRC_BASENAMES:%.cpp=$(OBJ_PREFIX)/%.o)

.PHONY: all clean

all: $(FUZZION2_BIN) $(KMERANK_BIN)

$(FUZZION2_BIN): $(FUZZION2_OBJS) | $(BIN_PREFIX)
	$(CXX) $^ $(CXXFLAGS) $(FUZZION2_LDLIBS) -o $@

$(FUZZION2_NAME): $(FUZZION2_BIN)

$(KMERANK_BIN): $(KMERANK_OBJS) | $(BIN_PREFIX)
	$(CXX) $^ $(CXXFLAGS) -o $@

$(KMERANK_NAME): $(KMERANK_BIN)

$(OBJ_PREFIX)/%.o: $(SRC_PREFIX)/%.cpp | $(OBJ_PREFIX)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_PREFIX) $(OBJ_PREFIX):
	mkdir -p $@

clean:
	rm -r $(BUILD_PREFIX)
