CXX=g++
CXXFLAGS=-std=c++11 -O3

SRC_PREFIX=src/cpp
BUILD_PREFIX=build
BIN_PREFIX=$(BUILD_PREFIX)/bin
OBJ_PREFIX=$(BUILD_PREFIX)/obj

FUZZION2_NAME=fuzzion2
FUZZION2_BIN=$(BIN_PREFIX)/$(FUZZION2_NAME)
FUZZION2_SRC_BASENAMES=fuzzion2.cpp bamread.cpp bin.cpp fastq.cpp hit.cpp kmer.cpp \
	match.cpp minimizer.cpp pairread.cpp pattern.cpp rank.cpp refgen.cpp \
        ubam.cpp util.cpp window.cpp
FUZZION2_SRCS=$(FUZZION2_SRC_BASENAMES:%.cpp=$(SRC_PREFIX)/%.cpp)
FUZZION2_OBJS=$(FUZZION2_SRC_BASENAMES:%.cpp=$(OBJ_PREFIX)/%.o)
FUZZION2_LDLIBS=-lhts -lpthread

FUZZALL_NAME=fuzzall
FUZZALL_BIN=$(BIN_PREFIX)/$(FUZZALL_NAME)
FUZZALL_SRC_BASENAMES=fuzzall.cpp hit.cpp kmer.cpp minimizer.cpp pattern.cpp \
	summary.cpp util.cpp window.cpp
FUZZALL_SRCS=$(FUZZALL_SRC_BASENAMES:%.cpp=$(SRC_PREFIX)/%.cpp)
FUZZALL_OBJS=$(FUZZALL_SRC_BASENAMES:%.cpp=$(OBJ_PREFIX)/%.o)

FUZZION2HTML_NAME=fuzzion2html
FUZZION2HTML_BIN=$(BIN_PREFIX)/$(FUZZION2HTML_NAME)
FUZZION2HTML_SRC_BASENAMES=fuzzion2html.cpp group.cpp hit.cpp kmer.cpp \
	minimizer.cpp pattern.cpp summary.cpp util.cpp window.cpp
FUZZION2HTML_SRCS=$(FUZZION2HTML_SRC_BASENAMES:%.cpp=$(SRC_PREFIX)/%.cpp)
FUZZION2HTML_OBJS=$(FUZZION2HTML_SRC_BASENAMES:%.cpp=$(OBJ_PREFIX)/%.o)

FUZZORT_NAME=fuzzort
FUZZORT_BIN=$(BIN_PREFIX)/$(FUZZORT_NAME)
FUZZORT_SRC_BASENAMES=fuzzort.cpp hit.cpp kmer.cpp minimizer.cpp \
	pattern.cpp util.cpp window.cpp
FUZZORT_SRCS=$(FUZZORT_SRC_BASENAMES:%.cpp=$(SRC_PREFIX)/%.cpp)
FUZZORT_OBJS=$(FUZZORT_SRC_BASENAMES:%.cpp=$(OBJ_PREFIX)/%.o)

FUZZUM_NAME=fuzzum
FUZZUM_BIN=$(BIN_PREFIX)/$(FUZZUM_NAME)
FUZZUM_SRC_BASENAMES=fuzzum.cpp group.cpp hit.cpp kmer.cpp minimizer.cpp \
	pattern.cpp summary.cpp util.cpp window.cpp
FUZZUM_SRCS=$(FUZZUM_SRC_BASENAMES:%.cpp=$(SRC_PREFIX)/%.cpp)
FUZZUM_OBJS=$(FUZZUM_SRC_BASENAMES:%.cpp=$(OBJ_PREFIX)/%.o)

KMERANK_NAME=kmerank
KMERANK_BIN=$(BIN_PREFIX)/$(KMERANK_NAME)
KMERANK_SRC_BASENAMES=kmerank.cpp bin.cpp kmer.cpp rank.cpp refgen.cpp util.cpp
KMERANK_SRCS=$(KMERANK_SRC_BASENAMES:%.cpp=$(SRC_PREFIX)/%.cpp)
KMERANK_OBJS=$(KMERANK_SRC_BASENAMES:%.cpp=$(OBJ_PREFIX)/%.o)

.PHONY: all clean

all: $(FUZZION2_BIN) $(FUZZALL_BIN) $(FUZZION2HTML_BIN) $(FUZZORT_BIN) $(FUZZUM_BIN) $(KMERANK_BIN)

$(FUZZION2_BIN): $(FUZZION2_OBJS) | $(BIN_PREFIX)
	$(CXX) $^ $(CXXFLAGS) $(FUZZION2_LDLIBS) -o $@
$(FUZZION2_NAME): $(FUZZION2_BIN)

$(FUZZALL_BIN): $(FUZZALL_OBJS) | $(BIN_PREFIX)
	$(CXX) $^ $(CXXFLAGS) -o $@
$(FUZZALL_NAME): $(FUZZALL_BIN)

$(FUZZION2HTML_BIN): $(FUZZION2HTML_OBJS) | $(BIN_PREFIX)
	$(CXX) $^ $(CXXFLAGS) -o $@
$(FUZZION2HTML_NAME): $(FUZZION2HTML_BIN)

$(FUZZORT_BIN): $(FUZZORT_OBJS) | $(BIN_PREFIX)
	$(CXX) $^ $(CXXFLAGS) -o $@
$(FUZZORT_NAME): $(FUZZORT_BIN)

$(FUZZUM_BIN): $(FUZZUM_OBJS) | $(BIN_PREFIX)
	$(CXX) $^ $(CXXFLAGS) -o $@
$(FUZZUM_NAME): $(FUZZUM_BIN)

$(KMERANK_BIN): $(KMERANK_OBJS) | $(BIN_PREFIX)
	$(CXX) $^ $(CXXFLAGS) -o $@
$(KMERANK_NAME): $(KMERANK_BIN)

$(OBJ_PREFIX)/%.o: $(SRC_PREFIX)/%.cpp | $(OBJ_PREFIX)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_PREFIX) $(OBJ_PREFIX):
	mkdir -p $@

clean:
	rm -r $(BUILD_PREFIX)
