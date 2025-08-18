SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .cc .o

CONFIG = -Wall -pipe -Wno-psabi
DEBUG = -g -O0
DEFINES = #-DBVH_SAH

OPTIMISE = -O3 -march=native -flto=auto -fopenmp

# default if not passed on command line
BUILD_TYPE ?= RELEASE

# compiler to use
COMPILER ?= gcc

ifeq ($(COMPILER), gcc)
CXX := g++-15
else ifeq ($(COMPILER), icc)
CXX := icpc
else ifeq ($(COMPILER), clang)
CXX := clang++
else
$(error Unrecognised compiler $(COMPILER))
endif

ifeq ($(BUILD_TYPE), DEBUG)
CONFIG += $(DEBUG)
else ifeq ($(BUILD_TYPE), RELEASE)
CONFIG += $(OPTIMISE)
DEFINES += -DNDEBUG
endif

SRCDIR = src
HDRDIR = include
EXTDIR = external
OBJDIR = obj
BINDIR = bin

EXEC = $(BINDIR)/rt_weekend

# Always add extra flags
CXXFLAGS += $(CONFIG)
CXXFLAGS += $(DEFINES)

# Add extra flags from command line
CXXFLAGS += $(EXTRA_CONFIG)
CXXFLAGS += $(EXTRA_DEFINES)

# Flag to emit header dependencies
MDFLAG = -MMD

# Flag to set language standard
STDFLAG = -std=c++17

# Arguments for linking profiler
# PROFFLAGS = -L$(shell brew --prefix gperftools)/lib -lprofiler

# Get the directory structure of the source tree
STRUCTURE := $(shell find $(SRCDIR) -type d)

# Get list of files within the source tree
CODEFILES := $(addsuffix /*,$(STRUCTURE))
CODEFILES := $(wildcard $(CODEFILES))

# Filter into categories
SRCFILES := $(filter %.cpp,$(CODEFILES))
OBJFILES := $(subst $(SRCDIR),$(OBJDIR),$(SRCFILES:%.cpp=%.o))
DEPFILES := $(OBJFILES:%.o=%.d)

.PHONY : build clean rebuild

# Default build target
build: $(EXEC)

# Main executable build recipe
$(EXEC): $(OBJFILES)
	@[ -d $(BINDIR) ] || mkdir $(BINDIR)
	$(CXX) $(CXXFLAGS) $(STDFLAG) $(PROFFLAGS) -o $@ $^

# Recipe for building object file from source
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(STDFLAG) $(PROFFLAGS) $(MDFLAG) -I$(HDRDIR) -I$(EXTDIR) -c $< -o $@

# Remove all object files
clean:
	-rm $(OBJFILES) $(DEPFILES) $(EXEC) $(LIB)

# .PHONY : print-%
# print-% :
# 	@echo $* = $($*)

rebuild: clean build

-include $(OBJDIR)/*.d
