#!/usr/bin/make -f

CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -fopenmp -Isrc -Iexternal -Og -g -march=native
LDFLAGS = -lsuperlu_mt_OPENMP -lblas_OPENMP -lHighFM -lfmt

# HighFM flags
CXXFLAGS += -I/home/andre_caseiro/intel/oneapi/2025.2/include -m64 -DHIGHFM_WITH_MKL -DHIGHFM_WITH_COMPLEX -DMKL_ILP64
LDFLAGS += -L/home/andre_caseiro/intel/oneapi/2025.2/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# Add dependency generation flags
DEPFLAGS = -MMD -MP

# Automatically find all cpp files
SRCS = $(shell find src -name "*.cpp")
OBJS = $(SRCS:src/%.cpp=bin/%.o)
DEPS = $(OBJS:.o=.d)   # Generate dependency files

# Find all unique directories in src and map them to bin
SRCDIRS = $(sort $(dir $(SRCS)))
OBJDIRS = $(SRCDIRS:src/%=bin/%)

TARGET = main

# Download CLI11
CLI11_HPP = external/CLI11.hpp
CLI11_URL = https://github.com/CLIUtils/CLI11/releases/download/v2.5.0/CLI11.hpp

# Download eigen
EIGEN_URL = https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
EIGEN_ARCHIVE = $(notdir $(EIGEN_URL))

.PHONY: all
all: external/Eigen/.dirstamp $(CLI11_HPP) $(TARGET)

cli: src/cli.cpp

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS) $(LDFLAGS)

bin/%.o: src/%.cpp | $(OBJDIRS)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIRS):
	mkdir -p $@

# Include generated dependency files
-include $(DEPS)

$(CLI11_HPP):
	mkdir -p external
	wget -O $(CLI11_HPP) $(CLI11_URL)

external/Eigen/.dirstamp:
	mkdir -p external
	wget -O $(EIGEN_ARCHIVE) $(EIGEN_URL)
	tar --bzip2 -xf $(EIGEN_ARCHIVE) --strip-components=1 -C external eigen-3.4.0/Eigen
	rm $(EIGEN_ARCHIVE)
	touch $@

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJS) $(DEPS)
