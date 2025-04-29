#!/usr/bin/make -f

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -fopenmp -Iexternal -O3

SRCS = src/cli.cpp src/matrix/csrd_matrix.cpp src/matrix/coo_matrix.cpp src/monte_carlo_ode_solver.cpp src/utils/read_vector.cpp src/spice/spice_parser.cpp src/spice/time_function.cpp
OBJS = $(SRCS:src/%.cpp=bin/%.o)

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
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

bin/%.o: src/%.cpp | $(OBJDIRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIRS):
	mkdir -p $@

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
	rm -f $(TARGET) $(OBJS)
