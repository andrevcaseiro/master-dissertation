#!/usr/bin/make -f

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -fopenmp -Iexternal

SRCS = src/cli.cpp src/matrix/csrd_matrix.cpp src/monte_carlo_ode_solver.cpp src/read_vector.cpp
OBJS = $(SRCS:src/%.cpp=bin/%.o)

# Find all unique directories in src and map them to bin
SRCDIRS = $(sort $(dir $(SRCS)))
OBJDIRS = $(SRCDIRS:src/%=bin/%)

TARGET = main

# Download CLI11
CLI11_HPP = external/CLI11.hpp
CLI11_URL = https://github.com/CLIUtils/CLI11/releases/download/v2.5.0/CLI11.hpp

all: $(TARGET) $(CLI11_HPP)

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

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJS)
