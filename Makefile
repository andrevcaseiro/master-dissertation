#!/usr/bin/make -f

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Iexternal/CLI11/include

SRCS = src/main.cpp src/matrix/CSRMatrix.cpp src/MCME.cpp
OBJS = $(SRCS:.cpp=.o)

TARGET = main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJS)
