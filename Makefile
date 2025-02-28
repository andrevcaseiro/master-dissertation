#!/usr/bin/make -f

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -fopenmp

SRCS = src/main.cpp src/matrix/csrd_matrix.cpp src/mcme.cpp
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
