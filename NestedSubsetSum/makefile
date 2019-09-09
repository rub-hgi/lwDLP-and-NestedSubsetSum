.PHONY: all clean
CC = g++

CFLAGS = -Ofast -g -Wextra -Wall -std=c++11 -march=native
LFLAGS = -lm -lstdc++

all: NestedRho


NestedRho: main.cpp main.h
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)


clean:
	rm -f NestedRho