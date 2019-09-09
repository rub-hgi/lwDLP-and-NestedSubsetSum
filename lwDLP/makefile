.PHONY: all clean
CC = /home/wiemer_sha/local/bin/x86_64-unknown-linux-gnu-gcc

CFLAGS = -Ofast -g -Wextra -Wall -std=c++11 -march=native
LFLAGS = -lm -lstdc++

all: RepsZ


RepsZ: main.cpp main.h
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)


clean:
	rm -f RepsZ