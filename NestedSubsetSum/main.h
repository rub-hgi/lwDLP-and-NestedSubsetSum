#pragma once
#ifndef MAIN
#define MAIN
#include <iostream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <fstream>
#include <string.h>
#include <random>
#include <stdio.h>
#include <stdint.h>
#include <omp.h>

//size of subset sum instance
int n = 24;
uint64_t modul=2017;

//for flavouring we choose previous_prime(modul);
uint64_t modul_p = 2017;
//48: 12271507
//40:658001
//32: 35951;
//24: 2017;
//16:  113;

int wt = n / 8;
#endif
