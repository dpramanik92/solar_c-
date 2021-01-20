#!/bin/bash

g++ main.cpp KamLAND_anti.cpp probability.cpp numerical.cpp read_files.cpp -L/usr/local/lib/ -std=c++14
./a.out
