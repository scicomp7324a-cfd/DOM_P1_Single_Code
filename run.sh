#!/bin/bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release &&
cmake --build build &&
cd build
./p1_app