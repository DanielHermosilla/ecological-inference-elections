#!/bin/bash

for i in {1..20}; do
    OPENBLAS_NUM_THREADS=1 ./build/util_exec "$i"
done
