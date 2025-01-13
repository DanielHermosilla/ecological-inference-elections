#!/bin/bash

# Enable ASan environment variables only for this execution
export ASAN_OPTIONS=detect_leaks=1:leak_check_at_exit=1
export DYLD_INSERT_LIBRARIES=/usr/local/Cellar/llvm/19.1.6/lib/clang/19/lib/darwin/libclang_rt.asan_osx_dynamic.dylib
export DYLD_FORCE_FLAT_NAMESPACE=1

# Run your C program
./build/util_exec "$@"
