cmake_minimum_required(VERSION 3.10)
project(util_exec C)

# -----------------------------------------------------------------------------
# 1. Compiler and C Standard
# -----------------------------------------------------------------------------
set(CMAKE_C_STANDARD 11)

# -----------------------------------------------------------------------------
# 2. Release Build Flags
# -----------------------------------------------------------------------------
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_FLAGS_RELEASE
    "-Wall -Wextra -O3 -march=native -DNDEBUG -ffast-math -mtune=native -Xpreprocessor -fopenmp"
)

# -----------------------------------------------------------------------------
# 3. Include Directories
# -----------------------------------------------------------------------------
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}          # For project headers (e.g., "globals.h")
    ${CMAKE_CURRENT_SOURCE_DIR}/utils    # For utility headers
    /usr/local/opt/openblas/include      # OpenBLAS
    /usr/local/Cellar/gsl/2.8/include    # GSL
    /usr/local/opt/libomp/include        # OpenMP
    /usr/local/opt/llvm/include          # LLVM
    /usr/local/Cellar/cjson/1.7.18/include # cJSON
)

# -----------------------------------------------------------------------------
# 4. Source Files
# -----------------------------------------------------------------------------
set(SRC_FILES
    main.c
    exact.c
    hitAndRun.c
    multinomial.c
    multivariate-cdf.c
    multivariate-pdf.c
    utils/fileutils.c
    utils/instanceGenerator.c
    utils/matrixUtils.c
    utils/memoizationUtil.c
    utils/multivariateUtils.c
    utils/combinations.c
)

# -----------------------------------------------------------------------------
# 5. Create the Executable
# -----------------------------------------------------------------------------
add_executable(util_exec ${SRC_FILES})

# -----------------------------------------------------------------------------
# 6. Linker Directories and Libraries
# -----------------------------------------------------------------------------
link_directories(
    /usr/local/opt/openblas/lib          # OpenBLAS
    /usr/local/Cellar/gsl/2.8/lib        # GSL
    /usr/local/opt/libomp/lib            # OpenMP
    /usr/local/Cellar/cjson/1.7.18/lib   # cJSON
)

# Link the necessary libraries
target_link_libraries(util_exec
    PUBLIC
	/usr/local/opt/openblas/lib/libopenblas.dylib
    /usr/local/Cellar/gsl/2.8/lib/libgsl.dylib     # GSL main library
    /usr/local/Cellar/gsl/2.8/lib/libgslcblas.dylib # GSL CBLAS backend
    /usr/local/opt/libomp/lib/libomp.dylib         # OpenMP
    /usr/local/Cellar/cjson/1.7.18/lib/libcjson.dylib # cJSON
	m                  # Math library
)

# -----------------------------------------------------------------------------
# 7. RPATH Configuration for macOS
# -----------------------------------------------------------------------------
set(CMAKE_BUILD_RPATH
    "/usr/local/opt/openblas/lib;/usr/local/Cellar/gsl/2.8/lib;/usr/local/opt/libomp/lib;/usr/local/Cellar/cjson/1.7.18/lib"
)
set(CMAKE_INSTALL_RPATH "${CMAKE_BUILD_RPATH}")

# -----------------------------------------------------------------------------
# Done!
# -----------------------------------------------------------------------------
message(STATUS "Configured util_exec with Release build options.")
