# Compiler and Flags
CC = clang

# Change -O1 to -O3 for efficiency
# Compiler Flags
CFLAGS = -Wall -Wextra -g -fsanitize=address -fsanitize=leak -O3 \
    -I$(shell R RHOME)/include \
    -I/usr/local/opt/openblas/include \
    -I/usr/local/Cellar/gsl/2.8/include \
    -I/usr/local/opt/libomp/include \
    -Xpreprocessor -fopenmp \
    -I/usr/local/opt/llvm/include

# Linker Flags for Shared Library
LDFLAGS_LIB = -L$(shell R RHOME)/lib \
    -L/usr/local/opt/openblas/lib \
    -L/usr/local/Cellar/gsl/2.8/lib \
    -L/usr/local/opt/libomp/lib \
    -shared -dynamiclib -fPIC \
    -lR -lopenblas -lgsl -lgslcblas -lomp \
    -Wl,-rpath,/usr/local/opt/libomp/lib

# Linker Flags for Executable
LDFLAGS_EXEC = -L$(shell R RHOME)/lib \
    -L/usr/local/opt/openblas/lib \
    -L/usr/local/Cellar/gsl/2.8/lib \
    -L/usr/local/opt/libomp/lib \
    -lR -lopenblas -lgsl -lgslcblas -lomp \
    -Wl,-rpath,/usr/local/opt/libomp/lib

# Files and Targets
LIB_TARGET = libutil.dylib  # Use libutil.so for Linux
EXEC_TARGET = util_exec     # Name for the executable
SRCS = util.c matrixUtils.c instanceGenerator.c multinomial.c fileutils.c exact.c memoizationUtil.c
OBJS = $(SRCS:.c=.o)
HEADERS = matrixUtils.h instanceGenerator.h multinomial.h fileutils.h memoizationUtil.h exact.h globals.h uthash.h

# Default Rule
all: $(LIB_TARGET) $(EXEC_TARGET)

# Build Shared Library
$(LIB_TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(LIB_TARGET) $(OBJS) $(LDFLAGS_LIB)

# Build Executable
$(EXEC_TARGET): util.o matrixUtils.o instanceGenerator.o multinomial.o fileutils.o exact.o memoizationUtil.o
	$(CC) $(CFLAGS) -o $(EXEC_TARGET) util.o matrixUtils.o instanceGenerator.o multinomial.o fileutils.o exact.o memoizationUtil.o $(LDFLAGS_EXEC)

# Build Object Files
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean Build Artifacts
clean:
	rm -f $(LIB_TARGET) $(EXEC_TARGET) $(OBJS)

