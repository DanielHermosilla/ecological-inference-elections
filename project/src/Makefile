# Compiler and Flags
CC = clang

# Compiler Flags
CFLAGS = -Wall -Wextra -g \
    -I$(shell R RHOME)/include \
    -I/usr/local/opt/openblas/include \
    -I/usr/local/Cellar/gsl/2.8/include \
    -I/usr/local/opt/libomp/include \
    -Xpreprocessor -fopenmp \
    -I/usr/local/opt/llvm/include

# Linker Flags
LDFLAGS = -L$(shell R RHOME)/lib \
    -L/usr/local/opt/openblas/lib \
    -L/usr/local/Cellar/gsl/2.8/lib \
    -L/usr/local/opt/libomp/lib \
    -shared -dynamiclib -fPIC \
    -lR -lopenblas -lgsl -lgslcblas -lomp \
    -Wl,-rpath,/usr/local/opt/libomp/lib

# Files and Targets
TARGET = util.so
SRCS = util.c matrixUtils.c
OBJS = $(SRCS:.c=.o)
HEADERS = matrixUtils.h

# Default Rule
all: $(TARGET)

# Build Shared Library
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Build Object Files
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean Build Artifacts
clean:
	rm -f $(TARGET) $(OBJS)

