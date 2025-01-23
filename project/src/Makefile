.PHONY: all clean debug release

# Compiler
CC = clang

# Base Compiler Flags
CFLAGS_BASE = -Wall -Wextra -I$(shell R RHOME)/include \
    -I/usr/local/opt/openblas/include \
    -I/usr/local/Cellar/gsl/2.8/include \
    -I/usr/local/opt/libomp/include \
    -Xpreprocessor -fopenmp \
    -I/usr/local/opt/llvm/include -Iutils \
    -I/usr/local/Cellar/cjson/1.7.18/include -I.

# Linker Flags
LDFLAGS_LIB_BASE = -L$(shell R RHOME)/lib \
    -L/usr/local/opt/openblas/lib \
    -L/usr/local/Cellar/gsl/2.8/lib \
    -L/usr/local/opt/libomp/lib \
    -fPIC \
    -lR -lopenblas -lgsl -lgslcblas -lomp \
    -Wl,-rpath,/usr/local/opt/libomp/lib \
    -L/usr/local/Cellar/cjson/1.7.18/lib \
    -lcjson

LDFLAGS_EXEC_BASE = -L$(shell R RHOME)/lib \
    -L/usr/local/opt/openblas/lib \
    -L/usr/local/Cellar/gsl/2.8/lib \
    -L/usr/local/opt/libomp/lib \
    -lR -lopenblas -lgsl -lgslcblas -lomp \
    -Wl,-rpath,/usr/local/opt/libomp/lib \
    -L/usr/local/Cellar/cjson/1.7.18/lib \
    -lcjson

# Debug-specific Flags
CFLAGS_DEBUG = -g -fsanitize=address -fsanitize=leak
LDFLAGS_DEBUG = $(LDFLAGS_LIB_BASE)

# Release-specific Flags
CFLAGS_RELEASE = -O3 -march=native -DNDEBUG -ffast-math -mtune=native
LDFLAGS_RELEASE = $(LDFLAGS_LIB_BASE)

# Directory for build artifacts
BUILDDIR = build

# Targets
LIB_TARGET  = libutil.dylib    # (on Linux, could be libutil.so)
EXEC_TARGET = util_exec

# List *all* .c files
SRCS = \
  main.c \
  exact.c \
  multivariate-pdf.c \
  multinomial.c \
  multivariate-cdf.c \
  hitAndRun.c \
  utils/matrixUtils.c \
  utils/instanceGenerator.c \
  utils/fileutils.c \
  utils/memoizationUtil.c \
  utils/multivariateUtils.c \
  utils/combinations.c \

# Convert each .c to a .o in build/ mirrored subdirs
OBJS = $(SRCS:%.c=$(BUILDDIR)/%.o)

# Dependency files
DEPS = $(OBJS:.o=.d)

# Default rule: build both debug and release
all: debug release

# Debug build
debug: CFLAGS = $(CFLAGS_BASE) $(CFLAGS_DEBUG)
debug: LDFLAGS = $(LDFLAGS_DEBUG)
debug: $(BUILDDIR)/$(LIB_TARGET) $(BUILDDIR)/$(EXEC_TARGET)

# Release build
release: CFLAGS = $(CFLAGS_BASE) $(CFLAGS_RELEASE)
release: LDFLAGS = $(LDFLAGS_RELEASE)
release: $(BUILDDIR)/$(LIB_TARGET) $(BUILDDIR)/$(EXEC_TARGET)

# Shared Library
$(BUILDDIR)/$(LIB_TARGET): $(OBJS)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# Executable
$(BUILDDIR)/$(EXEC_TARGET): $(OBJS)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# Pattern rule: build build/something.o from something.c 
# This also handles subdirs: e.g., utils/matrixUtils.c -> build/utils/matrixUtils.o
$(BUILDDIR)/%.o: %.c
	@mkdir -p $(dir $@)   # create build/utils/ if needed
	$(CC) $(CFLAGS) -MMD -c $< -o $@

# Include dependency files
-include $(DEPS)

clean:
	rm -rf $(BUILDDIR)
