.PHONY: all clean debug release bundle_libraries test

# Compiler
CC = clang

# Base Compiler Flags
CFLAGS_BASE = -Wall -Wextra \
    -I/usr/local/opt/openblas/include \
    -I/usr/local/Cellar/gsl/2.8/include \
    -I/usr/local/opt/libomp/include \
    -Xpreprocessor -fopenmp \
    -I/usr/local/opt/llvm/include -Iutils \
    -I/usr/local/Cellar/cjson/1.7.18/include -I.

# Linker Flags
LDFLAGS_LIB_BASE = -L/usr/local/opt/openblas/lib \
    -L/usr/local/Cellar/gsl/2.8/lib \
    -L/usr/local/opt/libomp/lib \
    -fPIC \
    -lopenblas -lgsl -lgslcblas -lomp \
    -Wl,-rpath,@executable_path/libs \
    -L/usr/local/Cellar/cjson/1.7.18/lib \
    -lcjson

LDFLAGS_EXEC_BASE = $(LDFLAGS_LIB_BASE)

# Debug-specific Flags
CFLAGS_DEBUG = -O1 -g -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer
LDFLAGS_DEBUG = $(LDFLAGS_EXEC_BASE)

# Release-specific Flags
CFLAGS_RELEASE = -O3 -march=native -DNDEBUG -ffast-math -mtune=native
LDFLAGS_RELEASE = $(LDFLAGS_EXEC_BASE)

# Directory for build artifacts
BUILDDIR = build
LIBSDIR = $(BUILDDIR)/libs

# Targets
LIB_TARGET  = libutil.dylib
EXEC_TARGET = util_exec
TEST_EXEC = $(BUILDDIR)/test-matrixUtils

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
all: debug release bundle_libraries

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

# Test Target for matrixUtils
test: CFLAGS = $(CFLAGS_BASE) $(CFLAGS_DEBUG)
test: LDFLAGS = $(LDFLAGS_DEBUG)
test: $(TEST_EXEC)
	@echo "Running matrix tests..."
	$(TEST_EXEC)

$(TEST_EXEC): tests/test-matrixUtils.c $(BUILDDIR)/utils/matrixUtils.o
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Pattern rule: build build/something.o from something.c
# This also handles subdirs: e.g., utils/matrixUtils.c -> build/utils/matrixUtils.o
$(BUILDDIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -c $< -o $@

# Include dependency files
-include $(DEPS)

# Clean build artifacts
clean:
	rm -rf $(BUILDDIR)

# Copy libraries into build/libs
bundle_libraries:
	@mkdir -p $(LIBSDIR)
	@echo "Bundling required libraries..."
	@otool -L $(BUILDDIR)/$(EXEC_TARGET) | awk '/\/usr\/local/ { print $$1 }' | while read lib; do \
	    cp "$$lib" $(LIBSDIR); \
	done
	@echo "Libraries bundled into $(LIBSDIR)"
