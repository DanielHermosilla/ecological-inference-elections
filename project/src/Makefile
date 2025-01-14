.PHONY: build

# Compiler
CC = clang

# Compiler & Link Flags
CFLAGS = -Wall -Wextra -g -fsanitize=address -fsanitize=leak -O3 \
    -I$(shell R RHOME)/include \
    -I/usr/local/opt/openblas/include \
    -I/usr/local/Cellar/gsl/2.8/include \
    -I/usr/local/opt/libomp/include \
    -Xpreprocessor -fopenmp \
    -I/usr/local/opt/llvm/include -Iutils -I.

LDFLAGS_LIB = -L$(shell R RHOME)/lib \
    -L/usr/local/opt/openblas/lib \
    -L/usr/local/Cellar/gsl/2.8/lib \
    -L/usr/local/opt/libomp/lib \
    -shared -dynamiclib -fPIC \
    -lR -lopenblas -lgsl -lgslcblas -lomp \
    -Wl,-rpath,/usr/local/opt/libomp/lib

LDFLAGS_EXEC = -L$(shell R RHOME)/lib \
    -L/usr/local/opt/openblas/lib \
    -L/usr/local/Cellar/gsl/2.8/lib \
    -L/usr/local/opt/libomp/lib \
    -lR -lopenblas -lgsl -lgslcblas -lomp \
    -Wl,-rpath,/usr/local/opt/libomp/lib

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
  utils/matrixUtils.c \
  utils/instanceGenerator.c \
  utils/fileutils.c \
  utils/memoizationUtil.c \
  utils/multivariateUtils.c

# Convert each .c to a .o in build/ mirrored subdirs
OBJS = $(SRCS:%.c=$(BUILDDIR)/%.o)

# Dependency files
DEPS = $(OBJS:.o=.d)

# Default rule: build both library + executable
all: $(BUILDDIR)/$(LIB_TARGET) $(BUILDDIR)/$(EXEC_TARGET)

# Shared Library
$(BUILDDIR)/$(LIB_TARGET): $(OBJS)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS_LIB)

# Executable
$(BUILDDIR)/$(EXEC_TARGET): $(OBJS)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS_EXEC)

# Pattern rule: build build/something.o from something.c 
# This also handles subdirs: e.g., utils/matrixUtils.c -> build/utils/matrixUtils.o
$(BUILDDIR)/%.o: %.c
	@mkdir -p $(dir $@)   # create build/utils/ if needed
	$(CC) $(CFLAGS) -MMD -c $< -o $@

# Include dependency files
-include $(DEPS)

clean:
	rm -f $(BUILDDIR)/*.o \
	      $(BUILDDIR)/*/*.o \
	      $(BUILDDIR)/*.d \
	      $(BUILDDIR)/*/*.d \
	      $(BUILDDIR)/$(LIB_TARGET) \
	      $(BUILDDIR)/$(EXEC_TARGET)
