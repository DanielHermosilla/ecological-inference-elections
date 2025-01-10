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

# Targets
BINDIR      = bin
LIB_TARGET  = libutil.dylib    # (on Linux, could be libutil.so)
EXEC_TARGET = util_exec

# List *all* .c files, no matter if top-level or in 'utils/' 
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

# Convert each .c to a .o in bin/ mirrored subdirs
OBJS = $(SRCS:%.c=$(BINDIR)/%.o)

# Dependency files
DEPS = $(OBJS:.o=.d)

# Default rule: build both library + executable
all: $(BINDIR)/$(LIB_TARGET) $(BINDIR)/$(EXEC_TARGET)

# Shared Library
$(BINDIR)/$(LIB_TARGET): $(OBJS)
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS_LIB)

# Executable
$(BINDIR)/$(EXEC_TARGET): $(OBJS)
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS_EXEC)

# Pattern rule: build bin/something.o from something.c 
# This also handles subdirs: e.g. utils/matrixUtils.c -> bin/utils/matrixUtils.o
$(BINDIR)/%.o: %.c
	@mkdir -p $(dir $@)   # create bin/utils/ if needed
	$(CC) $(CFLAGS) -MMD -c $< -o $@

# Include dependency files
-include $(DEPS)

clean:
	rm -f $(BINDIR)/*.o \
	      $(BINDIR)/*/*.o \
	      $(BINDIR)/*.d \
	      $(BINDIR)/*/*.d \
	      $(BINDIR)/$(LIB_TARGET) \
	      $(BINDIR)/$(EXEC_TARGET)
