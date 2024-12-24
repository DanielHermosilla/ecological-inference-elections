# Compiler and Flags
CC = clang
CFLAGS = -Wall -Wextra -g -I$(shell R RHOME)/include -I/usr/local/opt/openblas/include
LDFLAGS = -L$(shell R RHOME)/lib -L/usr/local/opt/openblas/lib -shared -dynamiclib -fPIC -lR -lopenblas -lgsl -lgslcblas

# Shared Library
TARGET = util.so
SRCS = util.c

# Default Rule
all: $(TARGET)

# Build Shared Library
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

# Clean Build Artifacts
clean:
	rm -f $(TARGET)

