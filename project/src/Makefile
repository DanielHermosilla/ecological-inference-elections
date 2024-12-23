# Compilador y flags
CC = clang
CFLAGS = -Wall -Wextra -g -I$(shell R RHOME)/include -I/usr/local/Cellar/gsl/2.8/include
LDFLAGS = -L$(shell R RHOME)/lib -L/usr/local/Cellar/gsl/2.8/lib -shared -fPIC -lR -lgsl -lgslcblas

# Biblioteca compartida
TARGET = util.so
SRCS = util.c

# Regla principal
all: $(TARGET)

# Compilar la biblioteca compartida
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

# Limpieza
clean:
	rm -f $(TARGET)

