# Compilador y flags
CC = clang
CFLAGS = -Wall -Wextra -g -I$(shell R RHOME)/include
LDFLAGS = -L$(shell R RHOME)/lib -shared -fPIC -lR

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

