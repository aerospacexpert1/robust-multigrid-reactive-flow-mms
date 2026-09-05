CC ?= gcc
CFLAGS ?= -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp
LDLIBS ?= -lm

TARGET := build/opposedflow_mg3v
SOURCE := src/opposedflow_mg3v.c

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(SOURCE)
	@mkdir -p build
	$(CC) $(CFLAGS) $(SOURCE) $(LDLIBS) -o $(TARGET)

clean:
	rm -rf build
