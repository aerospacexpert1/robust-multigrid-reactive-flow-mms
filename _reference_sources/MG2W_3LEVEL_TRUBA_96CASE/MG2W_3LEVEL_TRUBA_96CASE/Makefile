CC ?= gcc
CFLAGS ?= -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp
LDLIBS ?= -lm

TARGET := build/opposedflow_mg2w
SOURCE := src/opposedflow_mg2w.c

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(SOURCE)
	@mkdir -p build
	$(CC) $(CFLAGS) $(SOURCE) $(LDLIBS) -o $(TARGET)

clean:
	rm -rf build
