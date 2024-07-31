CC = gcc
CFLAGS = -std=c11
CFLAGS += -O3 -march=native -mtune=native -ffast-math -funroll-loops
CFLAGS += -Wno-unused-result -Wno-unused-variable -Wno-unused-function -Wno-unused-but-set-variable
SRCS   = suffix_array/engine.c suffix_array/main.c suffix_array/libsais.c suffix_array/libsais64.c
INCLUDES = -I./suffix_array
LIBS = -lc++ -lc++abi -L/usr/local/lib -lm -L./suffix_array
TARGET = bin/release/suffix_array

ifeq ($(shell uname),Darwin)
	LDFLAGS += -fopenmp
else
	LDFLAGS += -fopenmp
endif

all: install clean

run:
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) $(SRCS) $(LIBS) -o $(TARGET)
	./$(TARGET)

profile:
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $(TARGET) $(SRCS) $(LIBS) -pg
	./$(TARGET)
	gprof $(TARGET) gmon.out > analysis.txt

install:
	python -m pip install .

clean:
	rm -r build dist *.egg-info .cache

get_data:
	bash get_data.sh

demo:
	python searchapp_demo/main.py
