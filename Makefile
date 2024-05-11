## CXX = g++
CXX = clang++
CXXFLAGS = -std=c++17 -Wall -Wextra -Werror
CXXFLAGS += -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-unused-private-field
CXXFLAGS += -O3 -march=native -mtune=native -ffast-math -funroll-loops -fomit-frame-pointer -fno-exceptions -fno-rtti
SRCS = suffix_array/engine.cpp suffix_array/main.cpp
INCLUDES = -I./suffix_array
## OBJS = $(SRCS:.cpp=.o)
LIBS = -lc++ -lc++abi -L/usr/local/lib -lm -L./suffix_array
TARGET = bin/release/suffix_array

ifeq ($(shell uname),Darwin)
	LDFLAGS += -fopenmp
else
	LDFLAGS += -fopenmp
endif

ifeq ($(CXX),clang++)
	CXXFLAGS += -stdlib=libc++
endif

all: install clean

run:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $(SRCS) $(LIBS) -o $(TARGET)
	./$(TARGET)

profile:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $(TARGET) $(SRCS) $(LIBS) -pg
	./$(TARGET)
	gprof $(TARGET) gmon.out > analysis.txt

install:
	python -m pip install .

clean:
	rm -r build dist *.egg-info .cache

get_data:
	bash get_data.sh
