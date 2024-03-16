CXX = clang++
CXXFLAGS = -std=c++17 -Wall -Wextra -Werror -stdlib=libc++ 
CXXFLAGS += -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-unused-private-field
CXXFLAGS += -O3 -march=native -mtune=native -ffast-math
SRCS = $(wildcard suffix_array/*.cpp)
INCLUDES = -I./suffix_array -I./suffix_array/libsais/include
OBJS = $(SRCS:.cpp=.o)
LIBS = -lc++ -lc++abi -L/usr/local/lib -lm -L./suffix_array -llibsais
TARGET = suffix_array_exe

ifeq ($(shell uname),Darwin)
	LDFLAGS += -fopenmp
else
	LDFLAGS += -fopenmp
endif

all: install clean

run:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $(TARGET) $(SRCS) $(LIBS) 
	./$(TARGET)

install:
	python -m pip install .

clean:
	rm -r build dist *.egg-info .cache
