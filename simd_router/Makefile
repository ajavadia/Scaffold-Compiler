# the compiler: gcc for C program, define as g++ for C++
CXX=clang++
RM=rm -f

# compiler flags:
CPPFLAGS=-std=c++11 -O3

# define any directories containing header files other than /usr/include
INCLUDES=-I/usr/include

# define library paths in addition to /usr/lib
LFLAGS=-L/usr/lib

# define any libraries to link into executable:
LIBS=-lboost_serialization

# the build target executable:
TARGET=router

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(LFLAGS) $(LIBS) -o $(TARGET) $(TARGET).cpp

clean:
	$(RM) $(TARGET)
