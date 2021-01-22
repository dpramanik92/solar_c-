CC = gcc
CXX = g++
CXXFLAGS = -g -Wall -std=c++14
LIB = /usr/local/lib

TARGET = main


SRCS = numerical.cpp probability.cpp KamLAND_anti.cpp read_files.cpp interactions.cpp visible_anti.cpp

OBJS = $(TARGET).o numerical.o probability.o KamLAND_anti.o read_files.o interactions.o visible_anti.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o  $(TARGET) $(OBJS)

$(TARGET).o: $(TARGET).cpp $(SRCS)
	$(CXX)  $(CXXFLAGS) -c $(TARGET).cpp -L$(LIB)


numerical.o: numerical.cpp numerical.hpp
	$(CXX) $(CXXFLAGS)  -c numerical.cpp -L$(LIB)


probability.o: probability.cpp probability.hpp
	$(CXX) $(CXXFLAGS) -c probability.cpp


read_files.o: read_files.cpp read_files.hpp
	$(CXX) $(CXXFLAGS) -c read_files.cpp


KamLAND_anti.o: KamLAND_anti.cpp KamLAND_anti.hpp
	$(CXX) $(CXXFLAGS) -c KamLAND_anti.cpp -L$(LIB)


interactions.o: interactions.cpp interactions.hpp
	$(CXX) $(CXXFLAGS) -c interactions.cpp

visible_anti.o: visible_anti.cpp visible_anti.hpp
	$(CXX) $(CXXFLAGS) -c visible_anti.cpp


clear_obj:
	rm -f *.o


clean:
	rm -f *.o *.out $(TARGET)
