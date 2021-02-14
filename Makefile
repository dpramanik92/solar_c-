CC = gcc
CXX = g++
CXXFLAGS = -g -Wall -std=c++14
LIB = /usr/local/lib

TARGET = prob_main
TARGET1 = event_main

SRCS = numerical.cpp probability.cpp event.cpp read_files.cpp interactions.cpp visible_anti.cpp invisible_e.cpp convers_prob.cpp dec_probability.cpp

OBJS = numerical.o probability.o event.o read_files.o interactions.o visible_anti.o invisible_e.o convers_prob.o dec_probability.o

all: $(TARGET) $(TARGET1)

$(TARGET): $(TARGET).o $(OBJS)
	$(CXX) $(CXXFLAGS) -o  $(TARGET) $(TARGET).o $(OBJS)

$(TARGET).o: $(TARGET).cpp $(SRCS)
	$(CXX)  $(CXXFLAGS) -c $(TARGET).cpp -L$(LIB)


$(TARGET1): $(TARGET1).o $(OBJS)
	$(CXX) $(CXXFLAGS) -o  $(TARGET1) $(TARGET1).o $(OBJS)

$(TARGET1).o: $(TARGET1).cpp $(SRCS)
	$(CXX)  $(CXXFLAGS) -c $(TARGET1).cpp -L$(LIB)

numerical.o: numerical.cpp numerical.hpp
	$(CXX) $(CXXFLAGS)  -c numerical.cpp -L$(LIB)


probability.o: probability.cpp probability.hpp
	$(CXX) $(CXXFLAGS) -c probability.cpp


read_files.o: read_files.cpp read_files.hpp
	$(CXX) $(CXXFLAGS) -c read_files.cpp


event.o: event.cpp event.hpp
	$(CXX) $(CXXFLAGS) -c event.cpp -L$(LIB)


interactions.o: interactions.cpp interactions.hpp
	$(CXX) $(CXXFLAGS) -c interactions.cpp


visible_anti.o: visible_anti.cpp visible_anti.hpp
	$(CXX) $(CXXFLAGS) -c visible_anti.cpp


invisible_e.o: invisible_e.cpp invisible_e.hpp
	$(CXX) $(CXXFLAGS) -c invisible_e.cpp



dec_probability.o: dec_probability.cpp dec_probability.hpp
	$(CXX) $(CXXFLAGS) -c dec_probability.cpp

convers_prob.o: convers_prob.cpp convers_prob.hpp
	$(CXX) $(CXXFLAGS) -c convers_prob.cpp

clear_obj:
	rm -f *.o


clean:
	rm -f *.o *.out $(TARGET) $(TARGET1)
