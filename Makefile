CC = gcc
CXX = g++
CXXFLAGS = -g -fPIC -Wall -O2 -std=c++14
LIB = /usr/local/lib
LDFLAGS = -shared

TARGET_LIB = libchi2.so

TARGET = prob_main
TARGET1 = event_main
CHISQ = chi2_main
CHI_PY = chi2_py


TEST= match

SRCS = numerical.cpp probability.cpp event.cpp read_files.cpp interactions.cpp visible_anti.cpp invisible_e.cpp convers_prob.cpp dec_probability.cpp chisqmin.cpp minimizer.cpp

OBJS = numerical.o probability.o event.o read_files.o interactions.o visible_anti.o invisible_e.o convers_prob.o dec_probability.o chisqmin.o minimizer.o

.PHONY: all

all: $(TARGET) $(TARGET1) $(TEST) $(CHISQ) $(CHI_PY) $(TARGET_LIB)

$(TARGET_LIB): $(OBJS) $(CHI_PY).o 
	$(CXX) $(LDFLAGS) -o $@ $^


$(CHI_PY): $(CHI_PY).o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(CHI_PY) $(CHI_PY).o $(OBJS)

$(CHISQ): $(CHISQ).o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(CHISQ) $(CHISQ).o $(OBJS)


$(CHISQ).o: $(CHISQ).cpp $(OBJS)
	$(CXX) $(CXXFLAGS) -c $(CHISQ).cpp -L$(LIB)

$(TEST): $(TEST).o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TEST) $(TEST).o $(OBJS)


$(TEST).o: $(TEST).cpp $(OBJS)
	$(CXX) $(CXXFLAGS) -c $(TEST).cpp -L$(LIB)

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



chisqmin.o: chisqmin.cpp chisqmin.hpp
	$(CXX) $(CXXFLAGS) -c chisqmin.cpp



clear_obj:
	rm -f *.o


clean:
	rm -f *.o *.out $(TARGET) $(TARGET1)
