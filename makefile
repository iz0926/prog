CXX = g++
CXXFLAGS = -std=c++11 -O3

SOURCES = randmst.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = randmst

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
