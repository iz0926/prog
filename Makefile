CXX = g++
CXXFLAGS = -std=c++11

SOURCES = strassen.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = strassen

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)