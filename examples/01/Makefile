SHELL      = /bin/bash
CXX        = g++-10
CXXFLAGS   = -Wall -Wextra -pedantic -O3 -std=c++20 -fconcepts -pthread -I../../
TARGET     = example
SOURCES    = example.cc

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

run:
	./$(TARGET)

.PHONY: clean
clean:
	rm -f *~ $(TARGET)
