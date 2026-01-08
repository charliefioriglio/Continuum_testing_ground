CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp

SRC_DIR = src/cxx
SOURCES = src/cxx/tools.cpp src/cxx/gauss.cpp src/cxx/shell.cpp src/cxx/molecule.cpp src/cxx/dyson.cpp src/cxx/grid.cpp src/cxx/math_special.cpp src/cxx/cross_section.cpp src/cxx/angle_grid.cpp src/cxx/rotation.cpp src/cxx/continuum.cpp src/cxx/beta.cpp src/cxx/num_eikr.cpp src/cxx/beta_gen.cpp src/cxx/dyson_gen.cpp src/cxx/clebsch_gordan.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLES = dyson_gen beta_gen

all: $(EXECUTABLES)

dyson_gen: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o dyson_gen $(filter-out src/cxx/beta_gen.o, $(OBJECTS))

beta_gen: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o beta_gen $(filter-out src/cxx/dyson_gen.o, $(OBJECTS))

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRC_DIR)/*.o $(EXECUTABLE)
