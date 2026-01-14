CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp

SRC_DIR = src/cxx
LIB_SOURCES = src/cxx/tools.cpp src/cxx/gauss.cpp src/cxx/shell.cpp src/cxx/molecule.cpp src/cxx/dyson.cpp src/cxx/grid.cpp src/cxx/math_special.cpp src/cxx/cross_section.cpp src/cxx/angle_grid.cpp src/cxx/rotation.cpp src/cxx/continuum.cpp src/cxx/beta.cpp src/cxx/num_eikr.cpp src/cxx/point_dipole.cpp src/cxx/clebsch_gordan.cpp
LIB_OBJECTS = $(LIB_SOURCES:.cpp=.o)

EXECUTABLES = dyson_gen beta_gen dipole_verifier

all: $(EXECUTABLES)

dyson_gen: src/cxx/dyson_gen.o $(LIB_OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

beta_gen: src/cxx/beta_gen.o $(LIB_OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

dipole_verifier: src/cxx/dipole_verifier.o $(LIB_OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRC_DIR)/*.o $(EXECUTABLE)
