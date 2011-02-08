TARGET = modelCC83
CXX = gcc
CXXFLAGS = -Wall -lstdc++ -lgsl -lgslcblas
OBJ = Simulation.o Population.o Individual.o Chromosome.o
LINK = -L. -lTEs

all: libTEs.a $(TARGET)

libTEs.a: $(OBJ)
	rm -f $@
	ar rsc $@ $(OBJ)

modelCC83: modelCC83.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LINK)

clean:
	@find . -name '*~' -exec rm {} \;
	@find . -name '*.[oa]' -exec rm {} \;
	@if test -e $(TARGET); then rm -f $(TARGET); fi

test: test.cpp libTEs.a
	@if test -e $@; then rm $@; fi
	$(CXX) $(CXXFLAGS) $< -o $@ $(LINK)
