CXXFLAGS=-std=c++11 $(shell root-config --cflags)
LIBS=$(shell root-config --libs)

run: potSum
	@echo "Complied potSum..."
potSum: potSum.o
	g++ -o $@ $^ ${LIBS}

run: out_inspect
	@echo "Compiled out_inspect..."
out_inspect: out_inspect.o
	g++ -o $@ $^ ${LIBS}

%.o: %.cxx
	g++ ${CXXFLAGS} -o $@ -c $^
