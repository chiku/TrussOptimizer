CXX = g++
CC  = gcc

CCFLAGS += -O2 -Wall
CXXFLAGS += -O2 -Wall

SOURCE_INCLUDE_PATHS += -Ivendor/truss/include -Ivendor/truss/vendor

.PHONY: all
all: main

truss.o:
	@cd vendor/truss && make build/truss.o
	@cp -pv vendor/truss/build/truss.o .

main.o: main.cpp de.h
	${CXX} -c ${CXXFLAGS} ${SOURCE_INCLUDE_PATHS} -o $@ $<

de.o: de.cpp de.h
	${CXX} -c ${CXXFLAGS} ${SOURCE_INCLUDE_PATHS} -o $@ $<

main: main.o de.o truss.o
	${CXX} -o $@ $^ ${LDFLAGS}

.PHONY: clean
clean:
	rm -rf main *.o vendor/truss/build/truss.o
