CXX = g++
CC  = gcc

CCFLAGS += -O2 -Wall
CXXFLAGS += -O2 -Wall

SOURCE_INCLUDE_PATHS += -Ivendor/cmatrix

.PHONY: all
all: main

truss.o:
	@cd vendor/truss && make
	cp vendor/truss/build/truss.o .

main: main.cpp de.cpp truss.o
	${CXX} -o $@ $< truss.o ${SOURCE_INCLUDE_PATHS} ${CXXFLAGS} ${LDFLAGS}
	# ${CXX} -o $@ $^ ${SOURCE_INCLUDE_PATHS} ${CXXFLAGS} ${LDFLAGS}

.PHONE: clean
clean:
	rm -rf main *.o
