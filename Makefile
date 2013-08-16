CXX = g++
CC  = gcc

CCFLAGS += -O2 -Wall
CXXFLAGS += -O2 -Wall

SOURCE_INCLUDE_PATHS += -Ivendor/cmatrix

.PHONY: all
all: main

main: main.cpp de.cpp matrix.cpp truss.cpp
	${CXX} -o $@ $< ${SOURCE_INCLUDE_PATHS} ${CXXFLAGS} ${LDFLAGS}
	# ${CXX} -o $@ $^ ${SOURCE_INCLUDE_PATHS} ${CXXFLAGS} ${LDFLAGS}

.PHONE: clean
clean:
	rm -rf main *.o
