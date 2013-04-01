## Fetch source code

```shell
git clone https://github.com/chiku/TrussSolver.git
cd TrussSolver
git submodule update --init
```

## Compile

```shell 
g++ main.cpp -o main
```

## Execute

```shell
./main
```

## Details


The following files are needed for *truss optimizations*

|File|Purpose|
|----|-------|
|matrix.cpp|Implements matrix manipulations|
|truss.cpp|Truss solver|
|de.cpp|Differential Evolution|
|main.cpp|The main program|
|main|The main executable|

On executing main, these two files are searched for

|Name|Purpose|
|----|-------|
|de.dat|Stores some DE parameters|
|truss.dat|Stores the truss structure|


### Format for de.dat

```text
<lower limit for members> <upper limits for members>
<permissible stress> <permissible displacements in both x, y directions>
<constant penalty> 
<Population size>
<lower limit for F> <upper limit for F>
<lower limit for CR> <upper limit for CR>
```

[A sample **de.dat** is available here](https://raw.github.com/chiku/TrussSolver/master/de.dat).

### Format for truss.dat

```text
<total nodes>
{for each of the nodes the following are to be given}
	<x cor> <ycor>  <is force in x-dir known? (y/n)> <value if answer is y>
					<is force in y-dir known? (y/n)> <value if answer is y>
					<is disp. in x-dir known? (y/n)> <value if answer is y>
					<is disp. in y-dir known? (y/n)> <value if answer is y>
{for each node pair (i, j) j>i the following are to be given}
	<is node i connected to node j? (y/n)> <area if answer is y>
<Young's modulus>
```

[A sample **truss.dat** is available here](https://raw.github.com/chiku/TrussSolver/master/truss.dat).

## Notes

* Each node must either have force B.C.s or displacement B.C.s (B.C. stands for boundary condition)
* The members are numbered on the basis of nodes
   (e.g. member 1 connects nodes 1, 2 & member 2 connects nodes 1 & 3, & so on)
* Since DE is an evolution based algorithm, each run is likely to give different answers.
