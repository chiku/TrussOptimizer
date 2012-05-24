To compile the program
```shell
g++ main.cpp
```

To run the program
```shell
./a.out
```

The following files are needed for TRUSS OPTIMIZATION

<table>
  <tr>
    <th>File</th>
    <th>Purpose</th>
  </tr>
  <tr>
    <td>matrix.cpp</td>
    <td>Implements matrix manipulations</td>
  </tr>
  <tr>
    <td>truss.cpp</td>
    <td>Truss solver</td>
  </tr>
  <tr>
    <td>de.cpp</td>
    <td>Differential Evolution</td>
  </tr>
  <tr>
    <td>main.cpp</td>
    <td>The main program</td>
  </tr>
  <tr>
    <td>main.exe</td>
    <td>The main EXECUTABLE file</td>
  </tr>
</table>

On executing main.exe, the two files are searched for

<table>
  <tr>
    <td>de.dat</td>
    <td>Stores some DE parameters</td>
  </tr>
  <tr>
    <td>truss.dat</td>
    <td>Stores the truss structure</td>
  </tr>
</table>

Format for de.dat

```text
<lower limit for members> <upper limits for members>
<permissible stress> <permissible displacements in both x, y directions>
<constant penalty> 
<Population size>
<lower limit for F> <upper limit for F>
<lower limit for CR> <upper limit for CR>
```

Format for truss.dat

```text
<total nodes>
{for each of the nodes the following are to be given}
	<x cor> <ycor>  <is force in x-dir known(y/n)> <value if y>
					<is force in y-dir known(y/n)> <value if y>
					<is disp. in x-dir known(y/n)> <value if y>
					<is disp. in y-dir known(y/n)> <value if y>
{for each node pair (i, j) j>i the following are to be given}
	<is node i connected to node j(y/n)> <area if y>
<Youngs modulus>
```

Notes
--------

* Each node must either have force B.C.s or displacement B.C.s
* The members are numbered on the basis of nodes
   (e.g. member 1 connects nodes 1, 2 & member 2 connects nodes 1 & 3, & so on..)
* Since DE is an evolution based algo, each run is likely to give different answers.