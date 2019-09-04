# Percolation of CNT/polymer composites

Source codes for the three-dimensional simulation of an electrical percolation of a carbon nanotube based polymer composite. The Electromechanical properties of a given CNT/polymer composite under a strain can be simulated by modifying the "input.txt" file for the given simulation condition.

## Preparation of Simulation

### Before compilation 
Sources require Eigen Library as a matrix operation. Please ensure that Eigen library is configured to your environment. The library can be downloaded from the below link.

```
[Eigen](http://eigen.tuxfamily.org)
```

### Compilation
The execution file is named to "perc3d" in the given Makefile. You can change the execution file name. "perc3d" is used in this document  

### After making execution file
The simulation is prepared. Just type the execution file name followed by an input text file as below. 

```
$ perc3d input.txt
```

However you can change file name from input.txt to any given name as you wish. For example, if you want the file name to be "condition_1.txt", just type as below.

```
$ perc3d condition_1.txt
```



## Format of Input File

An example of the input file is as below.

```
# Input code for Percolation Network of CNT/Polymer Composite
STRUCTURE BOX PBC
DIMENSION x(-10.0 [um], +10.0 [um]) y(-10.0 [um], +10.0 [um]) z(-20.0 [um], +20.0 [um]) 
VOLTAGE 1.0 [V]
#LENGTH AVG = 10.0um, STD = 1.22um
LENGTH WEIBULL 10000 10.5 11
#RADIUS AVG = 6nm, STD = 0.5nm
RADIUS LOGNORMAL 10000 1.7883 0.083
TOLERANCE 1.0e-6
CNT noofcnt=1 wall=12 sigma=0.01 [S/um] node=10 tht_max=30.0 [dgr]
POLYMER strain=0.00 mass=1.19e-12 [g/um^3] poisson=0.49 deltaE=1.0 [eV] M=400 Dcutoff=0.0010 [um]
RANDTRIAL n=1 seed=100000
# if seed=-1 : equal spacing
# else if seed=-1 : time dependent random
# else random_seed = seed value
CALCULATE VOLF init=4000 noofsteps=1 delta=200
#CALCULATE STRN init=0.00 noofsteps=101 delta=0.01
END
```





