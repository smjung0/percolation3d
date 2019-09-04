# Percolation of CNT/polymer composites

Source codes for the three-dimensional simulation of an electrical percolation of a carbon nanotube based polymer composite. The Electromechanical properties of a given CNT/polymer composite under a strain can be simulated by modifying the "input.txt" file for the given simulation condition.



## Preparation of Simulation

### Before compilation 
Sources require Eigen Library as a matrix operation. Please ensure that Eigen library is configured to your environment. The library can be downloaded from the below link.

[Eigen](http://eigen.tuxfamily.org)

### Compilation
The execution file is named to "perc3d" in the given Makefile. You can change the execution file name. "perc3d" will be used in this document.  

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
DIMENSION x(-10.0 [um], +10.0 [um]) y(-10.0 [um], +10.0 [um]) z(-10.0 [um], +10.0 [um]) 
VOLTAGE 1.0 [V]
#LENGTH CONST 5.64
LENGTH WEIBULL 10000 10.5 11
#RADIUS CONST 0.025 
RADIUS LOGNORMAL 10000 1.7883 0.083
TOLERANCE 1.0e-7
CNT noofcnt=1 wall=12 sigma=0.01 [S/um] node=10 tht_max=30.0 [dgr]
POLYMER strain=0.00 mass=1.19e-12 [g/um^3] poisson=0.49 deltaE=1.0 [eV] M=400 Dcutoff=0.0010 [um]
RANDTRIAL n=1 seed=100000
# if seed=-1 : equal spacing
# else if seed=0 : time dependent random
# else random_seed = seed value
CALCULATE VOLF init=4000 noofsteps=1 delta=200
#CALCULATE STRN init=0.00 noofsteps=101 delta=0.01
END
```

- Comment line: started with # in the text file.

  

- To configure shape and boundary condition of calculation domain.

```
STRUCTURE [shape] [boundary_condition] 
```

Box structure: shape=BOX, Cylindrical structure: shape=CYL, Periodic boundary condition: boundary_condition=PBC, Free boundary condition: boundary_condition=FBC



- To set the range of calculation domain (for both BOX and CYL), just change the numbers in the following line.

```
DIMENSION x(-10.0 [um], +10.0 [um]) y(-10.0 [um], +10.0 [um]) z(-10.0 [um], +10.0 [um]) 
```



- To set the voltage difference between two electrodes, just change the number in the following line.
```
VOLTAGE 1.0 [V]
```



- To set the length of CNT, we have two options of constant length and following the Weibull distribution. If you want to use constant length, just use the following line. 

```
LENGTH CONST 5.64
```

Here, 5.64 is an example of the constant CNT length in micrometer. Or, if you want to use Weibull distribution, just use the following line by changing appropriate number for a given simulation condition.

```
LENGTH WEIBULL [rand_seed] [lambda] [k]
```

Here, [rand_seed] is the random seed number for Weibull distribution, [lambda] and [k] are parameters for the Weibull distribution. The relationships between the average length, m_l and standard deviation of the length of CNT, s_l are as follow.

<img src="https://latex.codecogs.com/svg.latex?m_l=\lambda\cdot\Gamma(1+{k}^{-1})" />

$$
m_l = \lambda\cdot  \Gamma (1+{k}^{-1}) \\
s_l = \lambda^2 \cdot \left[\Gamma \left(1+\frac{2}{k} \right) - \left(\Gamma \left(1+\dfrac{1}{k} \right) \right)^2\right]
$$


- To set the radius of CNT, we also have two options of constant radius and following the Log-normal distribution. If you want to use constant radius, just use the following line. 

```
RADIUS CONST 0.025
```

Here, 0.025 is an example of the constant CNT radius in micrometer. Or, if you want to use Log-normal distribution, just use the following line by changing appropriate number for a given simulation condition.

```
RADIUS LOGNORMAL [rand_seed] [mu] [sigma]
```

Here, [rand_seed] is the random seed number for Log-normal distribution, [mu] and [sigma] are parameters for the Log-normal distribution. The relationships between the average radius, m_r and standard deviation of the radius of CNT, s_r are as follow.

$$
m_r = exp \left( \mu + \frac{\sigma^2}{2}  \right) \\
s_r = \left\{ \left[ exp(\sigma^2) - 1 \right] exp(2 \mu + \sigma^2)\right\}^{1/2}   
$$


- To set tolerance of matrix solver, just change the number in the following line.

```
TOLERANCE 1.0e-7
```



- To set the material parameters of CNTs, just change the number in the following line.

```
CNT noofcnt=1 wall=12 sigma=0.01 [S/um] node=10 tht_max=30.0 [dgr]
```

Here, noofcnt is total number of CNTs in the composite. wall is the number of wall of given CNTs. sigma is a conductivity of CNTs, node is the number of initial node representing the curved CNTs. Finally, tht_max is the maximum deviation angle for representing the waviness of the curved CNTs.



- To set the material parameters of polymer and other properties, just change the number in the following line.

```
POLYMER strain=0.00 mass=1.19e-12 [g/um^3] poisson=0.49 deltaE=1.0 [eV] M=400 Dcutoff=0.0010 [um]
```

Here, strain is amount of deformation from initial size. 0 means no strain. mass is the mass density of the polymer, poisson is the Poisson's ratio of the given polymer. deltaE is energy barrier between CNTs and polymer and M is the number of conduction channel through tunneling path. Dcutoff is the cutoff distance representing a criteria for a union operation between two adjacent CNTs.



- To set the number of random trial and its random seed, change the following line.
```
RANDTRIAL n=1 seed=100000
```
Here, n is the number of random trial for single calculation step of volume fraction or strain. seed is the random seed for randomizing the position and shape of CNTs in the calculation domain. If seed<0, random seeds will be equally spaced from 1 to RAND_MAX for the number of random trial, n. Else if seed=0, time dependent random seed will be given, else (seed > 0) seed value will be given to random_seed.



- To calculate the conductance and conductivity as a function of volume fraction, use the following line.

```
CALCULATE VOLF init=4000 noofsteps=1 delta=200
```

Here, keyword VOLF is used for the calculation by volume fraction. init is the initial number of CNTs, noofstep is the number of steps for the calculation, and delta is the difference of the number of CNTs between each calculation step.



- To calculate the conductance and conductivity as a function of strain, use the following line.

```
CALCULATE STRN init=0.00 noofsteps=101 delta=0.01
```

Here, keyword STRN is used for the calculation by volume fraction. init is the initial strain of a composite, noofstep is the number of steps for the calculation, and delta is the difference of strain between each calculation step.



- Use "END" keyword to finish the script.























