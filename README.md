# Percolation of CNT/polymer composites

Source codes for the three-dimensional simulation of an electrical percolation of a carbon nanotube based polymer composite. The Electromechanical properties of a given CNT/polymer composite under a strain can be simulated by modifying the "input.txt" file for the given simulation condition.



## Author

- Sungmin Jung

  Nano-science and Technology Group, Electrical Engineering Division, 

  Department of Engineering, University of Cambridge

  

## Preparation of Simulation

### Before compilation 
Sources require Eigen Library as a matrix operation. Please ensure that Eigen library is configured to your environment. The library can be downloaded from the below link.

[Eigen](http://eigen.tuxfamily.org)

### Compilation
The execution file is named to "perc3d" in the given Makefile. You can change the execution file name. "perc3d" will be used in this document.  



## Simulation

### Execution file

The simulation is prepared. Just type the execution file name followed by an input text file as below. 

```
$ perc3d input.txt
```

However you can change file name from input.txt to any given name as you wish. For example, if you want the file name to be "condition_1.txt", just type as below.

```
$ perc3d condition_1.txt
```



### Format of Input File

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

Here, [rand_seed] is the random seed number for [Weibull distribution](http://www.cplusplus.com/reference/random/), [lambda] and [k] are parameters for the Weibull distribution. The relationships between the average length, $m_l$, standard deviation of the length of CNT, s_l and [lambda], [k] are as follows.

<img src="https://latex.codecogs.com/svg.latex?m_l=\lambda\cdot\Gamma(1+{k}^{-1})" />

<img src="https://latex.codecogs.com/svg.latex?s_l=\lambda^2\cdot\left[\Gamma\left(1+\frac{2}{k}\right)-\left(\Gamma\left(1+\dfrac{1}{k}\right)\right)^{2}\right]" />




- To set the radius of CNT, we also have two options of constant radius and following the Log-normal distribution. If you want to use constant radius, just use the following line. 

```
RADIUS CONST 0.025
```

Here, 0.025 is an example of the constant CNT radius in micrometer. Or, if you want to use Log-normal distribution, just use the following line by changing appropriate number for a given simulation condition.

```
RADIUS LOGNORMAL [rand_seed] [mu] [sigma]
```

Here, [rand_seed] is the random seed number for [Log-normal distribution](http://www.cplusplus.com/reference/random/), [mu] and [sigma] are parameters for the Log-normal distribution. The relationships between the average radius, m_r, standard deviation of the radius of CNT, s_r and [mu] and [sigma] are as follow.

<img src="https://latex.codecogs.com/svg.latex?m_r=exp\left(\mu+\frac{\sigma^2}{2}\right)" />

<img src="https://latex.codecogs.com/svg.latex?s_r=\sqrt{[exp(\sigma^2)-1]exp(2\mu+\sigma^{2})}" />




- To set tolerance of matrix solver, just change the number in the following line.

```
TOLERANCE 1.0e-7
```



- To set the material parameters of CNTs, just change the number in the following line.

```
CNT save=yes noofcnt=1 wall=12 sigma=0.01 [S/um] node=10 tht_max=30.0 [dgr]
```

Here, if save=yes, CNT data will be saved to './cntdata/' folder and if save=no, CNT data will not be saved.  noofcnt is total number of CNTs in the composite. wall is the number of wall of given CNTs. sigma is a conductivity of CNTs, node is the number of initial node representing the curved CNTs. Finally, tht_max is the maximum deviation angle for representing the waviness of the curved CNTs.



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



## Results

### Volume fraction vs. Conductivity(or Conductance) 

If you use the CALCULATE keyword with VOLF as below, The "conductance"_VOLF.txt" and "conductivity_VOLF.txt" files will be generated as results of the simulation. 
```
CALCULATE VOLF init=4000 noofsteps=1 delta=200
```



The file format of the result is as belows.

```
VF[%]		WF[%]		TRY(1)		TRY(2)		TRY(3)		TRY(4)
4.74E+00	7.26E+00	2.92E-04	2.74E-04	2.88E-04	4.69E-04
4.85E+00	7.42E+00	3.19E-04	3.11E-04	3.12E-04	5.51E-04
...
...
5.18E+00	7.92E+00	4.50E-04	4.24E-04	4.40E-04	7.32E-04
```

Here, the number i in TRY(i) is an index of single Monte-Carlo trial and 1 ≤ i ≤ n for the number of random trial. VF and WF are volume fraction and weight fraction of CNT in a composite. Unit of the data are simens [S] for "conductance_VOLF.txt" and simens per meter [S/m] for "conductivity_VOLF.txt". 



### Strain vs. Conductivity(or Conductance) 

If you use the CALCULATE keyword with STRN as below, The "conductance"_STRN.txt" and "conductivity_STRN.txt" files will be generated as results of the simulation. 

```
CALCULATE STRN init=0.00 noofsteps=10 delta=0.02
```



The file format of the result is as belows.

```
STRAIN[%]	VF[vol%]	WF[wt%]		TRY(1)		TRY(2)		TRY(3)		TRY(4)
0.00E+00	3.04E+00	4.24E+00	1.70E-05	1.57E-05	0.00E+00	4.92E-05
2.00E-02	3.04E+00	4.24E+00	1.80E-05	0.00E+00	0.00E+00	4.89E-05
4.00E-02	3.04E+00	4.24E+00	0.00E+00	0.00E+00	0.00E+00	4.85E-05
6.00E-02	3.03E+00	4.24E+00	0.00E+00	0.00E+00	0.00E+00	3.49E-05
8.00E-02	3.03E+00	4.24E+00	0.00E+00	0.00E+00	0.00E+00	3.48E-05
1.00E-01	3.03E+00	4.23E+00	0.00E+00	0.00E+00	0.00E+00	2.58E-05
1.20E-01	3.03E+00	4.23E+00	0.00E+00	0.00E+00	0.00E+00	2.37E-05
1.40E-01	3.03E+00	4.23E+00	1.20E-05	0.00E+00	0.00E+00	2.26E-05
1.60E-01	3.03E+00	4.23E+00	1.25E-05	0.00E+00	0.00E+00	2.25E-05
1.80E-01	3.03E+00	4.23E+00	1.25E-05	0.00E+00	0.00E+00	2.24E-05
2.00E-01	3.03E+00	4.23E+00	1.58E-05	0.00E+00	0.00E+00	2.27E-05
```

Here, the number i in TRY(i) is also an index of single Monte-Carlo trial and 1 ≤ i ≤ n for the number of random trial. Unit of the data are simens [S] for "conductance_VOLF.txt" and simens per meter [S/m] for "conductivity_VOLF.txt".



### CNT data

After finishing the calculation, all the distribution and the information of percolation network for single step  will be saved a folder "./cntdata". The format of file name is as below.


```
cnt_nXXX_sY.YY_rZZZZ.txt
```

Here, XXX is the number of CNTs, Y.YY is the amount of strain and ZZZZ is the random seed for this simulation. Take a look inside the data file as below.

```
DOMAIN: BOX
BNDC: PBC
MIN:(-10.000, -10.000, -10.000), MAX:(10.000, 10.000, 10.000)
No. of CNTs: 5704
Strain: 0.000
Poisson ratio: 0.490
Total No. of Perc. Network is 1
144

#CNT Data 

Address: 0
Parent:  144
Length: 4.050578, Radius: 0.024751
No. of Wall: 1
Total No. of Nodes: 11
Node: 	X	Y	Z(MEDIUM) 	L	Index	VTG	CND		Power
0: 	9.61	2.526	3.630(MED)	0.000	0	0.997	4.276e-05	3.105e-21
1: 	9.32	2.495	3.288(MED)	0.450	1	0.997	4.276e-05	1.385e-20
2: 	8.98	2.436	3.001(MED)	0.900	2	0.997	4.276e-05	8.794e-22
3: 	8.56	2.453	2.831(MED)	1.350	3	0.997	4.276e-05	1.034e-21
4: 	8.16	2.508	2.648(MED)	1.800	4	0.997	4.276e-05	1.311e-21
5: 	7.87	2.633	2.322(MED)	2.250	5	0.997	8.624e-05	2.689e-20
6: 	7.76	2.679	2.132(MED)	2.474	6	0.997	8.483e-05	1.513e-22
7: 	7.66	2.726	1.939(MED)	2.700	7	0.997	4.276e-05	1.163e-21
8: 	7.37 	2.761	1.591(MED)	3.150	8	0.997	4.276e-05	3.349e-22
9: 	7.16	2.950	1.244(MED)	3.601	9	0.997	4.276e-05	4.051e-22
10:	6.93	3.111	0.887(MED)	4.051	10	0.997

Address: 1
Parent:  1
Length: 2.262174, Radius: 0.025742
No. of Wall: 1
Total No. of Nodes: 11
Node: 	X	Y	Z(MEDIUM) 	L	Index	VTG	CND		Power
0: 	-8.68	8.51	-8.318(MED)	0.00	0	0.00	8.283e-05	0.000e+00
...
```



© Copyright 2019





















