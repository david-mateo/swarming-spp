## Running the examples
Some example programs using the `spp` library are provided under `examples/`. Each folder will contain
* One or more C++ source code files to run different calcualtions (for example, one file for metric interaction and one for topological interaction.)
* A Makefile to compile the source code, with may or may not require additional parameters to compile such as the noise level or the outdegree of the interaction.
* A Bash script to run the program for a range of possible parameters.

### Order parameter in Vicsek's model
This example computes the order parameter of a collection of 5000 self-propelling particles following the Vicsek model. The order parameter is defined as the norm of the average velocity divided by the average speed. The parameters of this calculation are chosen to reproduce as closely as possible the results presented in [Vicsek et al. PRL 75, 1226 (1995)](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.75.1226).


####Run one case
Navigate to `examples/order/` and type

```
  make vicsek_metric eta=0.1
```

to create an executable called `vicsek_metric` that computes the order parameter of a swarm using metric interaction (fixed interaction radius of 1.0) with noise level 0.1 (10% noise). To create the equivalent with topological interactioni (fixed outdegree of 7), type

```
  make vicsek_topo eta=0.1
  ```
  
Running `vicsek_metric` or `vicsek_topo` will run for `NITER`=10000 iterations and print on the screen the values of the order parameter every `OUTPUT`=10 iterations, after running for `TRANSIENT`=1000 iterations.

####Run all cases
To compute the order parameter for a range of noise values, run the script `run_metric_serial.sh`. This will compile and execute the program for noise levels 0.05, 0.10, 0.15 ... 1.0 . The results will be stored in `logs/metric_n{n}.res`, where `{n}`is the noise level.
To run the different noise levels in parallel via the `qsub` command, use `run_metric_pbs.sh` instead.

### Correlations and susceptibility
This example computes the correlations in velocity fluctuations in a collective of 2048 self-propelling particles following the Vicsek model. Because this calculation requires a large number of iterations to obtain statistically significant results, the computational cost of this is considerably higher than in other examples. This example showcases how to use the `Grid` class in conjunction with `Community` to significantly reduce the computation cost by storing information about which agents are 'in the neighborhood' (see `Grid` documentation for more info).
The correlation is computed following the framework presented in [Attanasi et al PLoS Comput Biol 10, e1003697 (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003697) using the `Community::correlation_histo`
function.

The output of this program requires non-trivial post-processing. The program prints a collection of 'correlation histograms' in the form of

```
distance    Total_Correlation    N
```

so that the correlation `C(r)` is `C(r=distance)=Total_Corelation/N`. Each histogram is preceded by a line containing `#Iteration %i` and followed by two blank lines.

####Run one case
Navigate to `examples/susceptibility/` and type

```
  make vicsek_metric_r{R}
```

where `{R}` is the desired value for the interaction radius. This will create an executable called `vicsek_metric_r{R}` that computes the correlation of a swarm using metric interaction with radius `{R}` and fixed noise level 0.05. To create the equivalent with topological interaction of outdegree `{K}`, type

```
  make vicsek_topo_k{K}
```

####Run all cases
The scripts `run_metric.sh` and `run_topo.sh` will submit jobs using `qsub` to compute the correlation for a range of radius and outdegrees respectively.
The results will be stored in `logs/metric_den1/correlation_r{R}.res` and `logs/topo_den1/correlation_k{K}.res`.

### Predator attack
This example simulates the attack of a predator in a swarm of 2048 self-propelling particles (preys). This calculation uses the class `HostileEnvironment`, an extension on `Community` that incorporates the predator-prey dynamic.

The behaviors of the predator and the preys are extensions of the Vicsek heading consensus. The predator behavior, implemented in the `Vicsek_predator` class, makes the predator move as fast as possible towards the nearest prey. If the predator catches a prey (meaning that the prey closer that the distance the predator can move in one iteration) the `HostileEnvironment` takes care of removing the prey from the swarm. The prey behavior, implemented in the `Vicsek_prey` class, makes the prey behave exactly as `Vicsek_consensus` would unless a predator is found in the prey's detection area. If a predator is detected, the prey aligns its velocity to flee from the predator.

The program outputs the number of iterations ellapsed between two consecutive "catches" of the predator (avoidance time).

####Run one case
Navigate to `examples/predator/` and type

```
  make predator_metric_r{R}
```

where `{R}` is the desired value for the interaction radius. This will create an executable called `predator_metric_r{R}` that simulates a predator attack on a swarm and outputs the avoidance times.
All the examples presented here allow the random seed to be passed as an argument on run, for example `./predator_metric_r1.4 53452345236`. In the case of the predator attack, it is mandatory to provide such argument. This is because the calculation of the mean avoidance time requires of a large sample of runs and it is imperative to have a good sampling of the initial configuration space for the whole swarm.

