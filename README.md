# virus_polymorphism_homophily

Example run for a simulation (single deme, used in most part of the paper)
```
gcc -lgsl simulation.c -o simulation.out
./simulation.out 100000 0.03 0.03 0.3 0.6 100000 1
```

Example run for a simulation with two demes
```
gcc -lgsl simulation_two_demes.c -o simulation_two_demes.out
./simulation_two_demes.out 50000 0.02 0.02 0.2 0.8 0.5 0.5 0.1 100000 1
```
