# tmQM_meth_methanol
Houses code for high-throughput screening of homogeneous catalysts for activating methane


Everybody is welcome to work on the tasks but we are segmenting in order to parallelize the execution.

# Yuri TODO List

- [x] Script that parses the denticity of the complex to identify the halogen singleton in the molecular graph: Filter out to 18 k candidates suitable activating methane
- [ ] Script that creates the open coordination at the halogen coordinate, placing the intermediates (oxo, hydroxyl, methanol, resting) generating a set of .xyz files available for further processing and electronic structure calculations. We need to apply certain geometric constraints so that the pre-optimization (with xTB) and the subsequent single point calculation (with PBE0/TZVP) is successful.
- [ ] We need to build an intelligent framework, using Fireworks, for parallelizing the relaxations + single point calculations against a Kubernetes cluster.


# Hoon TODO List
- [ ] Now that we have the 18k candidates, we need to vectorize their representation and using these fingerprints, figure out what the most appropriate basis to represent the complexes, which we will use to downsample to the set of candidates that will undergo single point calculations.
- [ ] PCA to downsample



# Nick TODO List
- [ ] Define the model type (GP, GNN?) architecture and target type (delta Oxo/ delta HAT)



