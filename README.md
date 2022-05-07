# tmQM_meth_methanol
Houses code for high-throughput screening of homogeneous catalysts for activating methane


Everybody is welcome to work on the tasks but we are segmenting in order to parallelize the execution.

# Yuri TODO List

- [x] 1. Script that parses the denticity of the complex to identify the halogen singleton in the molecular graph: Filter out to 18 k candidates suitable activating methane
- [x] 2. Script that creates the open coordination at the halogen coordinate, placing the intermediates (oxo, hydroxyl, methanol, resting) generating a set of .xyz files available for further processing and electronic structure calculations. We need to apply certain geometric constraints so that the pre-optimization (with xTB) and the subsequent single point calculation (with PBE0/TZVP) is successful.
- [x] 3. We need to build an intelligent framework, using Fireworks, for parallelizing the relaxations + single point calculations against a Kubernetes cluster.


# Hoon TODO List
- [x] 4. Now that we have the 18k candidates, we need to vectorize their representation and using these fingerprints, figure out what the most appropriate basis to represent the complexes, which we will use to downsample to the set of candidates that will undergo single point calculations.
- [x] 5. PCA to downsample 



# Nick TODO List
- [ ] 6. Define the model type (GP, GNN, multi-fidelity?) architecture and target type (delta Oxo/ delta HAT) (We can show the workflow on the original tmQM dataset with their calculated properties if finding the right charge and multiplicity is proving to be a bottleneck)
- [x] 7. For the 18k candidates, we need a script that figures out the correct multiplicity and charge to apply to the ORCA input files (Seems like there is not much difference)
- [ ] 8. Script to figure out which of the halogen atoms is most likely to be popped out in the case that there are more than one (Use low level theory like xTB to corroborate). This should be used in 2. (WIP: We have a script that uses Pymatgen to figure out the oxidation state of the metal center and together with the coordination (MND) and the spin_configuration ('LS', 'MS', 'HS'), is able to return the number of unpaired electrons in the metal center, which gets used to calculate spin. Molsimplify NN can help validate (Fe(II), Fe(III))



