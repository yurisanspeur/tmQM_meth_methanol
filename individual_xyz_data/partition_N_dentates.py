import autode as ade
import itertools
from networkx.algorithms.simple_paths import all_simple_paths
import sys
import glob


structs = glob.glob("*.xyz")

for struct in structs:

    try:
        # Load the molecule

        mol = ade.Molecule(struct, mul=1, charge=1)

        if mol.n_atoms > 100: # Too large
            continue

        # Identify all the immediate metal node degrees (MND)

        #mnds = [y for (x,y) in mol.graph.edges if x == 0] # Here we hardcoded the metal center but can get from is_metal

        #print(mnds)
        # Get the index of the TM metal center
        metal_center = [(index, atom) for index, atom in enumerate(mol.atoms) if atom.is_metal]
        assert len(metal_center) == 1
        neighbors = list(mol.graph.neighbors(metal_center[0][0])) #FIXME: Sometimes this does not pick up edges based on the default 1.3 x avg bond_length


        # We want to check which pairs of nodes are connected so as to determine whether the ligand is monodentate or bidentate,
        # which will help us decide how to create the open coordination
        pairs = itertools.combinations(neighbors, 2)


        connected_pairs = []
        connected_neighbor_paths = []


        for pair in pairs:
            simple_paths = [p for p in all_simple_paths(mol.graph, pair[0], pair[1]) if metal_center[0][0] not in p] # don't navigate through metal center
            # Do the union of simple_paths
            connected_neighbor_path = set().union(*simple_paths)
            connected_pair = connected_neighbor_path.intersection(set(neighbors))
            if connected_pair:
                connected_pairs.append(connected_pair)
                connected_neighbor_paths.append(connected_neighbor_path)


        # Now analyze the results to partition into monodentate, bidentate, tridentate etc...
        partitions = {}
        for i, connected_pair in enumerate(connected_pairs):
            intersection = set(connected_pair).intersection(set().union(*list(partitions.values())))
            if f"ligand-{i}" not in partitions and len(intersection) == 0: # Create a new ligand class
                partitions[f"ligand-{i}"] = connected_pair
            else:
                for k, v in partitions.items(): # Assign the intersection to the key and update the value of that key with the union of the connected_pair
                    if intersection.intersection(v):
                        # Update v with union
                        partitions[k] = v.union(connected_pair)
                


        partitions['monodentates'] = set(neighbors) - set().union(*list(partitions.values()))    
        #print(sys.argv[1], partitions, neighbors)
        labels = [mol.atoms[idx].label for idx in partitions['monodentates']]
        if partitions['monodentates'] and set(labels).intersection(['F', 'Cl', 'Br', 'I']):
            print(f"monodentates for {struct} are: {labels}")
        else:
            print("Discarding complex...")
    except:
        print("Error: Skipping to next structure...")
        continue








