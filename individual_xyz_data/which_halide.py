import autode as ade 
from autode.mol_graphs import split_mol_across_bond

ade.Config.n_cores = 3
#def get_halides(neighbors):

    #halides = ['F', 'Cl', 'Br', 'I']
    
    #neighbors = neighbors.replace("'", "")
    
    #neighbor_lst = neighbors.split(', ')
    
    #halides = [atom for atom in neighbor_lst if atom == 'F' or atom == 'Cl' or atom == 'Br' or atom == 'I']
    
    #return halides
    
def get_halides(neighbors):
    halides = []
    for Atom in neighbors:
        if Atom[0].atomic_symbol == 'F' or Atom[0].atomic_symbol == 'Cl' or Atom[0].atomic_symbol == 'Br' or Atom[0].atomic_symbol == 'I':
            halides.append(Atom)
    return halides



# Get list of selected Halides 
with open('selected_halides.txt') as f:
        lines = f.read().splitlines()



with open('halides.txt') as f:
    neighbor_information = f.read().splitlines()


for complex_name in lines:
    print(complex_name)
    with open(complex_name) as f:
        information = f.readlines()[1]
        
    # Extract multiplicity and formal charge information 
    q_index = information.index('q =') + 4
    S_index = information.index('S =') + 4
    q, S = int(information[q_index:q_index + 2]), int(information[S_index])
    


    for line in neighbor_information:
        if complex_name in line:
            neighbors = line[line.index('[') + 1:-1]
            break
    

    #halides = get_halides(neighbors)
    
    

    # Generate autodE molecule with xyz coordinates, multiplicity, charge 
    mol = ade.Reactant(complex_name, mul = 2 * S + 1, charge = q)
    metal_center = [(index, atom) for index, atom in enumerate(mol.atoms) if atom.is_metal]
    metal_center_index = metal_center[0][0]
    neighbors = list(mol.graph.neighbors(metal_center_index))
    
    
    
    atom_Lst = []
    
    for neighbor in neighbors:
        
        atom_Lst.append([mol.atoms[neighbor], neighbor])
    
    halide_Lst = get_halides(atom_Lst)
    
    if len(halide_Lst) > 1:
        try:
            print(halide_Lst)
            for halide in halide_Lst:

                complex_nodes, halide_nodes = split_mol_across_bond(mol.graph, bond=(metal_center_index, halide[1]))
                oxo_interm = ade.Product(name='oxo_interm', mult = 2 * (S + 1/2) + 1, charge = q, atoms=[mol.atoms[i] for i in complex_nodes])
                halide_interm = ade.Product(name='halide_ion', mult = 2 , charge = 0, atoms=[mol.atoms[i] for i in halide_nodes])
            
                #rxn = ade.Reaction(mol, oxo_interm, halide_interm)
                for species in (mol, oxo_interm, halide_interm):
                    species.single_point(method = ade.methods.XTB())
                #rxn.optimise_reacs_prods()
                #rxn.locate_transition_state()

                #energy = rxn.delta('E‡').to('kcal mol-1')
                
                print(627.5*(oxo_interm.energy + halide_interm.energy - mol.energy)) 
                
            
                
                #print(complex_nodes, halide_nodes)
        except:
            print("Error: Skipping to next structure...")
            continue


    

        
    

