(In progress, this should be removed when project working)

This project is designed to make use of the results from TreeGrafter (https://github.com/pantherdb/TreeGrafter). TreeGrafter when provided with protein information determines the PANTHER phylogenetic gene tree and the location in that tree at which a given protein best fits. 

The results output by TreeGrafter however don't take into account species. This program will update the location of the parent suggested by TreeGrafter to the parent node which represents the closest ancestral species. Additionally a list of all proteins / gene products in the given PANTHER tree which share an ancestral speciation node with the grafted protein will be output. 
