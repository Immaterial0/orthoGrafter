(In progress, this should be removed when project working)

This project is designed to make use of the results from TreeGrafter (https://github.com/pantherdb/TreeGrafter). TreeGrafter when provided with protein information determines the PANTHER phylogenetic gene tree and the location in that tree at which a given protein best fits. 

The results output by TreeGrafter however don't take into account species. This program will update the location of the parent suggested by TreeGrafter to the parent node which represents the closest ancestral species. Additionally a list of all proteins / gene products in the given PANTHER tree which share an ancestral speciation node with the grafted protein will be output. 

Updating the graft point to the closest related species requires inputting both a tree, initial graft point, and the species from which the protein comes (species can be input as either the short 5 character species names or as the ncbi taxonomy id number). The closest related ancestor specis is found using a precomputed taxon data table for child and parent. 

If no save location is provided, the file will be saved to updatedTreeGrafterPoint.out

Output results should be stable relative to the order of input results from TreeGrafter (note: TreeGrafter results are not necessarily stable themselves however). 
