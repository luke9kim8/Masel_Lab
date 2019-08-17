from ete3 import Tree
from Bio import Phylo
# t = Tree('input/SpeciesTree.nwk', format = 1)
# pruned_tree = t.prune(['Aspergillus_oryzae', 'Aspergillus_clavatus', 'Fusarium_fujikuroi'])
# t.write(format=1, outfile="new_tree.nw")

class TreePruner:
    def __init__(self, tree_name, format):
        self.tree = Tree(tree_name, format)

    def get_species_name(self, domain):
        """
        removes unique identifiers of domains and returns their species name
        :param domain:
        :return:
        """
        name = domain.name
        if name == None:
            species_name = ''
            return species_name
        elif '|' in name:
            bar_index1 = name.find(r'|')
            species_name = name[:bar_index1]
            return species_name
        else:
            species_name = domain.name
            return species_name

    def get_species_list(self, domain_tree_path):
        """

        :return: list of species in the domain_tree
        """
        domain_tree = Phylo.read(domain_tree_path,'newick')
        list = []
        for clade in domain_tree.get_terminals():
            name = self.get_species_name(clade)
            list.append(name)
            print(list)

        return list

    def prune_trees(self, domain_tree_path, output_path):
        """
        returns a newick file with a tree containing only the species in the list
        :param list_of_species:
        :param output_path:
        :return:
        """
        list_of_species = self.get_species_list(domain_tree_path)
        print(list_of_species)
        self.tree.prune(list_of_species)
        self.tree.write(format=1, outfile=output_path)

