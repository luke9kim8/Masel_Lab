from Bio import Phylo
from TreePruner import TreePruner
from Codeml_runner import Codeml_runner

class TreeCombinator:
    def __init__(self, domain_tree_path, pruned_species_tree_path, format):
        """
        Constructor for Treecombinator class
        :param domain_tree_path:
        :param pruned_species_tree_path:
        :param format: format of the tree
        """

        self.species_tree = Phylo.read(pruned_species_tree_path, format)
        self.domain_tree = Phylo.read(domain_tree_path, format)

    def get_species_list(self):
        """

        :return: list of species in the domain_tree
        """
        list = []
        for clade in self.domain_tree.get_terminals():
            name = self.get_species_name(clade)
            list.append(name)

        return list

    def get_domain_tree_with_speciation_bl(self):
        """

        :return: domain tree with speciation branch lengths obtained from the pruned species tree
        """

        self.add_terminal_branch_lengths()
        self.add_nonterminal_branch_lengths()
        Phylo.write(self.domain_tree, 'domain_tree_wSpeciation.nwk','newick')

    def get_speciation_dictionary(self):
        """

        :return: get a dictionary of [species name: speciation date]
        """
        species_dictionary = next(self.species_tree.find_clades()).depths()
        return species_dictionary

    def get_species_name(self, domain):
        """
        removes unique identifiers of domains and returns their species name
        :param domain:
        :return: the species name of the domain tree
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

    def get_domain_names(self, domain_tree_path):
        """

        :param domain_tree_path:
        :return: list of all domain names from the domain tree
        """
        tree = Phylo.read(domain_tree_path, 'newick')
        list_of_domain_names = []
        for terminals in tree.get_terminals():
            list_of_domain_names.append(terminals.name)
        return list_of_domain_names

    def add_terminal_branch_lengths(self):
        """
        add terminal branch lengths to the domain tree
        :return:
        """
        domain_terminals = self.domain_tree.get_terminals()

        for nonterminal in domain_terminals:
            species_name = self.get_species_name(nonterminal)
            for clade in self.get_speciation_dictionary():
                if clade.name == species_name:
                    nonterminal.branch_length = clade.branch_length

    def organize_species_list(self, clade_list):
        """
        organize species list so that there are no duplicates
        :param clade_list:
        :return:
        """
        species_list = []
        for clade in clade_list:
            species_name = self.get_species_name(clade)
            if species_name in species_list:
                pass
            else:
                species_list.append(species_name)
        species_list.sort()
        return species_list

    def duplication(self, nonterminal):
        """
        check if the node is duplication
        :param nonterminal:
        :return:
        """
        child1 = nonterminal.clades[0].get_terminals()
        child2 = nonterminal.clades[1].get_terminals()

        # domain_children = organize_species_list(nonterminal.get_terminals())

        child1_children = self.organize_species_list(child1)
        child2_children = self.organize_species_list(child2)

        for child1 in child1_children:
            for child2 in child2_children:
                if child1 == child2:
                    return True
                else:
                    return False

    def get_duplication_clades(self, tree):
        """

        :param tree:
        :return: list of clades(or nodeds) that are duplications
        """
        duplication_clade_list = []
        nonterminals = tree.get_nonterminals()

        for nonterminal in nonterminals:
            if self.duplication(nonterminal):
                duplication_clade_list.append(nonterminal)

        return duplication_clade_list

    def share_ancestor(self, children_list1, children_list2):
        """

        :param children_list1:
        :param children_list2:
        :return: true if children_list1 has the same children as children_list2
        """
        for child1 in children_list1:
            for child2 in children_list2:
                if child1 == child2:
                    return True
                else:
                    return False

    def add_nonterminal_branch_lengths(self):
        """

        :return: add nonterminal branch length to the domain tree
        """
        for nonterminal in self.domain_tree.get_nonterminals():
            child1 = nonterminal.clades[0].get_terminals()
            child2 = nonterminal.clades[1].get_terminals()

            # domain_children = organize_species_list(nonterminal.get_terminals())

            child1_children = self.organize_species_list(child1)
            child2_children = self.organize_species_list(child2)

            if self.duplication(nonterminal):
                pass
            else:
                for ancestor in self.species_tree.get_nonterminals():
                    ancestor_children1 = self.organize_species_list(ancestor.clades[0].get_terminals())
                    ancestor_children2 = self.organize_species_list(ancestor.clades[1].get_terminals())
                    if self.share_ancestor(ancestor_children1, child1_children) \
                            and self.share_ancestor(ancestor_children2, child2_children) \
                            or self.share_ancestor(ancestor_children1, child2_children) \
                            and self.share_ancestor(ancestor_children2, child1_children):
                        nonterminal.branch_length = ancestor.branch_length

        return self.domain_tree

    def add_default_branch_lengths(self):
        """
        for the clades(or nodes) without a branch length, set it to 30.0
        :return:
        """
        for clade in self.domain_tree.get_nonterminals():
            if clade.branch_length == None:
                clade.branch_length = 30.0

    def save_tree(self, path):
        """
        save the combined tree
        :param path:
        :return:
        """
        Phylo.write(self.domain_tree, path, 'newick')

    def eliminate_blank_space(self, text):
        """

        :param text:
        :return: returns a text without any spaces
        """
        result = ''
        for char in text:
            if char == ' ':
                result = result + ''
            else:
                result = result + char
        return result

    def get_tree_from_codeml_file(self, codeml_file_path, tree_path):
        """

        :param filename:
        :return newick_tree: read codeml output file and find the tree that matches with the domain tree
        """

        name_list = self.get_species_list()
        with open(codeml_file_path, 'r') as result_file:
            result_file_lines = result_file.readlines()
            for line in result_file_lines:
                line = self.eliminate_blank_space(line)
                result = all(elem in line for elem in name_list)
                if result:
                    with open(tree_path, 'w') as result_newick_tree_file:
                        result_newick_tree_file.write(self.eliminate_blank_space(line))

    def get_branch_length_ratio(self, clade):
        """

        :param clade:
        :return: ratio of branch lengths from the codeml output file
        """
        # get clade branch length
        clade_bl = clade.branch_length
        # get avg of 2 children clades' branch length
        children_bl = (clade.clades[0].branch_length + clade.clades[1].branch_length) / 2.0

        ratio = (clade_bl) / (children_bl + clade_bl)

        return ratio

    def show_domain_tree(self):
        Phylo.draw(self.domain_tree, branch_labels=lambda c: c.branch_length)

    def show_species_tree(self):
        Phylo.draw(self.species_tree, branch_labels=lambda c: c.branch_length)

    def resize_branch_length(self, paml_tree, domain_tree):
        """
        adds the branch length to duplication nodes in the domain tree
        :param paml_tree:
        :param domain_tree:
        :return:
        """
        paml_dup_list = self.get_duplication_clades(paml_tree)
        domain_dup_list = self.get_duplication_clades(domain_tree)

        print(paml_dup_list)
        print(domain_dup_list)

        if len(paml_dup_list) != len(domain_dup_list):
            raise ValueError('Number of duplication in PAML tree and Domain tree are not equal')
        else:
            num_of_dup = len(paml_dup_list)

        print(num_of_dup)
        for index in range(num_of_dup):
            if paml_dup_list[index].branch_length is None:
                index += 1
            else:
                bl_ratio = self.get_branch_length_ratio(paml_dup_list[index])
                print(bl_ratio)
                '''
                change duplication node branch length in the domain tree for clade in domain tree
                '''
                print('hey')
                for clade in domain_tree.get_nonterminals():
                    print(clade.clades[0].branch_length, clade.clades[1].branch_length)

                    if self.duplication(clade):
                        print(clade.clades[0].branch_length)
                        speciation_time = (clade.clades[0].branch_length + clade.clades[1].branch_length) / 2.0
                        clade.branch_length = speciation_time * bl_ratio
                        clade.clades[0].branch_length = speciation_time * (1 - bl_ratio)
                        clade.clades[1].branch_length = clade.clades[0].branch_length
                        print(clade.branch_length)
            index += 1

        return domain_tree


def main(domain_tree_path,  species_tree_path, seqfile_path, outfile_path, format):
    """
    takes a domain tree without dated branch lengths and add speciation and duplication dates
    based from codeml and speciation tree (units in MYA)
    :param domain_tree_path:
    :param species_tree_path:
    :param seqfile_path:
    :param format:
    :return:
    """
    pruned_species_tree_path = 'Input/pruned_domain_tree.nwk'
    treefile_path = domain_tree_path

    #prune_trees
    pruner = TreePruner(species_tree_path, 1)
    pruner.prune_trees(domain_tree_path,pruned_species_tree_path)

    #make domain_tree_with_speciation_branch_lengths
    tree = TreeCombinator(domain_tree_path, pruned_species_tree_path, format)
    tree.get_domain_tree_with_speciation_bl()
    domain_tree_wSpeciation = Phylo.read("domain_tree_wSpeciation.nwk", 'newick')

    #run codeml
    cml = Codeml_runner()
    cml.run_codeml(seqfile_path, treefile_path)
    tree.get_tree_from_codeml_file('Codeml/Codeml_files/codeml.file', 'Codeml/Codeml_trees/codeml.tree')

    #resize branch length based on codeml
    codeml_tree = Phylo.read('Codeml/Codeml_trees/codeml.tree','newick')
    result = tree.resize_branch_length(codeml_tree,domain_tree_wSpeciation)
    result_name = outfile_path
    Phylo.write(result, result_name, 'newick')


if __name__ == '__main__':
    domain_tree_path = 'Input/PF15897.reconciled'
    species_tree_path = 'Input/SpeciesTree.nwk'
    seqfile_path = 'Input/RAxML_bipartitions.PF15897.fasta.align'
    outfile_path = 'output.tree'
    format = 'newick'

    main(domain_tree_path, species_tree_path, seqfile_path, outfile_path, format)
# tree.add_branch_length()
# Phylo.draw(tree.domain_tree, branch_labels=lambda c: c.branch_length)