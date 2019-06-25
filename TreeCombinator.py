from Bio import Phylo

domain_tree = Phylo.read('PF00457.nwk', 'newick')              #input domain newick file
domain_terminals = domain_tree.get_terminals()                 #get list of terminal nodes(clades) in domain tree
domain_nonterminals = domain_tree.get_nonterminals()           #get list of domain nonterminal nodes

species_tree = Phylo.read('Pruned_SP_tree.nwk', 'newick')      #input species newick file
species_nonterminals = species_tree.get_nonterminals()         #get list of nonterminal nodes

species_dictionary = next(species_tree.find_clades()).depths() #get a dictionary of [species name: speciation date]

#trim the unique identifier of the domain
#EX) Aspergillus_oryzae|11096171-PF00457 -----> Aspergillus_oryzae
def trim_unique_identifier(domain):
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



#add species tree terminal branch lengths to domain tree terminals
for nonterminal in domain_terminals:
    species_name = trim_unique_identifier(nonterminal)
    for clade in species_dictionary:
        if clade.name == species_name:
            nonterminal.branch_length = clade.branch_length

#get list of species of the list of domain nodes
def organize_species_list(clade_list):
    species_list = []
    for clade in clade_list:
        species_name = trim_unique_identifier(clade)
        if species_name in species_list:
            pass
        else:
            species_list.append(species_name)
    species_list.sort()
    return species_list

#check if the node is domain duplication
def duplication(children1, children2):
    for child1 in children1:
        for child2 in children2:
            if child1== child2:
                return True
            else:
                return False

#check if the nodes shares the same ancestor (same as duplication)
def share_ancestor(children_list1, children_list2):
    for child1 in children_list1:
        for child2 in children_list2:
            if child1 == child2:
                return True
            else:
                return False

#add species tree nonterminal branch lengths to domain tree nonterminals
for nonterminal in domain_nonterminals:
    child1 = nonterminal.clades[0].get_terminals()
    child2 = nonterminal.clades[1].get_terminals()

    domain_children = organize_species_list(nonterminal.get_terminals())

    child1_children = organize_species_list(child1)
    child2_children = organize_species_list(child2)

    if duplication(child1_children, child2_children):
        pass
    else:
        for ancestor in species_nonterminals:
            ancestor_children1 = organize_species_list(ancestor.clades[0].get_terminals())
            ancestor_children2 = organize_species_list(ancestor.clades[1].get_terminals())
            if share_ancestor(ancestor_children1, child1_children) \
                    and share_ancestor(ancestor_children2, child2_children)\
                    or share_ancestor(ancestor_children1,child2_children) \
                    and share_ancestor(ancestor_children2, child1_children):
                nonterminal.branch_length = ancestor.branch_length


#for nodes that does not have a branch length, set it to 30.0 (just for visualization)
def add_default_branch_lengths():
    for clade in domain_tree.get_nonterminals():
        if clade.branch_length == None:
            clade.branch_length = None

add_default_branch_lengths()

Phylo.draw(domain_tree, branch_labels=lambda c: c.branch_length)
Phylo.write(domain_tree,"combined_tree.nwk",'newick')
