import os, sys, argparse
from time import gmtime
from time import strftime 
from typing import List

"""
taxonomy_tree.py

Step 1.)

Takes the nodes.dmp, names.dmp, and the accession2taxid as the input and creates a taxonomy tree
that we can later use to determine the relatives of a certain taxid.

References
----------
This code is inspired from Jennifer Lu's code linked below
https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py

How to run
----------

$ python src/taxonomy_tree.py

The above line is for if you want to call and test this script directly.

Otherwise, this script is designed to be part of the larger program in main.py,
so it will be automatically used in calls to main.py and other files in the program.

Authors
-------

Computational Genomics Team 47:
Dhruv Dubey
Mitra Harpale
Christopher Li
Jaeyoon Wang

"""

# TaxaTree data structure representing nodes in taxonomy
class TaxaTree:
    def __init__(self, tax_id : int, rank, parent = None, children = [], name = "", isRoot = False):
        # Node's own taxid, name, and rank
        self.tax_id : int = tax_id
        self.name = name
        self.rank = rank

        # Setting Parent
        self.parent = parent
        self.parentTaxId : int = -1

        # Children
        self.children = children
        self.childrenTaxId : List[int] = []


    def add_child(self, child_node):
        # Two-way connection between nodes!
        self.children.append(child_node)

    def isRoot(self):
        self.isRoot = True

    def get_parent(self):
        return self.parent

    def get_children(self):
        return self.children

    def get_tax_id(self):
        return self.tax_id

# for utility, for convenience
ranks_to_charRank = {
    "superkingdom": "SK",
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
    "subspecies": "SS",
    "subfamily": "SF",
    "no rank": ""
}

# initializing the main dictionaries which let us navigate
# and traverse the taxonomic tree

# tax id to node (where node is a TaxaTree object)
# https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py
taxonomy_id_to_node = {}

# a tax id mapped to its parent's tax id
# https://github.com/DerrickWood/kraken/blob/master/src/krakenutil.cpp
taxonomy_id_to_parent_id = {}

# new dictionary for the unlinked parents, the node's whose parents we haven't seen yet
# in chronological order as we go down the tree
parents_unseen = {}

# redundancy in the dictionaries let us navigate both up and down the tree
# following forward pointers and backpointers in the n-ary tree like a deque
# i.e. bidirectional links

def build_parent_map():
  """
  This method closely implements
  https://github.com/DerrickWood/kraken/blob/master/src/krakenutil.cpp
  
  Which basically extracts the first two columns from nodes.dmp and puts them into
  the taxonomy_id_to_parent_id map.
  
  Keyword arguments:
  argument -- description
  Return: return_description
  """
  
  nodes_dmp_file_handle = open('./taxonomy/nodes.dmp', 'r')
  # finish soon
  
  count_num_nodes = 0
  root_node = None # to be updated later with the first line containing the root node with taxonomy id 1
  # Iterating over the nodes.dmp (1st column in t)xId, 2nd column is the 1st col's parent's taxId
  for line in nodes_dmp_file_handle:
    # Increment total number of nodes we've seen so far
    count_num_nodes += 1

    if count_num_nodes % 1000 == 0:
       sys.stdout.write(f"\r\t{count_num_nodes} lines processed so far")

    # Spliting the line we are on
    tokens = line.strip().split('\t|\t')
    # print(tokens)

    # Storing info in variables
    current_taxonomy_id = tokens[0]
    parent_taxonomy_id = tokens[1]
    rank = tokens[2] # rank is one of "kingdom", "phylum", "class", "order", etc.
    
    abbreviated_rank = ""
    if rank in ranks_to_charRank.keys():
      abbreviated_rank = ranks_to_charRank[rank] # "K", "P", "C", "O", etc.

    # TODO: Is the current_node being created correctly? What redundancy can we remove here?
    # Create a node with the taxid and its rank
    current_node = TaxaTree(current_taxonomy_id, abbreviated_rank)
    current_node.parentTaxId = parent_taxonomy_id
    
    # nodes themselves have backpointers to their taxnonomy ids and their parent nodes
    # so we have redundancy in links like forward pointers and backpointers
    # in the n-ary tree structure

    # This is Jen's map - maps taxonomy ids to nodes
    taxonomy_id_to_node[current_taxonomy_id] = current_node

    # This is Derrick's map - maps taxonomy ids to parent taxonomy ids
    taxonomy_id_to_parent_id[current_taxonomy_id] = parent_taxonomy_id

    # three cases below:

    if current_taxonomy_id == "1":
      # root node
      current_node.isRoot()
      # save this
      root_node = current_node

    elif parent_taxonomy_id in taxonomy_id_to_node.keys():
      # If we've already seen and created the parent node before in the nodes.dmp file,
      # i.e. if we found the parent node at an earlier line
      # (a line above the current one in the nodes.dmp file)
      # then this is great because we can link up the nodes bidirectionally (like a deque)
      current_node.parent = taxonomy_id_to_node[parent_taxonomy_id] # query the dict for the parent node
      taxonomy_id_to_node[parent_taxonomy_id].add_child(current_node) # set current node as the child of the parent
      
    else:
      # If we haven't seen the parent node's taxnomy id in the nodes.dmp file,
      # i.e. if the parent comes somewhere later / below in the file at a later line
      # then we have to keep track of it as a separate disconnected component in the graph / forest
      # for now, and we will link up the disconnected islands at some point into one final
      # connected graph at a later step
      # i.e. we can only link up one direction for now
      # save this for later fixing
      parents_unseen[current_taxonomy_id] = current_node
      # add update the node's pointers as normal, but no opposite direction update
      # since the parent does not exist as of yet
      current_node.parentTaxId = parent_taxonomy_id
      # line below unneeded for now
      # current_node.tax_id = current_taxonomy_id

  nodes_dmp_file_handle.close()
  print("\nAll lines processed in the nodes.dmp file")

  # Create any parents that were unseen or appeared out of order in the nodes.dmp file
  for tax_id in parents_unseen:
    # get the current node whose parents were unseen (orphan)
    orphan_node = parents_unseen[tax_id]
    # check again whether the parent is now in the set (i.e. if they appeared
    # out of order in the nodes.dmp file)
    # if the parent node is actually in the set (i.e. was created and found later in the file),
    if orphan_node.parentTaxId in taxonomy_id_to_node:
      # then great, we found the parent, and the orphan is
      # no longer an orphan in the tree, in terms of nodes!
      # we can go ahead and update bidirectional links of the child
      # and parent nodes as in the normal case:
      orphan_node.parent = taxonomy_id_to_node[orphan_node.parentTaxId] # make the child no longer an orphan
      taxonomy_id_to_node[orphan_node.parentTaxId].add_child(orphan_node) # make the parent no longer childless      
    else:
      # If the parent node wasn't found anywhere in the file, even later, then this is
      # truly an orphan node and we canont resolve this issue, so we report it
      print(f"""A parent node of orphaned tax_id node {tax_id} was not found
             in the nodes.dmp file: the unfound parent is {orphan_node.parentTaxId}""")
  
  # Moving on to the next step, which is pruning unused branches of the taxonomy tree
  # i.e. the taxonomy tree is not deep (~7 ranks: kingdom, phylum, ..., genus, species)
  # but it is very, very wide, it is the taxonomic tree of the entire tree of life on Earth
  # The tree is so big that it would be advantageous to include only those root to leaf paths
  # that are actually used and hit by the sequences in the reference genomes (~20 genomes), so
  # this step may involve pruning the taxonomic tree to make it small enough for ideally
  # a smaller memory usage

  # TODO prune taxonomic tree and output pruned tree information to a file to be read later


  

  # and return the gathered tree information about the
  # taxonomic tree so that it doesn't have to be loaded in memory again
  # TODO make sure this is returning everything that needs to be returned
  return taxonomy_id_to_node, taxonomy_id_to_parent_id, root_node


def get_fasta_ncbi_accession_ids(database_directory : str) -> List[str]:
  """
  Goes through the genomes-of-common-contaminants directory and gets all of the accession IDs
  of the FASTA genome files.

  @param: database_directory (a string of the folder where the reference genomes are
    fed in through the command line argument --db)
  
  @return: a list of NCBI accession ID strings
  """
  ncbi_accession_ids = []
  
  # for each FASTA file
  for f in os.listdir(database_directory):
    
    # get the NCBI accession id for this FASTA file
    with open(os.path.join(database_directory, f), 'r') as fasta_file_handle:
      # getting the first token of the first line, which is the NCBI accession ID
      # e.g. NC_045512.2
      # is the NCBI accession ID for COVID-19 RefSeq complete genome
      first_token = fasta_file_handle.readline().strip().split()[0]      
      if first_token.startswith(">"):
         first_token = first_token[1:]
      else:
         print("Error: FASTA file first line of reference genome does not begin with '>'")
      cur_id = first_token
      ncbi_accession_ids.append(cur_id)

  # return the NCBI accession IDs of the reference genomes
  # in the database directory (which is by default "genomes-of-common-contaminants")
  return ncbi_accession_ids


def get_taxonomy_ids() -> List[int]:
  """
  Given a list of NCBI accession IDs (from the get_fasta_ncbi_accession_ids() method)
  return the corresponding taxonomy IDs of those genomes.
  """
  # get the taxonomy id for this FASTA file


def get_taxonomy_ids_from_file(taxonomy_ids_file) -> List[int]:
  """
  Load a custom list of taxonomy IDs from a file.

  This removes the need to search the large file
  nucl_wgs.accession2taxid.gz (~37 GB uncompressed)
  for just the 20 genomes that we want to classify contaminants in.
  """
  pass
   
   

def main():
  print(get_fasta_ncbi_accession_ids("genomes-of-common-contaminants"))
  # build_parent_map()

if __name__ == "__main__":
  main()