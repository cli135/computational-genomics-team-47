import os, sys, argparse
from time import gmtime
from time import strftime 

"""
Takes the nodes.dmp, names.dmp, and the accession2taxid as the input and creates a taxonomy tree
that we can later use to determine the relatives of a certain taxid.

References
----------
Code inspired from Jen Lu's code linked below
https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py
"""

# TaxaTree data structure representing nodes in taxonomy
class TaxaTree:
    def __init__(self, tax_id, parent = None, children = [], name = "", rank):
        # Node's own taxid, name, and rank
        self.tax_id = tax_id
        self.name = name
        self.rank = rank

        # Setting Parent
        self.parent = parent
        self.parentTaxId = -1

        # Children
        self.children = children
        self.childrenTaxId = []


    def add_child(self, child_node):
        # Two-way connection between nodes!
        child_node.parent = self
        self.children.append(child_node)

    def get_parent(self):
        return self.parent

    def get_children(self):
        return self.children

    def get_tax_id(self):
        return self.tax_id

# for utility, for convenience
map_ranks_to_abbreviations = {
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S"
}

# initializing the main dictionaries which let us navigate
# and traverse the taxonomic tree

# tax id to node (where node is a TaxaTree object)
# https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py
taxonomy_id_to_node = {}

# a tax id mapped to its parent's tax id
# https://github.com/DerrickWood/kraken/blob/master/src/krakenutil.cpp
taxonomy_id_to_parent_id = {}

# redundancy in the dictionaries let us navigate both up and down the tree
# following forward pointers and backpointers in the n-ary tree like a deque

def build_parent_map():
  """
  This method closely implements
  https://github.com/DerrickWood/kraken/blob/master/src/krakenutil.cpp
  
  Which basically extracts the first two columns from nodes.dmp and puts them into
  the taxonomy_id_to_parent_id
  
  Keyword arguments:
  argument -- description
  Return: return_description
  """
  
  file_handle = open('./taxonomy/nodes.dmp', 'r')
  # finish soon
  


