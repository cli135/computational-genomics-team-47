import os
import taxonomy_tree
from taxonomy_tree import TaxaTree
from typing import Dict

"""
kmer_to_lca_mapping.py

Step 2.)

Takes the taxonomy tree (taxonomy_id_to_node, taxonomy_id_to_parent_id, and root_node
from Step 1.) taxonomy_tree.py) as the input and creates a dictionary mapping, kmer_to_lca_mapping, 
that we can later use to lookup the least common ancestor taxonomic node of all the taxonomic nodes
that a particular given k-mer appears in.

Notably, this step involves actually going over the ~20 reference genomes of common contaminants
(i.e. Mycoplasma, E.coli lambda page phiX174, etc.) and doing offline-preprocessing over all kmers
in the reference genomes there.

The genomes total in about ~90 MB in file size so this should be manageable in terms of memory and speed.

References
----------
This code is inspired from Jennifer Lu's code linked below
https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py

How to run
----------

$ python src/kmer_to_lca_mapping.py

The above line is for if you want to call and test this script directly.
This currently has no functionality (as it is not necessary to call it directly for now)
but might be implemented later if we need to do this.

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

# Step 2. After the parent map (i.e. taxonomy tree) is built in taxonomy_tree.py,
# We will build the database with actual cross-references to kmers and lcas
def build_database() -> Dict[str, int]:
  """
  Given a directory of ~20 FASTA files with genomes of common contaminants,
  we want to traverse all kmers in the FASTA files.
  
  For each kmer update the kmer_to_lca dictionary, we want to set

    kmer_to_lca[kmer] = lca(kmer_to_lca[kmer], curr_taxonomy_id)
  
  the least commmon ancestor of the current value in kmer_to_lca dictionary
  and the current taxonomy id of the node that the kmer was just found in.

  This is like accumulating the LCA iteratively as we go across all the kmers.

  @return: the kmer_to_lca_mapping dictionary
  """
  # TODO partial code below but unfinished

  # for each FASTA file
  for f in os.listdir("./"):
    # get the NCBI accession id for this FASTA file
    pass
    # get the taxonomy id for this FASTA file
    pass
    # for each distinct k-mer
      # update the LCA dictionary with the LCA of the existing value and this value
  pass

  # TODO implement me!
  raise NotImplementedError()



def lca(taxonomy_id_to_parent_id, first_taxonomy_id : int, second_taxonomy_id : int) -> int:
  """
  Compute the least common ancestor of the nodes in the tree
  with taxonomy id first_taxonomy_id and second_taxonomy_id.

  Some notes / invariants about this function:

    lca(map, a, b) = lca(map, b, a) # symmetry of a and b
    lca(map, 1, a) = lca(map, a, 1) = 1 # anything lca'ed with the root is root
    lca(map, 0, a) = lca(map, a, 0) = a # 0 is not a valid taxonomy_id,
      so this means to return the other taxonomy_id unchanged

  This code is inspired from Derrick Wood's krakenutil.cpp code at:
  https://github.com/DerrickWood/kraken/blob/master/src/krakenutil.cpp
  """
  # This method deals with integers only
  # So the primary use will be the dictionary taxonomy_id_to_parent_id
  # from the taxonomy_tree file

  # TODO implement me!
  raise NotImplementedError()




def lca(root_node : TaxaTree, first_node : TaxaTree, second_node : TaxaTree) -> TaxaTree:
  """
  Compute the least common ancestor of two nodes, first_node
  and second_node, in the TaxaTree rooted at root_node

  ***Same functionality as the other lca(a, b) method above,
  but the type signature is overloaded this time to deal with TaxaTree
  nodes directly instead of cycling through integer dictionaries,
  i.e. this time the n-ary tree structure is explicitly stored in
  memory instead of chasing integer references***

  Some notes / invariants about this function:

    lca(root_node, a_node, b_node) = lca(root_node, b_node, a_node) # symmetry of a_node and b_node
    lca(root_node, root_node, a) = lca(root_node, a, root_node) = 1 # anything lca'ed with the root is root
    lca(root_node, 0, a) = lca(root_node, a, 0) = a # 0 is not a valid taxonomy_id,
      so this means to return the other taxonomy_id unchanged

  This code is inspired from Derrick Wood's krakenutil.cpp code at:
  https://github.com/DerrickWood/kraken/blob/master/src/krakenutil.cpp
  """
  # The implementation for this function may be similar to
  # finding the intersection of two linked lists on a leaf-to-root path
  # or any other algorithm for finding the least common ancestor of two nodes
  # in a binary tree

  # TODO implement me!
  raise NotImplementedError()

