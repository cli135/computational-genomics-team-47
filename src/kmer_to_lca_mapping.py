import taxonomy_tree
from taxonomy_tree import TaxaTree

# Step 2. After the parent map (i.e. taxonomy tree) is built in taxonomy_tree.py,
# We will build the database with actual cross-references to kmers and lcas
def build_database():
  """
  Given a directory of ~20 FASTA files with genomes of common contaminants,
  we want to traverse all kmers in the FASTA files.
  
  For each kmer update the kmer_to_lca dictionary, we want to set

    kmer_to_lca[kmer] = lca(kmer_to_lca[kmer], curr_taxonomy_id)
  
  the least commmon ancestor of the current value in kmer_to_lca dictionary
  and the current taxonomy id of the node that the kmer was just found in.

  This is like accumulating the LCA iteratively as we go across all the kmers.
  """
  # for each FASTA file
  for f in os.listdir("./"):
    # get the NCBI accession id for this FASTA file
    pass
    # get the taxonomy id for this FASTA file
    pass
    # for each distinct k-mer
      # update the LCA dictionary with the LCA of the existing value and this value
  pass


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
  pass
  # This method deals with integers only
  # So the primary use will be the dictionary taxonomy_id_to_parent_id
  # from the taxonomy_tree file


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