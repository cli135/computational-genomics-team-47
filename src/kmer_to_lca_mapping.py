import os
import taxonomy_tree
from taxonomy_tree import TaxaTree
from typing import Dict, List, Set
from collections import defaultdict

"""
kmer_to_lca_mapping.py

Step 2.)

Takes the taxonomy tree (pruned_taxonomy_id_to_node, pruned_taxonomy_id_to_parent_id, and pruned_tree_root_node
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
def build_database(
    file_directory: str,
    custom_taxonomy_ids_filename : str,
    k: int,
    taxonomy_id_to_parent_id : Dict[str, str]) -> Dict[str, str]:
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

  # Dictionary mapping k-mers to their LCA taxonomy IDs
  kmers_to_lca = {}

  # Get the NCBI accession id to tax id mapping
  ncbi_accession_id_to_tax_id_mapping = \
    make_ncbi_accession_id_to_tax_id_mapping(
      custom_taxonomy_ids_filename
    )
  
  # File that we are searching through
  file_count = 1

  # For each FASTA file
  for f in os.listdir(f"./{file_directory}"):

    # opening the file
    file_path = os.path.join(file_directory, f)

    with open(file_path, "r") as reference_genome_assembly:
      # reading in the first line of the file
      first_line = reference_genome_assembly.readline()
        
      # get the accession id for this FASTA file
      accession_id = first_line.split()[0][1:]

      # check if the accession id is in the accesion id to tax id mapping
      if (accession_id in ncbi_accession_id_to_tax_id_mapping.keys()):
        # We are searching this file
        print(f"searching file {file_count} with accession_id {str(accession_id)}")

        # if it is, then get the tax id
        tax_id = ncbi_accession_id_to_tax_id_mapping[accession_id]
        
        # Split the genome assembly into kmers
        reference_genome_assembly_sequence = ""
        for line in reference_genome_assembly:
            if line.startswith(">"):
                # Skip header lines in FASTA format
                continue
            # Add this line to the genome assembly sequence
            reference_genome_assembly_sequence += line.strip()

        for i in range(len(reference_genome_assembly_sequence) - k + 1):
            kmer = reference_genome_assembly_sequence[i:i + k]

            # If it is a kmer we haven't seen before, then set it to the tax_id corresponding
            # to the accession id of this FASTA file
            if (kmer not in kmers_to_lca.keys()):
              # kmers_to_lca[kmer] = [tax_id]
              kmers_to_lca[kmer] = tax_id
            else:
              # If it is a kmer we have seen before, then update the LCA of this kmer
              kmers_to_lca[kmer] = lca(taxonomy_id_to_parent_id, kmers_to_lca[kmer], tax_id)
      else:
        continue
  
  return kmers_to_lca

def make_ncbi_accession_id_to_tax_id_mapping(
    custom_taxonomy_ids_filename : str) -> List[str]:
  """
  Load a custom taxonomy IDs file and 

  This removes the need to search the large file
  nucl_wgs.accession2taxid.gz (~37 GB uncompressed)
  for just the 20 genomes that we want to classify contaminants in.
  """
  ncbi_accession_id_to_tax_id = {}
  # using the custom seqid2taxid instead of searching the whole 37 GB mapping at
  # https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
  with open(custom_taxonomy_ids_filename, 'r') as fp:
    # getting the second column of each line,
    # which is the taxonomy ID of that accession ID
    for line in fp.readlines():
      tokens = line.strip().split()
      ncbi_accession_id = tokens[0]
      tax_id = tokens[1]
      ncbi_accession_id_to_tax_id[ncbi_accession_id] = tax_id
  # return the dictionary
  return ncbi_accession_id_to_tax_id

def lca(taxonomy_id_to_parent_id, first_taxonomy_id : str, second_taxonomy_id : str) -> str:
  """
  Compute the least common ancestor of the nodes in the tree
  with taxonomy ids first_taxonomy_id and second_taxonomy_id.

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
  # raise NotImplementedError()
  a = first_taxonomy_id
  b = second_taxonomy_id
  # it seems we just return the other if one is '0'
  # and '0' is not the root taxonomy id, '1' is the root
  if a == '0' or b == '0':
    return a if a else b
  
  path_from_a_to_root : Set[str] = set()
  # collect taxonomy ids of the ndoes
  # on the path from a towards root in a set
  while a != '1':
    path_from_a_to_root.add(a)
    a = taxonomy_id_to_parent_id[a]
  
  # ok, now we track up from b leaf node
  # and find the first point of intersection
  # using the hashset built above
  while b not in path_from_a_to_root:
    b = taxonomy_id_to_parent_id[b]
  
  # at this point, b is the first node found
  # that is also on the path from a to root.
  # this means we found the lca, the earliest
  # point of intersection of the path from a to root
  # and the path from b to root
  return b

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
  # raise NotImplementedError()


# def main():
#   output = build_database("gocc-shortened", "taxonomy/custom_taxonomy_ids.txt", 31)
#   print(output)
#   lca(None, None, None)


# if __name__ == "__main__":
#   main()