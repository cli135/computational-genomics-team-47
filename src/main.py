import os

# helper files
import taxonomy_tree
import kmer_to_lca_mapping
import get_kmer_hit_counts
import root_to_leaf_paths


def main():
  # Step 1. Build the taxonomy
  # This method is found in the taxonomy_tree.py file
  taxonomy_tree.build_parent_map()

  # Step 2. After the parent map (i.e. taxonomy tree) is built in taxonomy_tree.py,
  # We will build the database with actual cross-references to kmers and lcas
  # This method is found in the kmer_to_lca_mapping.py file
  kmer_to_lca_mapping.build_database()

  # Step 3. Scan through the query pseudoreads and count how many times each k-mer
  # is hit (matches exactly) with a kmer in the database of contaminants.
  # This method is found in the get_kmer_hit_counts.py file
  hit_counts = get_kmer_hit_counts.get_kmer_hit_counts_with_database_from_psuedoreads()

  # Step 4. Root to leaf paths to find the most likely contaminant
  # This method is found in the root_to_leaf_paths.py file
  path_of_likely_contaminant = root_to_leaf_paths.find_likely_contaminant()

  # Step 5. print data and summary below of contaminants found
  print_summary_contaminants_found()

def print_summary_contaminants_found():
  """
  Prints to stdout a summary of the contamination information
  found, after the pseudoreads have been checked for in the database 
  """
  pass
  

if __name__ == "__main__":
  main()