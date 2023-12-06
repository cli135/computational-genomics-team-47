import os
import sys
import argparse
from typing import List, Dict

# helper files
import taxonomy_tree
import kmer_to_lca_mapping
import get_kmer_hit_counts
import root_to_leaf_paths

def parse_args():
  parse = argparse.ArgumentParser(
    description="Find contaminants in a input query genome by looking up kmers in a known contaminant database"
  )

  parse.add_argument(
    "--db",
    default="genomes-of-common-contaminants",
    help="Name of the directory containing the database of known contaminants \
      which you want to cross-check the input sequence for (default: genomes-of-common-contaminants)"
  )

  parse.add_argument(
    "--taxonomy",
    default="taxonomy",
    help="Name of the directory containing the taxonomy (including names.dmp and nodes.dmp) \
      (default: taxonomy)"
  )

  parse.add_argument(
    "--taxonomy-ids",
    default="taxonomy/custom_taxonomy_ids.txt",
    help="Name of the plaintext file (default: taxonomy/custom_taxonomy_ids.txt) which contains taxonomy ids sourced \
      from NCBI corresponding to the NCBI accession IDs of the FASTA files in \
      the database"
  )

  parse.add_argument(
    "--input-query",
    required=True, # this argument is required since we need to know which sequence to search for contaminants in
    help="Filename of the sequence in which the program will search for contaminants (required)"
  )

  parse.add_argument(
    "--output",
    default=sys.stdout,
    help="Output stream to output results of contamination detection and classification to \
        (default: sys.stdout)"
  )

  return parse.parse_args()

def main():

  # parse command line arguments
  args = parse_args()

  # print out command line arguments entered
  print("Database:", args.db)
  print("Input query sequence:", args.input_query)
  print("Output:", args.output)

  # ====================================================
  # Begin the contamination detection and classification 
  # ====================================================

  # Step 1. Build the taxonomy
  # This method is found in the taxonomy_tree.py file
  pruned_taxonomy_id_to_node, pruned_taxonomy_id_to_parent_id, pruned_tree_root_node = \
    taxonomy_tree.build_parent_map(
      taxonomy_directory=args.taxonomy,
      custom_taxonomy_ids_filename=args.taxonomy_ids
    )

  # Step 2. After the parent map (i.e. taxonomy tree) is built in taxonomy_tree.py,
  # We will build the database with actual cross-references to kmers and lcas
  # This method is found in the kmer_to_lca_mapping.py file
  kmer_to_lca : Dict[str, int] = \
    kmer_to_lca_mapping.build_database(args.db, args.taxonomy_ids)

  # Step 3. Make the pseudoreads from the query sequence
  # pseudoreads = make_pseudoreads(query_sequence)

  # Step 4. Scan through the query pseudoreads and count how many times each k-mer
  # is hit (matches exactly) with a kmer in the database of contaminants.
  # This method is found in the get_kmer_hit_counts.py file
  hit_counts = get_kmer_hit_counts.get_kmer_hit_counts_with_database_from_psuedoreads()

  # Step 5. Root to leaf paths to find the most likely contaminant
  # This method is found in the root_to_leaf_paths.py file
  path_of_likely_contaminant = root_to_leaf_paths.find_likely_contaminant()

  # Step 6. print data and summary below of contaminants found
  print_summary_contaminants_found()

  # The program has finished at this point
  exit(0)


# Step 6. print data and summary below of contaminants found
# TODO need to add necessary arguments and parameters
def print_summary_contaminants_found():
  """
  Prints to stdout a summary of the contamination information
  found, after the pseudoreads have been checked for in the database 
  """
  # TODO implement me!
  raise NotImplementedError()


if __name__ == "__main__":
  main()