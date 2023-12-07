import os
import sys
import argparse
from typing import List, Dict

# helper files
import taxonomy_tree
import kmer_to_lca_mapping
import get_kmer_hit_counts
import root_to_leaf_paths
import pseudoreads

def parse_args():
  parse = argparse.ArgumentParser(
    description="Find contaminants in a input query genome by looking up kmers in a known contaminant database"
  )

  parse.add_argument(
    "--db",
    # TODO swap out which default database we want to use
    default="genomes-of-common-contaminants",
    # we'll use the shorter one
    # default="gocc-shortened",
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
  # Dict[str, str]
  kmer_to_lca = \
    kmer_to_lca_mapping.build_database(
      args.db,
      args.taxonomy_ids,
      31,
      pruned_taxonomy_id_to_parent_id
    )
  
  # Step 3. Make the pseudoreads from the query sequence
  # TODO I think the below syntax might need fixing since I'm not sure if I did it right
  pseudoreads_list = pseudoreads.split_genome_into_pseudo_reads_from_fasta("covid-assemblies/covid_assembly_fasta.txt")
  
  # Step 4. Scan through the query pseudoreads and count how many times each k-mer
  # is hit (matches exactly) with a kmer in the database of contaminants.
  # This method is found in the get_kmer_hit_counts.py file
  pseudoread_to_hit_counts = {}
  # total_accumulated_hit_counts = {}
  for pseudoread in pseudoreads_list:
    # Feed each psuedoread to the function to get the hit counts
    hit_counts = get_kmer_hit_counts.get_kmer_hit_counts_with_database_from_psuedoreads(pseudoread, kmer_to_lca, 31)
    pseudoread_to_hit_counts[pseudoread] = hit_counts
    # .update() will add all entries of one Python dictionary to another
    # total_accumulated_hit_counts.update(hit_counts)
  print(pseudoread_to_hit_counts)
  # print(total_accumulated_hit_counts)

  # Step 5. Root to leaf paths to find the most likely contaminant
  # This method is found in the root_to_leaf_paths.py file
  # path_of_likely_contaminant = root_to_leaf_paths.find_likely_contaminant()

  # Step 6. print data and summary below of contaminants found
  print_pseudoreads_classified(pseudoread_to_hit_counts)
  print_kmers_classified(pseudoread_to_hit_counts)

  # The program has finished at this point
  exit(0)


# Step 6. print data and summary below of contaminants found
def print_pseudoreads_classified(pseudoread_to_hit_counts):
  """
  Prints to stdout a summary of what percentage of pseudoreads
  were mapped to each contaminant
  """
  tax_count = {}
  for pseudoread, hit_counts in pseudoread_to_hit_counts.items():
      # Find what the pseudoread mapped to the most
      taxonomy_id_with_max_hits = \
        max(hit_counts.keys(), key=lambda x : hit_counts[x])
      # Keep track of what taxonomy ids have been mapped to so far and how often
      if taxonomy_id_with_max_hits in tax_count:
          tax_count[taxonomy_id_with_max_hits] += 1
      else:
          tax_count[taxonomy_id_with_max_hits] = 1
  print("Psuedoreads classified:")
  total_hit_count = sum(tax_count.values())
  for taxonomy_id, count in tax_count.items():
      print(f"{round(count/total_hit_count*100, 2)}% of reads mapped to Taxonomy ID {taxonomy_id}")
      
  
  for tax_id, count in tax_count.items():
      print(f"Tax ID: {tax_id}, Number of Pseudoreads: {count}")


# Step 6. print data and summary below of kmer hits
def print_kmers_classified(reads_dict):
  """
  Prints to stdout a summary of how many kmers from the pseudoreads
  were mapped to each taxid inputted into the database
  """
  tax_count = {}
  
  for tax_info in reads_dict.values():
      for tax_id, count in tax_info.items():
          if tax_id in tax_count:
              tax_count[tax_id] += count
          else:
              tax_count[tax_id] = count
  
  for tax_id, total_count in tax_count.items():
      print(f"Tax ID: {tax_id}, Total Accumulated Hit Counts: {total_count}")


if __name__ == "__main__":
  main()