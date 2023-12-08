import os
import sys
import argparse
from typing import List, Dict
import time

# helper files
import taxonomy_tree
import kmer_to_lca_mapping
import get_kmer_hit_counts
import pseudoreads

# Command line option parsing
def parse_args():
  parse = argparse.ArgumentParser(
    description="Find contaminants in a input query genome by looking up kmers in a known contaminant database"
  )

  parse.add_argument(
    "--db",
    # TODO swap out which default database we want to use
    default="genomes-of-common-contaminants-size-halved",

    # ------------------------
    # FOR TUTORIAL PURPOSES:
    # Use: default="gocc-shortened"
    # 
    # ***NOTE:*** Using default="genomes-of-common-contaminants" as the reference database
    #             requires significant computing power (more than typical ugrad cluster). So we
    #            have provided a shortened version of the database for to use for a tutorial.
    # ------------------------

    help="Name of the directory containing the database of known contaminants \
      which you want to cross-check the input sequence for (default: genomes-of-common-contaminants)"
  )

    # ------------------------
    # PLEASE ENSURE THAT "names-and-nodes.zip" in the "taxonomy" directory
    # is unzipped before running the program
    # ------------------------

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

    # ------------------------
    # THIS IS THE ONLY REQUIRED INPUT VIA COMMAND LINE ARGUMENTS
    # Here, the genome assembly to search for contamination in is specified
    # ------------------------

  parse.add_argument(
    "--input-query",
    required=True, # this argument is required since we need to know which sequence to search for contaminants in
    help="Filename of the sequence in which the program will search for contaminants (required)"
  )

  parse.add_argument(
    "--k",
    default=31,
    type=int,
    help="k, the length of the kmer (default: k = 12, which runs on the ugrad machines well (within memory constraints). k = 31 is ideal if the computer has enough memory.)"
  )

  # parse.add_argument(
  #   "--output",
  #   default=sys.stdout,
  #   help="Output stream to output results of contamination detection and classification to \
  #       (default: sys.stdout)"
  # )

  return parse.parse_args()

def main():

  # parse command line arguments
  args = parse_args()

  # print out command line arguments entered
  print("Database:", args.db)
  print("Input query sequence:", args.input_query)
  print("Taxonomy:", args.taxonomy)
  print("Seq ID to Taxonomy ID Mapping:", args.taxonomy_ids)
  # print("Output:", args.output)

  # ====================================================
  # Begin the contamination detection and classification 
  # ====================================================

  start_time = time.time()

  # Step 0. Pick k
  # the kmer length
  k = args.k
  print("kmer length, k:", k)
  print("k, the length of the kmer (default: k = 12, which runs on the ugrad machines well (within memory constraints). k = 31 is ideal if the computer has enough memory.)")

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
      k,
      pruned_taxonomy_id_to_parent_id
    )

  # Step 3. Make the pseudoreads from the query sequence
  pseudoreads_list = pseudoreads.split_genome_into_pseudo_reads_from_fasta(args.input_query)
  
  # Step 4. Scan through the query pseudoreads and count how many times each k-mer
  # is hit (matches exactly) with a kmer in the database of contaminants.
  # This method is found in the get_kmer_hit_counts.py file
  pseudoread_to_hit_counts = {}
  # total_accumulated_hit_counts = {}
  for pseudoread in pseudoreads_list:
    # Feed each psuedoread to the function to get the hit counts
    hit_counts = get_kmer_hit_counts.get_kmer_hit_counts_with_database_from_psuedoreads(pseudoread, kmer_to_lca, k)
    pseudoread_to_hit_counts[pseudoread] = hit_counts

  # Step 5. print data and summary below of contaminants found

  # Dictionary of taxonomy ids to assembly name
  genome_data = {
      '511145': "Escherichia coli str. K-12 substr. MG1655, complete genome",
      '208964': "Pseudomonas aeruginosa PAO1, complete genome 6,264,404 bp circular DNA",
      '198214': "Shigella flexneri 2a str. 301 chromosome, complete genome 4,607,202 bp circular DNA",
      '99287': "Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome 4,857,450 bp circular DNA",
      '386585': "Escherichia coli O157:H7 str. Sakai DNA, complete genome",
      '224308': "Bacillus subtilis subsp. subtilis str. 168 complete genome",
      '192222': "Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819 chromosome, complete genome",
      '227882': "Streptomyces avermitilis MA-4680 = NBRC 14893, complete sequence",
      '340047': "Mycoplasma capricolum subsp. capricolum ATCC 27343, complete sequence",
      '93061': "Staphylococcus aureus subsp. aureus NCTC 8325 chromosome, complete genome",
      '871585': "Acinetobacter pittii PHEA-2 chromosome, complete genome",
      '83332': "Mycobacterium tuberculosis H37Rv, complete genome",
      '1125630': "Klebsiella pneumoniae subsp. pneumoniae HS11286 chromosome, complete genome",
      '2886930': "Escherichia phage phiX174, complete genome",
      '32604': "Human herpesvirus 6B, complete genome",
      '60550': "Burkholderia pyrrocinia strain DSM 10685 chromosome 1, complete sequence",
      '10376': "Human gammaherpesvirus 4, complete genome",
      '28449': "Neisseria subflava strain ATCC 49275 chromosome, complete genome",
      '858423': "Bradyrhizobium arachidis strain CCBAU 051107 chromosome, complete genome",
      '735': "Haemophilus parahaemolyticus strain FDAARGOS_1199 chromosome, complete genome",
      '2842456': "Ralstonia wenshanensis strain 56D2 chromosome, complete genome",
      '38310': "Rhodococcus coprophilus strain NCTC10994 chromosome 1, complete sequence",
      '655813': "Streptococcus oralis ATCC 35037 strain NCTC 11427 chromosome 1, complete sequence",
      '2697049': "Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome"
  }
  
  print("#############################################")  
  print("############## SUMMARY ######################")
  print("#############################################")  
  print()
  print("############## SEQUENCE CLASSIFICATION ######")
  print("#############################################")  
  print_pseudoreads_classified(pseudoread_to_hit_counts, genome_data=genome_data)
  print()
  print("###### KMER MATCHES FOR CLASSIFICATION ######")
  print("#############################################")  
  print_kmers_classified(pseudoread_to_hit_counts, genome_data=genome_data)
  print()

  end_time = time.time()
  print("############## TIME TAKEN ###################")
  print(f"Total time taken: {end_time - start_time} seconds")

  # The program has finished at this point
  exit(0)


# Step 6. print data and summary below of contaminants found
def print_pseudoreads_classified(pseudoread_to_hit_counts, genome_data):
  """
  Prints to stdout a summary of what percentage of pseudoreads
  were mapped to each contaminant
  """
  tax_count = {}
  for pseudoread, hit_counts in pseudoread_to_hit_counts.items():
      # Find what the pseudoread mapped to the most
      # add a bounds check to avoid an error in
      # the case where hit_counts is empty
      if len(hit_counts) != 0:
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
    if taxonomy_id in genome_data:
      print(f"{round(count/total_hit_count*100, 2)}% of reads mapped to Taxonomy ID {taxonomy_id}, {genome_data[taxonomy_id]}")
    else:
      print(f"{round(count/total_hit_count*100, 2)}% of reads mapped to Taxonomy ID {taxonomy_id}")
      
  print()
  for tax_id, count in tax_count.items():
    if taxonomy_id in genome_data:
      print(f"Tax ID: {tax_id}, {genome_data[taxonomy_id]}, Number of Pseudoreads: {count}")
    else:
      print(f"Tax ID: {tax_id}, Number of Pseudoreads: {count}")

  print()

# Step 6. print data and summary below of kmer hits
def print_kmers_classified(reads_dict, genome_data):
  """
  Prints to stdout a summary of how many kmers from the pseudoreads
  were mapped to each taxid inputted into the database
  """
  tax_count = {}
  overal = 0
  for tax_info in reads_dict.values():
      for tax_id, count in tax_info.items():
          if tax_id in tax_count:
              tax_count[tax_id] += count
              overal += count
          else:
              tax_count[tax_id] = count
              overal += count
  
  for tax_id, total_count in tax_count.items():
    if tax_id in genome_data:
      print(f"Tax ID: {tax_id}, {genome_data[tax_id]}\n")
    else:
      print(f"Tax ID: {tax_id}\n")

    print(f"Total Accumulated Hit Counts: {total_count} out of {overal}\n")

if __name__ == "__main__":
  main()