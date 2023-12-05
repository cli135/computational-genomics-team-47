# Step 3. Scan through the query pseudoreads and count how many times each k-mer
# is hit (matches exactly) with a kmer in the database of contaminants.

def get_kmer_hit_counts_with_database_from_psuedoreads():
  """
  Scanning through all kmers in the pseudoreads and finding which kmers
  in the reads hit (match exactly with) a kmer in the contaminant database,
  and counting how many times these matches occur.

  This involves:

  for all kmers
    if kmer is found in database library
      get the lca corresponding to that kmer and add it to a hit_counts map
  compute the highest weighted root_to_leaf path, given the hit_counts map

  An identified contaminant is the species that has the highest weighted root to leaf path.  

  @param: kmer_to_lca dictionary

  @return: the hit_counts dictionary, which is a key-value mapping of taxonomy ids (ints) to counts (also ints)
    of how many times they were 'hit' by an exact match with one of the kmers in the pseudoreads.
  """
  # psuedocode below (feel free to change anything):
  # call the pseudoreads function to split up the query genome into pseudoreads, probably FASTQ format
 
  # for all kmers in the psuedoreads:
    
    # if kmer in kmer_to_lca (i.e. if the kmer is found in the precomputed library database
    #   based on the the ~20 reference genomes that are common contaminants)
      
      # then we found an exact match hit with a contaminant! To remember which contaminant,
      # we add 1 count to that node corresponding to that taxonomy id in the tree, i.e.

      # lca_node_taxonomy_id = kmer_to_lca[kmer] # getting the taxonomy id of the least common ancestor node
      # hit_counts[lca_node_taxonomy_id] += 1

  # return the hit_counts array which is the goal of computing this method.
  # This hit_counts array will be fed into Step 4. find_highest_weighted_root_to_leaf_path()

