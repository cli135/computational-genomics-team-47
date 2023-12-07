def get_kmer_hit_counts_with_database_from_psuedoreads(pseudoreads, kmer_to_lca, kmer_length):
  """
  Scan through all kmers in the pseudoreads and find which kmers
  in the reads hit (match exactly with) a kmer in the contaminant database,
  and count how many times these matches occur.

  :param pseudoreads: A string representing the pseudoreads.
  :param kmer_to_lca: Dictionary mapping k-mers to their LCA taxonomy IDs.
  :param kmer_length: The length of each k-mer.
  :return: Dictionary of hit counts, mapping taxonomy IDs to counts.
  """
  # Function to split the pseudoreads into k-mers
  def split_into_kmers(sequence, k):
      return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

  # Split the pseudoreads into k-mers
  kmers = split_into_kmers(pseudoreads, kmer_length)

  # Initialize a dictionary to count hits for each taxonomy ID
  hit_counts = {}

  # Iterate over each k-mer in the pseudoreads
  for kmer in kmers:
      # Check if the k-mer is in the contaminant database
      if kmer in kmer_to_lca:          
          # TODO: Maybe fix the 0th index here in case of ties
          # Get the LCA taxonomy ID for this k-mer
          #   lca_node_taxonomy_id = kmer_to_lca[kmer][0]
          if type(kmer_to_lca[kmer]) != str:
              continue
          lca_node_taxonomy_id = kmer_to_lca[kmer]

          # Increment the hit count for this taxonomy ID
          if lca_node_taxonomy_id in hit_counts:
              hit_counts[lca_node_taxonomy_id] += 1
          else:
              hit_counts[lca_node_taxonomy_id] = 1

  return hit_counts

