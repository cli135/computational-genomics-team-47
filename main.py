import os


# Step 1. Build the taxonomy
#   this method is found in the taxonomy_tree.py file


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
  
  Keyword arguments:
  argument -- description
  Return: return_description
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


def main():
  pass

if __name__ == "__main__":
  main()