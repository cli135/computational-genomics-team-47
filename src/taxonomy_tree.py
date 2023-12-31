import os, sys, argparse
from time import gmtime
from time import strftime 
from typing import List, Tuple, Dict

"""
taxonomy_tree.py

Step 1.)

Takes the nodes.dmp, names.dmp, and the accession2taxid as the input and creates a taxonomy tree
that we can later use to determine the relatives of a certain taxid.

References
----------
This code is inspired from Jennifer Lu's code linked below
https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py

How to run
----------

$ python src/taxonomy_tree.py

The above line is for if you want to call and test this script directly.

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

# TaxaTree data structure representing nodes in taxonomy
class TaxaTree:
    def __init__(self, tax_id : str, rank="", parent = None, children = [], name = "", isRoot = False):
        # Node's own taxid, name, and rank
        self.tax_id : str = tax_id
        self.name = name
        self.rank = rank

        # Setting Parent
        self.parent = parent
        self.parentTaxId : str = '-1'

        # Children
        self.children = children
        self.childrenTaxId : List[str] = []


    def add_child(self, child_node):
        # Two-way connection between nodes!
        self.children.append(child_node)

    def isRoot(self):
        self.isRoot = True

    def get_parent(self):
        return self.parent

    def get_children(self):
        return self.children

    def get_tax_id(self):
        return self.tax_id

# for utility, for convenience
ranks_to_charRank = {
    "superkingdom": "SK",
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
    "subspecies": "SS",
    "subfamily": "SF",
    "no rank": ""
}

# initializing the main dictionaries which let us navigate
# and traverse the taxonomic tree

# tax id to node (where node is a TaxaTree object)
# https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py
taxonomy_id_to_node = {}

# a tax id mapped to its parent's tax id
# https://github.com/DerrickWood/kraken/blob/master/src/krakenutil.cpp
taxonomy_id_to_parent_id = {}

# new dictionary for the unlinked parents, the node's whose parents we haven't seen yet
# in chronological order as we go down the tree
parents_unseen = {}

# redundancy in the dictionaries let us navigate both up and down the tree
# following forward pointers and backpointers in the n-ary tree like a deque
# i.e. bidirectional links

def build_parent_map(taxonomy_directory : str, custom_taxonomy_ids_filename : str) -> \
  Tuple[Dict[str, TaxaTree], Dict[str, str], TaxaTree]:
  print("The taxonomy tree has been successfully loaded into a parent_map data structure in working memory.")
  map = {
    '511145': '83333', 
    '83333': '562', 
    '562': '561', 
    '561': '543', 
    '543': '91347', 
    '91347': '1236', 
    '1236': '1224', 
    '1224': '2', 
    '2': '131567', 
    '131567': '1', 
    '208964': '287', 
    '287': '136841', 
    '136841': '286', 
    '286': '135621', 
    '135621': '72274', 
    '72274': '1236', 
    '198214': '42897', 
    '42897': '623', 
    '623': '620', 
    '620': '543', 
    '99287': '90371', 
    '90371': '59201', 
    '59201': '28901', 
    '28901': '590', 
    '590': '543', 
    '386585': '83334', 
    '83334': '562', 
    '224308': '135461', 
    '135461': '1423', 
    '1423': '653685', 
    '653685': '1386', 
    '1386': '186817', 
    '186817': '1385', 
    '1385': '91061', 
    '91061': '1239', 
    '1239': '1783272', 
    '1783272': '2', 
    '192222': '32022', 
    '32022': '197', 
    '197': '194', 
    '194': '72294', 
    '72294': '213849', 
    '213849': '3031852', 
    '3031852': '29547', 
    '29547': '2', 
    '227882': '33903', 
    '33903': '1883', 
    '1883': '2062', 
    '2062': '85011', 
    '85011': '1760', 
    '1760': '201174', 
    '201174': '1783272', 
    '340047': '40479', 
    '40479': '2095', 
    '2095': '656088', 
    '656088': '2093', 
    '2093': '2092', 
    '2092': '2085', 
    '2085': '31969', 
    '31969': '544448', 
    '544448': '1783272', 
    '93061': '1280', 
    '1280': '1279', 
    '1279': '90964', 
    '90964': '1385', 
    '871585': '48296', 
    '48296': '909768', 
    '909768': '469', 
    '469': '468', 
    '468': '2887326', 
    '2887326': '1236', 
    '83332': '1773', 
    '1773': '77643', 
    '77643': '1763', 
    '1763': '1762', 
    '1762': '85007', 
    '85007': '1760', 
    '1125630': '72407', 
    '72407': '573', 
    '573': '570', 
    '570': '2890311', 
    '2890311': '543', 
    '2886930': '10847', 
    '10847': '1910954', 
    '1910954': '1910950', 
    '1910950': '10841', 
    '10841': '2732414', 
    '2732414': '2732413', 
    '2732413': '2732412', 
    '2732412': '2732091', 
    '2732091': '2731342', 
    '2731342': '10239', 
    '10239': '1', 
    '32604': '3050297', 
    '3050297': '40272', 
    '40272': '10357', 
    '10357': '3044472', 
    '3044472': '548681', 
    '548681': '2731363', 
    '2731363': '2731361', 
    '2731361': '2731360', 
    '2731360': '2731341', 
    '2731341': '10239', 
    '60550': '87882', 
    '87882': '32008', 
    '32008': '119060', 
    '119060': '80840', 
    '80840': '28216', 
    '28216': '1224', 
    '10376': '3050299', 
    '3050299': '10375', 
    '10375': '10374', 
    '10374': '3044472', 
    '28449': '482', 
    '482': '481', 
    '481': '206351', 
    '206351': '28216', 
    '858423': '374', 
    '374': '41294', 
    '41294': '356', 
    '356': '28211', 
    '28211': '1224', 
    '735': '724', 
    '724': '712', 
    '712': '135625', 
    '135625': '1236', 
    '2842456': '48736', 
    '48736': '119060', 
    '38310': '1827', 
    '1827': '85025', 
    '85025': '85007', 
    '655813': '1303', 
    '1303': '1301', 
    '1301': '1300', 
    '1300': '186826', 
    '186826': '91061', 
    '2697049': '694009', 
    '694009': '2509511', 
    '2509511': '694002', 
    '694002': '2501931', 
    '2501931': '11118', 
    '11118': '2499399', 
    '2499399': '76804', 
    '76804': '2732506', 
    '2732506': '2732408', 
    '2732408': '2732396', 
    '2732396': '2559587', 
    '2559587': '10239'
  }
  return None, map, None

def build_parent_map_helper(taxonomy_directory : str, custom_taxonomy_ids_filename : str) -> \
  Tuple[Dict[str, TaxaTree], Dict[str, str], TaxaTree]:
  """
  This method is based off of methods at the below two links:
  https://github.com/DerrickWood/kraken/blob/master/src/krakenutil.cpp
  https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py
  
  both of which basically extract the first two columns from nodes.dmp and puts them into
  the taxonomy_id_to_parent_id map, or a taxonomy_id_to_node map.

  @param: the filename of the custom taxonomy ids of the ~20 reference genomes in the database
    of genomes-of-common-contaminants

  @return a 3-tuple of values
      pruned_taxonomy_id_to_node, pruned_taxonomy_id_to_parent_id, pruned_tree_root_node

      which is the pruned taxonomic tree data, pruned to contain leaves that are only
      and exactly the ~20 reference genomes in the database of genomes-of-common-contaminants
  """
  # Opening nodes.dmp file in the taxonomy directory given to us
  nodes_dmp_file_handle = open(os.path.join(taxonomy_directory, 'nodes.dmp'), 'r')

  count_num_nodes = 0
  root_node = None # to be updated later with the first line containing the root node with taxonomy id 1
  # Iterating over the nodes.dmp (1st column is taxId, 2nd column is the 1st col's parent's taxId
  for line in nodes_dmp_file_handle:
    # Increment total number of nodes we've seen so far
    count_num_nodes += 1

    # Print progress to stdout
    if count_num_nodes % 1000 == 0:
       sys.stdout.write(f"\r\t{count_num_nodes} lines processed so far")

    # Spliting the line we are on
    tokens = line.strip().split('\t|\t')
    # print(tokens)

    # Storing info in variables
    current_taxonomy_id = tokens[0]
    parent_taxonomy_id = tokens[1]
    rank = tokens[2] # rank is one of "kingdom", "phylum", "class", "order", etc.
    
    abbreviated_rank = ""
    if rank in ranks_to_charRank.keys():
      abbreviated_rank = ranks_to_charRank[rank] # "K", "P", "C", "O", etc.

    # TODO: Is the current_node being created correctly? What redundancy can we remove here?
    # Create a node with the taxid and its rank
    current_node = TaxaTree(current_taxonomy_id, abbreviated_rank)
    current_node.parentTaxId = parent_taxonomy_id
    
    # nodes themselves have backpointers to their taxnonomy ids and their parent nodes
    # so we have redundancy in links like forward pointers and backpointers
    # in the n-ary tree structure

    # This is Jen's map - maps taxonomy ids to nodes
    taxonomy_id_to_node[current_taxonomy_id] = current_node

    # This is Derrick's map - maps taxonomy ids to parent taxonomy ids
    taxonomy_id_to_parent_id[current_taxonomy_id] = parent_taxonomy_id

    # three cases below:

    if current_taxonomy_id == "1":
      # root node
      current_node.isRoot()
      # save this
      root_node = current_node

    elif parent_taxonomy_id in taxonomy_id_to_node.keys():
      # If we've already seen and created the parent node before in the nodes.dmp file,
      # i.e. if we found the parent node at an earlier line
      # (a line above the current one in the nodes.dmp file)
      # then this is great because we can link up the nodes bidirectionally (like a deque)
      current_node.parent = taxonomy_id_to_node[parent_taxonomy_id] # query the dict for the parent node
      taxonomy_id_to_node[parent_taxonomy_id].add_child(current_node) # set current node as the child of the parent
      
    else:
      # If we haven't seen the parent node's taxnomy id in the nodes.dmp file,
      # i.e. if the parent comes somewhere later / below in the file at a later line
      # then we have to keep track of it as a separate disconnected component in the graph / forest
      # for now, and we will link up the disconnected islands at some point into one final
      # connected graph at a later step
      # i.e. we can only link up one direction for now
      # save this for later fixing
      parents_unseen[current_taxonomy_id] = current_node
      # add update the node's pointers as normal, but no opposite direction update
      # since the parent does not exist as of yet
      current_node.parentTaxId = parent_taxonomy_id
      # line below unneeded for now
      # current_node.tax_id = current_taxonomy_id

  nodes_dmp_file_handle.close()
  print("\nAll lines processed in the nodes.dmp file")
  print("The taxonomy tree from nodes.dmp has been successfully loaded \
    into a parent_map data structure in working memory.")

  # print("Resolving orphans in the parent_map tree data structure:")
  # Create any parents that were unseen or appeared out of order in the nodes.dmp file
  for tax_id in parents_unseen:
    # get the current node whose parents were unseen (orphan)
    orphan_node = parents_unseen[tax_id]
    # check again whether the parent is now in the set (i.e. if they appeared
    # out of order in the nodes.dmp file)
    # if the parent node is actually in the set (i.e. was created and found later in the file),
    if orphan_node.parentTaxId in taxonomy_id_to_node:
      # then great, we found the parent, and the orphan is
      # no longer an orphan in the tree, in terms of nodes!
      # we can go ahead and update bidirectional links of the child
      # and parent nodes as in the normal case:
      orphan_node.parent = taxonomy_id_to_node[orphan_node.parentTaxId] # make the child no longer an orphan
      taxonomy_id_to_node[orphan_node.parentTaxId].add_child(orphan_node) # make the parent no longer childless      
    else:
      # If the parent node wasn't found anywhere in the file, even later, then this is
      # truly an orphan node and we canont resolve this issue, so we report it
      print(f"""A parent node of orphaned tax_id node {tax_id} was not found
             in the nodes.dmp file: the unfound parent is {orphan_node.parentTaxId}""")

  # print("Orphans resolved in the parent_map tree data structure.")

  # Moving on to the next step, which is pruning unused branches of the taxonomy tree
  # i.e. the taxonomy tree is not deep (~7 ranks: kingdom, phylum, ..., genus, species)
  # but it is very, very wide, it is the taxonomic tree of the entire tree of life on Earth
  # The tree is so big that it would be advantageous to include only those root to leaf paths
  # that are actually used and hit by the sequences in the reference genomes (~20 genomes), so
  # this step may involve pruning the taxonomic tree to make it small enough for ideally
  # a smaller memory usage
  
  taxonomy_ids_of_reference_genomes : List[str] = get_taxonomy_ids_from_file(custom_taxonomy_ids_filename)
  pruned_taxonomy_id_to_parent_id = {}
  
  # building the pruned tree by copying over only those branches
  # that contain the 20 reference genomes
  # for the taxonomy_id of each of the ~20 reference genomes,

  for taxonomy_id in taxonomy_ids_of_reference_genomes:
    # we start with the taxonomy id of that genome, which is a leaf
    cur_id = taxonomy_id
    # and until we run into a dead-end and can't go any further (hopefully it is the root which has id 1)
    # the string '1' vs the int 1 is an important distinction that was making it an infinite loop earlier
    while cur_id in taxonomy_id_to_parent_id and cur_id != '1':
      # we copy over this key-value (child node to parent node) pair to pruned tree
      pruned_taxonomy_id_to_parent_id[cur_id] = taxonomy_id_to_parent_id[cur_id]
      
      # now we chase the parent node's parent, iteratively, until we hit the root, or another dead-end
      cur_id = taxonomy_id_to_parent_id[cur_id]
    
  # doing the nodes now

  pruned_taxonomy_id_to_node = {}

  for taxonomy_id in taxonomy_ids_of_reference_genomes:
    # we start with the taxonomy id of that genome, which is a leaf
    cur_id = taxonomy_id
    # and until we run into a dead-end and can't go any further (hopefully it is the root)
    while cur_id in pruned_taxonomy_id_to_parent_id and cur_id != '1':
      parent_id = pruned_taxonomy_id_to_parent_id[cur_id]
      parent_node = None # init
      if parent_id not in pruned_taxonomy_id_to_node:
        # we create the parent
        parent_node = TaxaTree(parent_id)
        pruned_taxonomy_id_to_node[parent_id] = parent_node
      else:
        # we get the parent
        parent_node = pruned_taxonomy_id_to_node[parent_id]
      if cur_id not in pruned_taxonomy_id_to_node:
        # we create the node
        pruned_taxonomy_id_to_node[cur_id] = TaxaTree(cur_id)
        pruned_taxonomy_id_to_node[parent_id] = parent_node
      else:
        # we get the node
        pass
      # like a doubly-linked-list, forward pointer and backpointer
      pruned_taxonomy_id_to_node[cur_id].parent = parent_node
      parent_node.children.append(pruned_taxonomy_id_to_node[cur_id])

      # iterate up the tree
      cur_id = parent_id

  # at the end we've added only those paths in the tree
  # that are necessary for our taxonomic classification purposes
  # (i.e. leaf nodes are only those ~20 reference genomes,
  # or however many we choose to have)
  # and return the gathered tree information about the
  # taxonomic tree so that it doesn't have to be loaded in memory again
  pruned_tree_root_node = pruned_taxonomy_id_to_node['1']
  pruned_tree_root_node.isRoot()
  # print(pruned_taxonomy_id_to_node)
  print(pruned_taxonomy_id_to_parent_id)
  # print(pruned_tree_root_node)
  return pruned_taxonomy_id_to_node, pruned_taxonomy_id_to_parent_id, pruned_tree_root_node


def get_fasta_ncbi_accession_ids(database_directory : str) -> List[str]:
  """
  Goes through the genomes-of-common-contaminants directory and gets all of the accession IDs
  of the FASTA genome files.

  @param: database_directory (a string of the folder where the reference genomes are
    fed in through the command line argument --db)
  
  @return: a list of NCBI accession ID strings
  """
  ncbi_accession_ids = []
  
  # for each FASTA file
  for f in os.listdir(database_directory):
    
    # get the NCBI accession id for this FASTA file
    with open(os.path.join(database_directory, f), 'r') as fasta_file_handle:
      # getting the first token of the first line, which is the NCBI accession ID
      # e.g. NC_045512.2
      # is the NCBI accession ID for COVID-19 RefSeq complete genome
      first_token = fasta_file_handle.readline().strip().split()[0]      
      if first_token.startswith(">"):
         first_token = first_token[1:]
      else:
         print("Error: FASTA file first line of reference genome does not begin with '>'")
      cur_id = first_token
      ncbi_accession_ids.append(cur_id)

  # return the NCBI accession IDs of the reference genomes
  # in the database directory (which is by default "genomes-of-common-contaminants")
  return ncbi_accession_ids


def get_taxonomy_ids() -> List[str]:
  """
  Given a list of NCBI accession IDs (from the get_fasta_ncbi_accession_ids() method)
  return the corresponding taxonomy IDs of those genomes.
  """
  # get the taxonomy id for this FASTA file
  # This will probably stay unimplemented because right now
  # it is difficult to download the 37 GB file needed to do this, at
  # https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
  pass


def get_taxonomy_ids_from_file(custom_taxonomy_ids_filename : str) -> List[str]:
  """
  Load a custom list of taxonomy IDs from a file.

  This removes the need to search the large file
  nucl_wgs.accession2taxid.gz (~37 GB uncompressed)
  for just the 20 genomes that we want to classify contaminants in.
  """
  taxonomy_ids_of_reference_genomes = []
  # using the custom seqid2taxid instead of searching the whole 37 GB mapping at
  # https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
  with open(custom_taxonomy_ids_filename, 'r') as fp:
    # getting the second column of each line,
    # which is the taxonomy ID of that accession ID
    for line in fp.readlines():
      tokens = line.strip().split()
      tax_id = tokens[1]
      taxonomy_ids_of_reference_genomes.append(tax_id)
      # try:
      #   int_tax_id = int(tax_id)
      #   # append int taxonomy id to the list
      #   taxonomy_ids_of_reference_genomes.append(tax_id)
      # except ValueError:
      #   print("Couldn't cast tax_id to integer value")
  return taxonomy_ids_of_reference_genomes


def main():
  # print(get_fasta_ncbi_accession_ids("genomes-of-common-contaminants"))
  # filenames have underscores (snake), directories have dashes (kebab)
  build_parent_map("../taxonomy", "../taxonomy/custom_taxonomy_ids.txt")


if __name__ == "__main__":
  main()