# computational-genomics-team-47
computational-genomics-team-47 Fall 2023 - Dr. Langmead's Computational Genomics: Sequences

## Specific Reproducible Runs of the Program

Below we describe specific runs of the program, reproducible on the ugrad machines, and their expected results.

- 


## How to run the program / Usage
- `unzip taxonomy/names-and-nodes.zip -d taxonomy/`
  - This unzips the `names.dmp` and `nodes.dmp` file into the `taxonomy` directory.
- `python src/main.py --input-query covid-assemblies/covid-assembly-1.txt`
  - The program can be run on any of the 9 input assemblies, e.g.
  - `python src/main.py --input-query covid-assemblies/covid-assembly-9.txt`

### Additional arguments
Here is an invocation of the program providing all optional arguments that are available:

- `python3 src/main.py --db new-tutorial-reference-database/ --input-query covid-assemblies/covid-assembly-1.txt --taxonomy taxonomy --taxonomy-ids taxonomy/custom_taxonomy_ids.txt --k 31`

The arguments are outlined below:
- `--db`
  - Name of the directory containing the database of known contaminants which you want to cross-check the input sequence for (default: `new-tutorial-reference-database`)
- `--input-query`
  - Filename of the sequence in which the program will search for contaminants (required)
- `--taxonomy`
  - Name of the directory containing the taxonomy (including `names.dmp` and `nodes.dmp`) (default: `taxonomy`)
- `--taxonomy-ids`
  - Name of the plaintext file (default: `taxonomy/custom_taxonomy_ids.txt`) which contains taxonomy ids sourced from NCBI corresponding to the NCBI accession IDs of the FASTA files in the database
- `--k`
  - k, the length of the kmer (default: k = 12, which runs on the ugrad machines well (within memory constraints). k = 31 is ideal if the computer has enough memory.)


## Overview

Our program has 5 steps, which are identified in the below picture:

![Image summary of potential plan for how our program works](images/summary_of_planned_program.png)

Each step roughly corresponds to a separate `.py` file.

The `main.py` file is the one that runs the entire contamination classification program.

## Repository Folder Structure Summary
- `src` contains the Python program that classifies contaminants in reads user a exact matching k-mer database (and the least common ancestor algorithm and root-to-leaf paths).
- The program is structured chronologically like below:
  - (0.) `main.py` (main method starts here)
  - (1.) `taxonomy_tree.py` 
  - (2.) `kmer_to_lca_mapping.py`
  - (3.) `psuedoreads.py`
  - (4.) `get_kmer_hit_counts.py`
  - (5.) `print_summary_contaminants_found()` (in `main.py`)
- `genomes-of-common-contaminants` contains approximately 20 genomes (totalling ~90 MB) of bacteria and viruses that common contaminate DNA sequences (Mycoplasma, Eschericia lambda phage phiX174, etc.)
  - `gocc-shortened` is a condensed version of this database.
- `covid-assemblies` contains 9 input sequences that are fed into the program and checked for contamination.
- `taxonomy` contains:
  - `nodes.dmp` (~243 MB) which is a file of the taxonomy of all organisms in the tree of life, each represented with a unique integer taxonomy id.
  - `names.dmp` (~183 MB) which maps each integer taxonomy id to the common plaintext name of the species.

## Dependencies
### `python` 3
### `pip`
- None so far.

## Authors
- Dhruv Dubey
- Mitra Harpale
- Christopher Li
- Jaeyoon Wang
