# computational-genomics-team-47
computational-genomics-team-47 Fall 2023 - Dr. Langmead's Computational Genomics: Sequences

## File Summary
- `src` contains the Python program that classifies contaminants in reads user a exact matching k-mer database (and the least common ancestor algorithm and root-to-leaf paths).
- `genomes-of-common-contaminants` contains approximately 20 genomes (totalling ~90 MB) of bacteria and viruses that common contaminate DNA sequences (Mycoplasma, Eschericia lambda phage phiX174, etc.)
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
