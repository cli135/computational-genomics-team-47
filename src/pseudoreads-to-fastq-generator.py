def read_fasta_file(fasta_file_path):
    # Reads a FASTA file and returns the sequence as a string.
    # Parameters:
    # fasta_file_path (str): The path to the FASTA file.
    # Returns:
    # str: The genome sequence from the FASTA file.
    sequence = ""
    with open(fasta_file_path, 'r') as file:
        for line in file:
            # Skip the header line
            if line.startswith('>'):
                continue
            # Add the sequence line to the cumulative sequence
            sequence += line.strip()
    return sequence

def split_genome_into_pseudo_reads_from_fasta(fasta_file_path, output_file, read_length=100, overlap=50):
    # Splits a genome sequence from a FASTA file into overlapping pseudo-reads.
    # Parameters:
    # fasta_file_path (str): The path to the FASTA file containing the genome sequence.
    # read_length (int): The length of each pseudo-read.
    # overlap (int): The length of the overlap between consecutive reads.
    # Returns:
    # list: A list of pseudo-reads.
    
    # Read the genome sequence from the FASTA file
    genome_sequence = read_fasta_file(fasta_file_path)
    # Calculate the step size for the next read (read length minus overlap)
    step_size = read_length - overlap
    # Initialize an empty list to store the pseudo-reads
    pseudo_reads = []
    # Iterate over the genome sequence to create pseudo-reads
    for i in range(0, len(genome_sequence) - read_length + 1, step_size):
        # Extract the pseudo-read from the genome sequence
        read = genome_sequence[i:i + read_length]
        # Add the pseudo-read to the list
        pseudo_reads.append(read)


    # Write the pseudo-reads to a FASTQ file
    # USED FOR TESTING COMPARISONS VIA KRAKEN

    with open(output_file, 'w') as output_file:
    # Iterate over the pseudo-reads to write them to the FASTQ file
        for i, read in enumerate(pseudo_reads):
            # Generate a unique identifier for each pseudo-read
            read_id = f"@Read{i + 1}"
            # Write the read information (ID, sequence, and quality scores) to the FASTQ file
            output_file.write(f"{read_id}\n{read}\n+\n{'I' * len(read)}\n")  # Assuming quality scores are all 'I'

    return pseudo_reads


# pseudo_reads = 
split_genome_into_pseudo_reads_from_fasta('covid-assemblies/covid-assembly-1.txt',"kraken-fastq-covid-assembly-1.fastq")
split_genome_into_pseudo_reads_from_fasta('covid-assemblies/covid-assembly-2.txt',"kraken-fastq-covid-assembly-2.fastq")
split_genome_into_pseudo_reads_from_fasta('covid-assemblies/covid-assembly-3.txt',"kraken-fastq-covid-assembly-3.fastq")
split_genome_into_pseudo_reads_from_fasta('covid-assemblies/covid-assembly-4.txt',"kraken-fastq-covid-assembly-4.fastq")
split_genome_into_pseudo_reads_from_fasta('covid-assemblies/covid-assembly-5.txt',"kraken-fastq-covid-assembly-5.fastq")
split_genome_into_pseudo_reads_from_fasta('covid-assemblies/covid-assembly-6.txt',"kraken-fastq-covid-assembly-6.fastq")
split_genome_into_pseudo_reads_from_fasta('covid-assemblies/covid-assembly-7.txt',"kraken-fastq-covid-assembly-7.fastq")
split_genome_into_pseudo_reads_from_fasta('covid-assemblies/covid-assembly-8.txt',"kraken-fastq-covid-assembly-8.fastq")
split_genome_into_pseudo_reads_from_fasta('covid-assemblies/covid-assembly-9.txt',"kraken-fastq-covid-assembly-9.fastq")
