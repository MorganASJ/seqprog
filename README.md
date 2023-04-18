# SEQPROG

Requires Python 3.6+ with following packages:

* Argparse
* BioPython
* Pandas

Reads and writes sequences in common format

## Using the program

>   Usage: python3 seqprog.py [-h] [-f PATH | -df [PATH]] [-e EMAIL] [-o PATH]
                    [--allow-sub-species] [--no-allow-sub-species]
                    [--remove_ids] [--keep_ids]
                    [-rmdup [Optimal_Length Max_Length [Optimal_Length Max_Length ...]]]
                    [-s [SEARCH_TERM]]

    Process NCBI GenBank Records into custom format

    optional arguments:
    -h, --help            show this help message and exit
    -f PATH, --input PATH
                            Input file of taxrecords
    -df [PATH], --Taxonomy_info_csv [PATH]
                            Produce a CSV file of taxonomic information
    -e EMAIL, --email EMAIL
                            Input email address
    -o PATH, --output PATH
                            Path for output file
    --allow-sub-species   Allow sub-species in taxonomic records
    --no-allow-sub-species
                            Do not allow sub-species in taxonomic records
                            (default)
    --remove_ids          Remove accession numbers from FASTA file (default)
    --keep_ids            Keep accession numbers in FASTA file
    -rmdup [Optimal_Length Max_Length [Optimal_Length Max_Length ...]], --remove_duplicates [Optimal_Length Max_Length [Optimal_Length Max_Length ...]]
                            Remove duplicate taxa records by optimal length
    -s [SEARCH_TERM], --search_term [SEARCH_TERM]
                            Clade/Group term for NCBI GenBank Search (Required if
                            -f not used)
>
