# SEQPROG

Requires Python 3.7+ with following packages:

* Argparse
* BioPython
* Pandas

Reads and writes sequences in common format

## Usage

>   Usage: seqprog.py [-h] [-f PATH | -df PATH] -e EMAIL [-o PATH] [-rmdup]

    Process NCBI GenBank Records into custom format

    optional arguments:
    -h, --help            show this help message and exit
    -f PATH, --input PATH
                            Input file of taxrecords
    -df PATH, --Taxonomy_info_csv PATH
                            Produce a CSV file of taxonomic information
    -e EMAIL, --email EMAIL
                            Input email address
    -o PATH, --output PATH
                            Path for output file
    -rmdup, --remove_duplicates
                            Remove duplicate taxa records by optimal length
>