#!/usr/bin/python
import os, argparse, re
from Bio import Entrez, SeqIO, SeqRecord
import pandas as pd

# Variables
LOCI="COI" # Irrelevant for now - edit in get_entrez_sequences_iterator and SeqRecords_iterator_to_TaxaRecords

ALLOW_SUB_SPECIES = False
REMOVE_IDS = False

class TaxaRecord:
    def __init__(self, taxid, taxonomy):
        global FAMILIES
        self.id = taxid
        self.taxonomy = taxonomy
        self.family = "" # Not implemented
        self.seqs = []
    
    def __str__(self):
        global REMOVE_IDS
        out = ""
        # if self.family:
        #     out += f">{self.family}_{self.id}_{self.seqs[0][0]}\n{self.seqs[0][1]}"
        # else:

        if REMOVE_IDS:
            out += f">{self.id}\n{self.seqs[0][1]}"

            for seqpair in self.seqs[1:]:
                out += f"\n>{self.id}\n{seqpair[1]}"
        else:
            out += f">{self.id}_{self.seqs[0][0]}\n{self.seqs[0][1]}"

            for seqpair in self.seqs[1:]:
                out += f"\n>{self.id}_{seqpair[0]}\n{seqpair[1]}"

        
        return out

    def add_seq(self, seq):
        self.seqs.append(seq)

    def get_seq_ids_str(self):
        return "".join([self.seqs[0][0]] + [f", {seq[0]}" for seq in self.seqs[1:]])

    def get_taxonomy_str(self):
        # print(self.taxonomy)
        return "".join([self.taxonomy[0]] + [f", {taxon}" for taxon in self.taxonomy[1:]])

    def get_genus(self):
        return self.id.split("_")[0]
    
    def get_species(self):
        return self.id.split("_")[1]

    def calc_family(self):
        print(self.id.replace('_', ' '))
        newFamily = None
        handle = Entrez.esearch(db="taxonomy", term=self.id.replace('_', ' '))
        record = Entrez.read(handle)
        if record["Count"] == "0":
            print(f"No family was found for {self.id}.")
            return None
        id_list = record["IdList"]
        handle = Entrez.efetch(db="taxonomy", id=id_list[0], retmode="xml")
        record = Entrez.read(handle)
        for entry in record[0]["LineageEx"]:
            if entry["Rank"] == "family":
                newFamily = entry["ScientificName"]
        if newFamily:
            self.family = newFamily
            return newFamily
        else:
            print(f"No family was found for {self.id}.")
            return None
            
# Get SeqRecords from NCBI Genbank
def get_entrez_sequences_iterator(term, retstart):

    # Search for all the COI sequences in the family
    handle = Entrez.esearch(db="nucleotide", term=f"{term}[Organism] AND CytB[Gene]", retstart=retstart*500, retmax=500)

    # Retrieve the list of matching sequence IDs
    record = Entrez.read(handle)
    id_list = record["IdList"]

    # Download the sequences in GB format
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb")
    records = SeqIO.parse(handle, "genbank")

    return records

# Parse SeqRecord list as custom TaxaRecord objs
def SeqRecords_iterator_to_TaxaRecords(records, taxrecs):

    global ALLOW_SUB_SPECIES
  
    counter = 0

    for record in records:

        seq = None

        # Extract COI feature(s)
        for feature in record.features:
            # Find features with CDS type and COI as gene name or product
            try:
                if feature.type == "CDS" and (feature.qualifiers["gene"][0].upper() in ["CYTB", "CYTOCHROME B"] and feature.qualifiers["product"][0].lower() in ["cytochrome b", "cytochrome-b"]):

                    # Extract the COI feature and create a new SeqRecord object
                    seq = str(feature.extract(record.seq))

                    if len(seq) < 20:
                        seq = None
                        

                    break

                if "organism" in feature.qualifiers:
                    organism = feature.qualifiers["organism"][0]

            except:
                continue
        
        # Extract taxon info if new entry

        if not seq:
            print("Record does not contain a COI feature. ID=" + record.id)
            continue

        if not organism:
            print("Record does not show organism source. ID=" + record.id)
            continue

        if len(organism.split()) > 2: 
            if ALLOW_SUB_SPECIES:
                organism = organism.split()[0] + " " + organism.split()[1] + "_" + organism.split()[2]
            else:
                print("Record has unexpected name. Taking first 2 words. " + organism + " ID=" + record.id)
                organism = organism.split()[0] + " " + organism.split()[1]

        taxid = organism.replace(" ", "_")
    
        # print(f"{taxid} : {record.id}")

        if taxid in taxrecs:
            taxrecs[taxid].add_seq((record.id, seq))
            counter += 1
        else:
            # Parse taxonomic information
            try:
                taxonomy = record.annotations['taxonomy']
                taxrecs[taxid] = TaxaRecord(taxid, taxonomy)
                taxrecs[taxid].add_seq((record.id, seq))
                counter += 1
            except:
                print("Record does not show taxonomic source. ID=" + record.id)
                taxrecs[taxid] = TaxaRecord(taxid, ["None"])
                taxrecs[taxid].add_seq((record.id, seq))
                counter += 1
        if counter % 500 == 0:
            print(counter)

    return (taxrecs, counter)

# Collect for search term into taxrecs obj
def collect_and_parse_records(term):
    addition = 0
    retstart = 0
    taxrecs = {}
    # Collect and parse records
    while True:
        addition = 0
        recordsIterator = get_entrez_sequences_iterator(term, retstart)
        taxrecs, addition = SeqRecords_iterator_to_TaxaRecords(recordsIterator, taxrecs)
        print(f"Parsed {str(addition)} records from GenBank.")
        retstart += 1

        print(f"RUN: {retstart}. Parsed {addition} new records")

        if not addition: break
    
    return taxrecs

# Create a datafield from taxa data
def get_df_for_taxa_data(taxrecs, term, filepath=None):

    # if not os.path.exists(filepath):
    #     return False

    filepath = os.path.join(filepath, f"{term}_taxa_data.csv")

    d={'family' : [], 'genus' : [], 'species' : [], 'taxonomy' : [], 'ids' : []}
    
    for taxid, rec in taxrecs.items():
        d['family'].append(rec.family)
        d['genus'].append(rec.get_genus())
        d['species'].append(rec.get_species())
        d['taxonomy'].append(rec.get_taxonomy_str())
        d['ids'].append(rec.get_seq_ids_str())

    df = pd.DataFrame(data=d)

    if filepath:
        df.to_csv(filepath)
    
    return df

# Write taxa to a FASTA file
def write_taxarecs_to_fasta(taxrecs, filepath=None):

    if filepath == None:
        filepath = "out.fas"

    with open(filepath, 'w') as f:
        for taxid, rec in taxrecs.items():
            f.write(str(rec) + "\n")
    
    return True

# Solve families WIP
def attempt_to_solve_families(taxrecs):

    for taxid, rec in taxrecs.items():
        rec.calc_family()
    
    return taxrecs

# Read taxa records from a custom FASTA file created with this program
def read_taxrecs_from_file(filepath):

    taxrecs = {}
    counter = 0

    with open(filepath, 'r') as f:
        
        # Define collection variables
        currentTaxID = ""
        currentSeqID = ""
        currentSeq = ""

        # Read lines in file
        for line in f:
            # If header
            if line.startswith(">"):

                

                #If new header add last entry to list
                if currentTaxID:
                    counter += 1  
                    if currentTaxID in taxrecs:
                        taxrecs[currentTaxID].add_seq((currentSeqID, currentSeq))
                    else:
                        # print("New Taxa")
                        taxrecs[currentTaxID] = TaxaRecord(currentTaxID, ["Not Specified"])
                        taxrecs[currentTaxID].add_seq((currentSeqID, currentSeq))
                
                currentSeq = ""                             
                
                #Get header
                header = line.strip().lstrip(">")
                fields = header.split("_")
                
                #Collect details
                if len(fields) == 4:
                    currentTaxID = fields[1] + '_' + fields[2]
                    currentSeqID = fields[3]
                    # print(currentTaxID + " " + str(counter))
                else:
                    #skip this taxa
                    currentTaxID=""

            else:
                #Read sequence
                currentSeq += line.strip()
        
        # add the last entry
        if currentTaxID:
            counter += 1  
            if currentTaxID in taxrecs:
                taxrecs[currentTaxID].add_seq((currentSeqID, currentSeq))
            else:
                # print("New Taxa")
                taxrecs[currentTaxID] = TaxaRecord(currentTaxID, ["Not Specified"])
                taxrecs[currentTaxID].add_seq((currentSeqID, currentSeq))

    print(f"Successfully parsed {str(counter)} records from {filepath}.")

    return taxrecs

# Prune by largest sequence
def prune_by_largest_coverage(taxrecs, length):

    rm = []
    rmNno = 0

    maxlen = length[1]
    optlen = length[0]

    for taxid, rec, in taxrecs.items():
    
        pruned = [seq for seq in rec.seqs if len(seq[1]) < maxlen and len(seq[1]) > 0.15*maxlen]

        if len(pruned) == 0:
            rm.append(taxid)
            # del taxrecs[taxid]
            continue

        remN = [seq for seq in pruned if not any(letter in seq[1] for letter in ["N", "X"])]

        if len(remN) == 0:
            print("No choice but to use seq with missing data")
            print(pruned)
            rmNno += 1
        else:
            pruned = remN


        closest = None
        closest_distance = float('inf')
        for tup in rec.seqs:
            seq_len = len(tup[1])
            distance = abs(seq_len - optlen)
            if distance < closest_distance:
                closest = tup
                closest_distance = distance

        rec.seqs = [closest]

    for rmkey in rm:
        del taxrecs[rmkey]

    print("Number of removed taxa: " + str(len(rm)))
    print("Removed Taxa: ")
    print(rm)
    print(f"Number of sequences kept with missing data = {rmNno}")

    return taxrecs
        
def validate_email(email):
    # check if the email address is valid
    if not re.match(r"[^@]+@[^@]+\.[^@]+", email):
        raise argparse.ArgumentTypeError(f"{email} is not a valid email address")
    return email   

# Define the argument parser
parser = argparse.ArgumentParser(description='Process NCBI GenBank Records into custom format')

# Mutually exclusive arguments
group = parser.add_mutually_exclusive_group()
group.add_argument('-f', '--input', dest='input_file', metavar='PATH',
                    help='Input file of taxrecords')
group.add_argument('-df', '--Taxonomy_info_csv', dest='taxonomy_file', metavar='PATH', nargs='?', const='.',
                    help='Produce a CSV file of taxonomic information')
# Required argument
parser.add_argument('-e', '--email', dest='email',
                    type=validate_email, help='Input email address')

# Optional argument
parser.add_argument('-o', '--output', dest='output_file', metavar='PATH',
                    help='Path for output file')
parser.add_argument('--allow-sub-species', dest='ALLOW_SUB_SPECIES', action='store_true',
                    help='Allow sub-species in taxonomic records')
parser.add_argument('--no-allow-sub-species', dest='ALLOW_SUB_SPECIES', action='store_false',
                    help='Do not allow sub-species in taxonomic records (default)')
parser.add_argument('--remove_ids', dest='REMOVE_IDS', action='store_true',
                    help='Remove accession numbers from FASTA file (default)')
parser.add_argument('--keep_ids', dest='REMOVE_IDS', action='store_false',
                    help='Keep accession numbers in FASTA file')
parser.add_argument('-rmdup', '--remove_duplicates', nargs='*', type=int, 
                    metavar='Optimal_Length Max_Length', action="store",
                    help='Remove duplicate taxa records by optimal length')


parser.add_argument('-s', '--search_term', dest='search_term', type=str, nargs='?',
                    help='Clade/Group term for NCBI GenBank Search (Required if -f not used)')
parser.set_defaults(allow_sub_species=False) 
parser.set_defaults(allow_sub_species=False) 

args = parser.parse_args()  

ALLOW_SUB_SPECIES = args.ALLOW_SUB_SPECIES
REMOVE_IDS = args.REMOVE_IDS

# Set email
Entrez.email = args.email

# Define the tax dict
taxrecs = {}

# call the appropriate function based on the command line arguments
if args.input_file:
    if os.path.exists(args.input_file):
        print("Reading taxa records from file: " + str(args.input_file))
        taxrecs = read_taxrecs_from_file(str(args.input_file))
else:
    if not args.email:
        parser.error('Email is required to download from NCBI GenBank.')
    if not args.search_term:
        parser.error('Search term is required if -f is not provided.')
    if len(str(args.search_term).split()) > 1:
        parser.error('Please provide only one search term')
    print("Collecting records from NCBI Genbank with search term: " + str(args.search_term))
    taxrecs = collect_and_parse_records(args.search_term)

if args.remove_duplicates is not None:
    if len(args.remove_duplicates) not in (0, 2):
        parser.error('Specify two numbers or none for defaults (1100, 1300)')
    if len(args.remove_duplicates) == 0:
        print("Removing duplicate taxa records based on length: opt=1100, max=1300")
        taxrecs = prune_by_largest_coverage(taxrecs, [1100, 1300])
    else:
        print(f"Removing duplicate taxa records based on length: opt={args.remove_duplicates[0]} max={args.remove_duplicates[1]}")
        taxrecs = prune_by_largest_coverage(taxrecs, args.remove_duplicates)

if args.output_file:
    write_taxarecs_to_fasta(taxrecs, str(args.output_file))
else:
    if args.input_file:
        write_taxarecs_to_fasta(taxrecs, str(args.input_file) + "_o.fas")
    else:
        write_taxarecs_to_fasta(taxrecs, "seqprog_out.fas")

if args.taxonomy_file:
    df = get_df_for_taxa_data(taxrecs, args.search_term, args.taxonomy_file)

    


    