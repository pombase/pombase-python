#!/usr/bin/python3

# Read a contig file and write a EMBL format file for submission to ENA
#
# Run as:
#   cd pombe-embl
#   for contig in *.contig
#   do
#      PATH_TO_POMBASE_PYTHON/make_embl_file.py ./ftp_site/pombe/names_and_identifiers/PomBase2UniProt.tsv contig > embl_resubmission/$contig
#   done

import sys
import io
import csv

from Bio import SeqIO

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def die(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.exit(1)

if len(sys.argv) <= 2:
    sys.exit('needs 2 arguments')

uniprot_mapping_file_name = sys.argv[1]

pombe_uniprot_map = {}

with open(uniprot_mapping_file_name, mode='r') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        pombe_uniprot_map[row[0]] = row[1]

file_name = sys.argv[2]

eprint('processing: ' + file_name)

types_to_keep = set(['CDS', "3'UTR", "5'UTR", 'LTR', 'gap', 'tRNA', 'rRNA',
                     'ncRNA', 'lncRNA', 'sncRNA', 'snoRNA', 'snRNA'])

def remove_non_embl_qualfiers(feature):
    for key in ['SO', 'colour', 'controlled_curation', 'note', 'EC_number',
                'gene', 'obsolete_name', 'GO', 'shared_id']:
        feature.qualifiers.pop(key, None)

def add_dbxrefs(feature):
    if 'locus_tag' in feature.qualifiers:
        systematic_id = feature.qualifiers['locus_tag'][0]
        db_xrefs = feature.qualifiers.pop('db_xref', [])
        db_xrefs.append('PomBase:' + systematic_id)
        if systematic_id in pombe_uniprot_map:
            uniprot_id = pombe_uniprot_map[systematic_id]
            db_xrefs.append('UniProtKB/Swiss-Prot:' + uniprot_id)
        feature.qualifiers['db_xrefs'] = db_xrefs

def process_qualifers(feature):
    remove_non_embl_qualfiers(feature)

    qualifiers = feature.qualifiers
    if 'systematic_id' in qualifiers:
        sys_ids = qualifiers['systematic_id']
        if len(sys_ids) > 1:
            die('too many /systematic_id qualifiers' + str(feature))
        del qualifiers['systematic_id']
        qualifiers['locus_tag'] = sys_ids
    else:
        if feature.type not in ['gap']:
            eprint('no /systematic_id in:')
            sys.exit(feature)

    primary_name = qualifiers.pop('primary_name', None)
    if primary_name != None:
        qualifiers['standard_name'] = primary_name

    synonym = qualifiers.pop('synonym', None)
    if synonym != None:
        qualifiers['gene_synonym'] = synonym

    add_dbxrefs(feature)

def process_ltr(feature):
    feature.type = 'repeat_region'
    feature.qualifiers['rpt_type'] = 'long_terminal_repeat'

def process_rna(feature):
    feature_type = feature.type
    feature.type = 'ncRNA'
    feature.qualifiers['ncRNA_class'] = feature_type

with open(file_name) as contig_in:
    contig = SeqIO.read(file_name, 'embl')

    new_features = []

    for feature in contig.features:
        if feature.type not in types_to_keep:
            continue

        process_qualifers(feature)

        if feature.type == 'LTR':
            process_ltr(feature)

        if feature.type in ['lncRNA', 'snoRNA', 'snRNA']:
            process_rna(feature)

        new_features.append(feature)

    contig.features = new_features

    SeqIO.write(contig, file_name + '.embl_fix', 'embl')
