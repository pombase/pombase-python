#!/usr/bin/python3

# Read a contig file and write a EMBL format file for submission to ENA
#
# Run as:
#   cd pombe-embl
#   for i in 1 2 3
#   do
#      PATH_TO_POMBASE_PYTHON/make_embl_file.py \
#          ./ftp_site/pombe/names_and_identifiers/PomBase2UniProt.tsv \
#          embl_resubmission/CU*_chromosome$i.embl chromosome$i.contig \
#           embl_resubmission/chromosome$i.embl
#   done

import sys
import io
import csv
import re

from Bio import SeqIO

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def die(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.exit(1)

if len(sys.argv) != 5:
    sys.exit('needs 4 arguments: uniprot_mapping entry_from_embl contig_file output_file')

uniprot_mapping_file_name = sys.argv[1]

pombe_uniprot_map = {}

with open(uniprot_mapping_file_name, mode='r') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        pombe_uniprot_map[row[0]] = row[1]

entry_from_embl_file_name = sys.argv[2]
input_file_name = sys.argv[3]
output_file_name = sys.argv[4]

eprint('processing: ' + input_file_name)

protein_id_map = {}

entry_from_embl = SeqIO.read(entry_from_embl_file_name, 'embl')

for feature in entry_from_embl.features:
    qualifiers = feature.qualifiers

    if 'locus_tag' in qualifiers:
        locus_tags = qualifiers['locus_tag']
        if 'protein_id' in qualifiers:
            protein_ids = qualifiers['protein_id']
            protein_id_map[locus_tags[0]] = protein_ids[0]

types_to_keep = set(['CDS', "3'UTR", "5'UTR", 'LTR', 'gap', 'tRNA', 'rRNA',
                     'ncRNA', 'lncRNA', 'sncRNA', 'snoRNA', 'snRNA'])

def remove_non_embl_qualfiers(feature):
    for key in ['SO', 'colour', 'controlled_curation', 'note', 'EC_number',
                'partial', 'feature_source',
                'gene', 'obsolete_name', 'GO', 'shared_id']:
        feature.qualifiers.pop(key, None)

def add_dbxrefs(feature):
    if 'locus_tag' in feature.qualifiers:
        systematic_id = feature.qualifiers['locus_tag'][0]
        db_xref = feature.qualifiers.pop('db_xref', [])
        db_xref.append('PomBase:' + systematic_id)
        if systematic_id in pombe_uniprot_map:
            uniprot_id = pombe_uniprot_map[systematic_id]
            db_xref.append('UniProtKB/Swiss-Prot:' + uniprot_id)
        feature.qualifiers['db_xref'] = db_xref

def process_product(qualifiers):
    if 'product' in qualifiers:
        new_product = re.sub(' \(predicted\)$', '', qualifiers['product'][0])
        qualifiers['product'][0] = new_product

def process_qualifers(feature):
    remove_non_embl_qualfiers(feature)

    qualifiers = feature.qualifiers
    if 'systematic_id' in qualifiers:
        sys_ids = qualifiers['systematic_id']
        if len(sys_ids) > 1:
            die('too many /systematic_id qualifiers' + str(feature))
        del qualifiers['systematic_id']
        qualifiers['locus_tag'] = sys_ids
        sys_id = sys_ids[0]
        if sys_id in protein_id_map:
            qualifiers['protein_id'] = protein_id_map[sys_id]
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

    process_product(qualifiers)

    add_dbxrefs(feature)

def process_ltr(feature):
    feature.type = 'repeat_region'
    feature.qualifiers['rpt_type'] = 'long_terminal_repeat'

def process_rna(feature):
    feature_type = feature.type

    if 'SO' in feature.qualifiers:
        so_term = feature.qualifiers['SO'][0]
        if so_term == 'SO:0002247':
            feature_type = 'other'

    feature.type = 'ncRNA'
    if feature_type != None:
        feature.qualifiers['ncRNA_class'] = feature_type

contig = SeqIO.read(input_file_name, 'embl')

new_features = []

for feature in contig.features:
    if feature.type not in types_to_keep:
        continue

    if feature.type == 'LTR':
        process_ltr(feature)

    if feature.type in ['lncRNA', 'snoRNA', 'sncRNA', 'snRNA']:
        process_rna(feature)

    process_qualifers(feature)

    new_features.append(feature)

contig.features = new_features

SeqIO.write(contig, output_file_name, 'embl')
