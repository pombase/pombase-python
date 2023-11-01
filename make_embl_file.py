#!/usr/bin/python3

# Read a contig file and write a EMBL format file for submission to ENA
#
# # Run as:
#   cd pombe-embl
#   for i in 1 2 3
#   do
#      PATH_TO_POMBASE_PYTHON/make_embl_file.py \
#          ./ftp_site/pombe/names_and_identifiers/PomBase2UniProt.tsv \
#          embl_resubmission/CU*_chromosome$i.embl chromosome$i.contig \
#           embl_resubmission/chromosome$i.embl
#   done
#
# # then in pombe-embl/embl_resubmission run:
#    ./chr_submit.sh



import sys
import io
import csv
import re
import warnings

warnings.filterwarnings('ignore', '.*Non-standard molecule type: genomic DNA.*',)

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

def process_dbxrefs(feature, systematic_id):
    db_xrefs = feature.qualifiers.pop('db_xref', [])
    db_xrefs = [x for x in db_xrefs if x.startswith('PMID:')]

    if 'locus_tag' in feature.qualifiers:
        db_xrefs.append('PomBase:' + systematic_id)
        if systematic_id in pombe_uniprot_map and feature.type not in ["5'UTR", "3'UTR"]:
            uniprot_id = pombe_uniprot_map[systematic_id]
            db_xrefs.append('UniProtKB/Swiss-Prot:' + uniprot_id)

    feature.qualifiers['db_xref'] = db_xrefs

def process_product(qualifiers):
    if 'product' in qualifiers:
        new_product = re.sub(' \(predicted\)$', '', qualifiers['product'][0])
        qualifiers['product'][0] = new_product

def is_utr(feature):
    return feature.type in ["5'UTR", "3'UTR"]

def process_qualifers(feature):
    remove_non_embl_qualfiers(feature)
    sys_id = None

    qualifiers = feature.qualifiers
    if 'systematic_id' in qualifiers:
        sys_ids = qualifiers['systematic_id']
        if len(sys_ids) > 1:
            die('too many /systematic_id qualifiers' + str(feature))
        del qualifiers['systematic_id']
        sys_id = sys_ids[0]

        original_sys_id = sys_id

        if feature.type in ['CDS', 'mRNA', "5'UTR", "3'UTR"]:
            sys_id = re.sub(r"\.\d$", '', sys_id)

        is_first_transcript = original_sys_id == sys_id or original_sys_id.endswith('.1')

        qualifiers['locus_tag'] = ['SPOM_' + sys_id]
        if is_first_transcript and sys_id in protein_id_map and not is_utr(feature):
            qualifiers['protein_id'] = protein_id_map[sys_id]
    else:
        if feature.type not in ['gap']:
            eprint('no /systematic_id in:')
            sys.exit(feature)

    primary_name = qualifiers.pop('primary_name', None)
    if primary_name != None and not is_utr(feature):
        qualifiers['gene'] = primary_name[0]

    synonym = qualifiers.pop('synonym', None)
    if synonym != None and not is_utr(feature):
        qualifiers['gene_synonym'] = synonym

    # avoid /pseudo="" in the output:
    if 'pseudo' in qualifiers:
        qualifiers['pseudo'] = None

    if 'ribosomal_slippage' in qualifiers:
        qualifiers['ribosomal_slippage'] = None

    process_product(qualifiers)

    process_dbxrefs(feature, sys_id)

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

seen_feature_locs = {}

for feature in contig.features:
    if feature.type not in types_to_keep:
        continue

    if feature.type == 'LTR':
        process_ltr(feature)

    if feature.type in ['lncRNA', 'snoRNA', 'sncRNA', 'snRNA']:
        process_rna(feature)

    process_qualifers(feature)

    key = feature.type + '--' + str(feature.location)

    if key in seen_feature_locs and (feature.type == "5'UTR" or feature.type == "3'UTR"):
        print('warning: ígnoring duplicate: ' + key)
    else:
        if key in seen_feature_locs:
            print("warning: found duplicate feature that isn't a UTR: " +
                  f"{feature.type} {feature.location}")

        seen_feature_locs[key] = True
        new_features.append(feature)

contig.features = new_features

SeqIO.write(contig, output_file_name, 'embl')
