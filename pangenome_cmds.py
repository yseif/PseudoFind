

import pandas as pd
import os
from collections import Counter, defaultdict, OrderedDict
from glob import glob
import urllib
import re
import pickle
from collections import Counter, defaultdict
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
from glob import glob

def get_bbh(BR1out, BR2out, PID_threshold, output_file_name = '' ):
    
    BR2 = pd.read_table(BR2out, names= ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore'], index_col =['gene'])
    BR2 = BR2.loc[(BR2['PID'] > PID_threshold) & (BR2['eVal'] < 0.000001)] 

    BR1 = pd.read_table(BR1out , names= ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore'], index_col =['subject'])
    BR1 = BR1.loc[(BR1['PID'] > PID_threshold) & (BR1['eVal'] < 0.000001)]

    bbh = {}
    for gene in set(BR1.index) & set(BR2.index):
        hits_1 = BR1.loc[[gene]]
        subject_1 = hits_1.loc[hits_1['PID'] == max(hits_1['PID'])].gene.values[0]

        hits_2 = BR2.loc[BR2['subject'] == subject_1]
        if len(hits_2) == 0:
            continue
        gene_2 = hits_2.loc[hits_2['PID'] == max(hits_2['PID'])].index.values[0]

        if gene == gene_2:
            bbh[gene] = subject_1

    if output_file_name != '':
        pickle.dump(bbh, open(output_file_name, 'wb'))
    return bbh


def get_prokka2PATRIC_fna_mapping(prokka_fna, patric_fna):
    
    with open(prokka_fna, 'r') as f:
        string = f.read()
    gene_seqs_prokka = {x.split('\n')[0].split('Prokka|')[1].replace('>',''):''.join(x.split('\n')[1:]).lower() for x in string.split('\n>')}

    with open(patric_fna, 'r') as f:
        string = f.read()
    gene_seqs_patric = {x.split(' ')[0].replace('>','').split('\n')[0]:''.join(x.split('\n')[1:]).lower() for x in string.split('\n>')}

    for contig, seq in gene_seqs_patric.items():
        for letter in set(seq)- {'a','c','t','g'}:
            gene_seqs_patric[contig] = gene_seqs_patric[contig].replace(letter,'n')

    gene_seqs_patric_rev = {y:x for x,y in gene_seqs_patric.items()}

    patric2prokka = {gene_seqs_patric_rev[seq]:prokka_contig for prokka_contig, seq in gene_seqs_prokka.items() if seq in gene_seqs_patric_rev.keys()}

    return patric2prokka

def get_FT_from_gb(gb_file, output_file = ''):
    
    rows = []
    for seqrecord in SeqIO.parse(gb_file, 'genbank'):
        for feature in seqrecord.features:
            if feature.type == 'CDS':
                d = {x:y[0] for x,y in feature.qualifiers.items()}
                d.update({'start':int(feature.location.start), 'end':int(feature.location.end), 'strand':feature.location.strand, 'contig':seqrecord.id})
                rows.append(d)
    FT = pd.DataFrame(rows)
    
    if output_file != '':
        FT.to_csv(output_file, index = False)
    else:
        return FT

def get_fna_from_PATRIC(list_of_strains, output_directory):
    for gid in list_of_strains:
        F_genome_url = 'ftp://ftp.patricbrc.org/genomes/%s/%s.fna'%(gid,gid)
        F_page = urllib2.urlopen(F_genome_url, 'lxml')
        string = F_page.read()
        with open('%s/%s.fna'%(output_directory,gid), 'w') as y:
            y.write(string)

def annotate_fna_with_prokka(list_of_fna_files, output_directory, **kwargs):     
    for fa in list_of_fna_files:
        gid = fa.split('/')[-1].split('.fna')[0].replace('PATRIC_','')
        prokka_cmd = '/home/yara/prokka/bin/prokka --kingdom Bacteria --outdir %s/%s --locustag %s --prefix %s --evalue 0.001 --compliant %s'%(output_directory, gid, gid, gid, fa)
        print(gid, os.system(prokka_cmd))       
        
        
# annotation_platform = prokka
# file_extension = faa
# list_of_fasta_files is obtained from the prokka output

def reconfigure_fasta_for_blast(list_of_fasta_files, annotation_platform, file_extension):

    for fasta_file in list_of_fasta_files:

        with open(fasta_file, 'r') as f:
            string = f.read()
        
        if annotation_platform == 'patric':
            gene_seqs = {'fig|'+x.split('|')[1].split(' ')[0]:''.join(x.split('\n')[1:]) for x in string.split('\n>')}
        elif annotation_platform == 'prokka':
            gene_seqs = {x.split(' ')[0].replace('>',''):''.join(x.split('\n')[1:]) for x in string.split('\n>')}
        
        
        with open(fasta_file.split('.%s'%file_extension)[0]+'_m.%s'%file_extension, 'w') as y:
            for gene, seq in gene_seqs.items():
                if seq != '':
                    y.write('>'+gene+'\n'+seq+'\n')    

                    
                    
# build pan genome database from prokka output

def generate_pangenome_database(list_of_fasta_files, output_directory):
    with open(output_directory, 'w') as w:
        for faa_directory in faa_directories:
            with open(faa_directory, 'r') as f:
                string = f.read()
            gene_seqs = {x.split(' ')[0].replace('>',''):''.join(x.split('\n')[1:]) for x in string.split('\n>')}
            for gene, seq in gene_seqs.items():
                if seq != '':
                    w.write('>'+gene+'\n'+seq+'\n') 
    print('Done')
    
    
def run_bidirectional_blast(reference_fasta, list_of_fasta_files, blast_directory):
    
    if len(glob(blast_directory)) == 0:
        os.makedirs(blast_directory)
    if len(glob('%s/BR1'%blast_directory)) == 0:
        os.makedirs('%s/BR1'%blast_directory)
    if len(glob('%s/BR2'%blast_directory)) == 0:
        os.makedirs('%s/BR2'%blast_directory)    
        
    makeblastdb = 'makeblastdb -in %s -parse_seqids -dbtype prot'%reference_fasta
    o = os.system(makeblastdb)
    
    if o != 0:
        print('Warning: the fasta file for %s was not correctly configured. Blast database could not be constructed'%reference_fasta)
        return
    
    errors = defaultdict(dict)    
    for fasta in list_of_fasta_files:
        
        gid = fasta.split('/')[-1].split('.f')[0]
        print('Blasting: %s...'%gid)
        makeblastdb = 'makeblastdb -in %s -parse_seqids -dbtype prot'%fasta
        o = os.system(makeblastdb)
        if o != 0:
            errors[gid]['BlastDB'] = 1
            print('Warning: the fasta file for %s was not correctly configured. Blast database could not be constructed'%gid)
        
        BR1_out = '%s/BR1/%s'%(blast_directory, gid)
        blastp = 'blastp -query %s -db %s -outfmt 6 -out %s'%(reference_fasta, fasta, BR1_out)
        o = os.system(blastp)
        if o != 0:
            errors[gid]['Blastp1'] = 1
            print('Error: the fasta file for %s was not correctly configured or did not have a blast database'%gid)
                
        BR2_out =  '%s/BR2/%s'%(blast_directory, gid)
        blastp = 'blastp -query %s -db %s -outfmt 6 -out %s'%(fasta, reference_fasta, BR2_out)
        o = os.system(blastp)
        if o != 0:
            errors[gid]['Blastp2'] = 1
            print('Error: the fasta file for %s was not correctly configured or did not have a blast database'%gid)
    return errors     


def reconfigure_fasta_for_blast(list_of_fasta_files, annotation_platform, file_extension):

    for fasta_file in list_of_fasta_files:

        with open(fasta_file, 'r') as f:
            string = f.read()
        
        if annotation_platform == 'patric':
            gene_seqs = {'fig|'+x.split('|')[1].split(' ')[0]:''.join(x.split('\n')[1:]) for x in string.split('\n>')}
        elif annotation_platform == 'prokka':
            gene_seqs = {x.split(' ')[0].replace('>',''):''.join(x.split('\n')[1:]) for x in string.split('\n>')}
        
        
        with open(fasta_file.split('.%s'%file_extension)[0]+'_m.%s'%file_extension, 'w') as y:
            for gene, seq in gene_seqs.items():
                if seq != '':
                    y.write('>'+gene+'\n'+seq+'\n')    
 

def read_blast_results(blast_directory, PID_threshold, output_file_name):
    
    rows = {}
    for BR1_out in glob('%s/BR1/*'%blast_directory):
        
        BR2_out = BR1_out.replace('BR1','BR2')
        gid = BR1_out.split('/')[-1]
        BR2 = pd.read_table(BR2_out, names= ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore'], index_col =['gene'])
        BR2 = BR2.loc[BR2['PID'] > PID_threshold] # you can change your threshold here
        BR2 = BR2.loc[BR2['eVal'] < 0.000001]

        BR1= pd.read_table(BR1_out , names= ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore'], index_col =['subject'])
        BR1 = BR1.loc[BR1['PID'] > PID_threshold]
        BR1 = BR1.loc[BR1['eVal'] < 0.000001]


        bbh = {}
        for gene in set(BR1.index):
            hits_1 = BR1.loc[[gene]]
            subject_1 = hits_1.loc[hits_1['bitScore'] == max(hits_1['bitScore'])].gene.values[0]

            if gene not in list(BR2.index):
                continue

            hits_2 = BR2.loc[[gene]]
            subject_2 = hits_2.loc[hits_2['bitScore'] == max(hits_2['bitScore'])].subject.values[0]

            if subject_1 == subject_2:
                bbh[gene] = subject_1

            elif len(hits_1) > 1 and len(hits_2) > 1 and len(hits_1) == len(hits_2):
                double_hits = [x for x,y in Counter(list(hits_1['gene'])+list(hits_2['subject'])).items() if y > 1]
                if len(double_hits) > 1:
                    bbh[gene] = double_hits[0]
        rows[gid] = bbh    
    pickle.dump(rows, open('%s/%s.p'%(blast_directory,output_file_name),'w'))
    
    
def run_cdhit(fasta_files, pangenome_directory, threshold):
    
    with open(pangenome_directory, 'w') as y:
        for fa in fasta_files:
            with open(fa, 'r') as f:
                y.write(f.read())

    INPUT = pangenome_directory.split('.f')[0]            
    cdhit = 'cdhit -i %s.fa -o %s_cdhit -c %f -M 9000 -n 5 -d 0'%(INPUT, INPUT, threshold)
    print(os.system(cdhit))
    
                    
def get_cdhit_clusters(cdhit_directory):

    cdhit_clusters = defaultdict(set)

    with open(cdhit_directory, 'r') as f:
        for line in f.readlines():
            if 'Cluster' in line:
                cluster = line.split('>')[1].split('\n')[0]
            else:
                gene_id = line.split(' >')[-1].split('...')[0]
                cdhit_clusters[cluster].update({gene_id})
    return cdhit_clusters

def get_gene_family_matrix(cdhit_directory):
    cdhit_clusters = get_cdhit_clusters('%s.clstr'%cdhit_directory)
    cdhit_clusters_rev = {y:x for x,z in cdhit_clusters.items() for y in z}

    gid2clusters = defaultdict(dict)
    for cluster_id, genes in cdhit_clusters.items():
        for gene in genes:
            gid2clusters['_'.join(gene.split('_')[:-1])].update({cluster_id:1})

    gene_family_matrix = pd.DataFrame(gid2clusters).fillna(0)
    gene_family_matrix = gene_family_matrix.rename(columns = {x:x.replace('-','_') for x in gene_family_matrix.columns})
    
    gid2clusters2 = defaultdict(list)

    for cluster_id, genes in cdhit_clusters.items():
        for gene in genes:
            gid = '_'.join(gene.split('|')[-1].split('.pe')[0].split('_')[:-1])
            gid2clusters2[gid].append(cluster_id)

    duplicated = {GF for gid, clusters in gid2clusters2.items() for GF, count in Counter(clusters).items() if count > 1}
    
    return gene_family_matrix, duplicated

from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from collections import Counter, defaultdict, OrderedDict
import re


def get_all_mutations(ref_seq, sequences):
    mutations = defaultdict(list)
    for seq_2 in sequences.keys():
        found = 0
        count = 0
        aln = pairwise2.align.globalxx(sequences[ref_seq],sequences[seq_2])
        for char in aln[0][0]:
            if aln[0][1][count] != char:
                pos = len(aln[0][0][:count].replace('-',''))
                mutations[seq_2].append('%s%d%s'%(char,pos, aln[0][1][count]))
            count +=1
    return mutations

def get_insertions(muts_right):
    insertions_right = {}
    excluded = []
    for key, mut in muts_right.items():
        if key not in excluded:
            for ind in range(key, max(muts_right.keys())+1):
                if len(set(muts_right.keys()) & set(range(key, ind+1))) == len(set(range(key, max(muts_right.keys())+1))):
                    insertions_right.update({key:''.join([muts_right[i] for i in range(key,ind+1)])})
                    excluded += range(key,ind+1)
                    break
                if len(set(muts_right.keys()) & set(range(key, ind+1))) != len(set(range(key, ind+1))):
                    insertions_right.update({key:''.join([muts_right[i] for i in range(key,ind)])})
                    excluded += range(key,ind)
                    break
    return insertions_right

def get_all_insertions(seq_id, mutations):
    muts_left = OrderedDict(sorted({int(re.findall('[0-9]+', mut)[0]):mut[0] for mut in mutations[seq_id] if mut.endswith('-')}.items(), key = lambda a: a[0]) )
    muts_right = OrderedDict(sorted({int(re.findall('[0-9]+', mut)[0]):mut[-1] for mut in mutations[seq_id] if mut.startswith('-')}.items(), key = lambda a: a[0]) )
    insertions_left = get_insertions(muts_left)
    insertions_right = get_insertions(muts_right)
    r = defaultdict(str)
    l = defaultdict(str)
    for mut in mutations[seq_id]:
        if mut.startswith('-'):
            r[int(re.findall('[0-9]+', mut)[0])] += re.findall('[A-Z]+', mut)[0]
        elif mut.endswith('-'):
            l[int(re.findall('[0-9]+', mut)[0])] += re.findall('[A-Z]+', mut)[0]
            
    insertions_right.update({x:y for x,y in r.items() if len(y) > 1}) 
    insertions_left.update({x:y for x,y in l.items() if len(y) > 1}) 
    
    return {'deletions': insertions_left, 'insertions': insertions_right}

def get_point_mutations(all_insertions):
    point_mutations = []
    insertions = all_insertions['insertions'].copy()
    deletions = all_insertions['deletions'].copy()
    for loc, mut in all_insertions['deletions'].items():
        if len(mut) < 4:
            if loc+len(mut) in all_insertions['insertions'].keys():
                point_mutations.append('%s%d%s'%(mut,loc+1,all_insertions['insertions'][loc+len(mut)]))
                del deletions[loc]
                del insertions[loc+len(mut)]
    return point_mutations, deletions, insertions