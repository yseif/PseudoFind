import pickle
import pandas as pd
import cobra
from collections import defaultdict, OrderedDict, Counter
from cobra import Reaction, Metabolite, Model, Gene
import time
from glob import glob

from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches
import sympy
from sympy import to_dnf, Add
from Bio import SeqIO
import re
import os
import scipy
# from auxotrophy_script import *
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import pairwise2


def run_blastx(fasta, output_dir, blastdb, sabir):
    BR1_out = '%s/%s'%(output_dir, sabir)
    blastn = 'blastx -query %s -db %s -max_target_seqs 5000 -outfmt 6 -out %s -reward 2 -penalty -3'%(fasta, blastdb ,BR1_out)
    print(os.system(blastn))
    
def run_blastn(fasta, output_dir, blastdb, sabir):
    BR1_out = '%s/%s'%(output_dir, sabir)
    blastn = 'blastn -query %s -db %s -max_target_seqs 5000 -outfmt 6 -out %s -word_size 11 -evalue 10 -reward 2 -penalty -3 -gapopen 5 -gapextend 2'%(fasta, blastdb ,BR1_out)
    print(os.system(blastn))
    

def check_seq(strand, gene_seqs, contig, BR_s, i, K1, K2):
    if strand == '+':
        seq2m = gene_seqs[contig][BR_s.iloc[i]['queryStart'] - 1 - K1: BR_s.iloc[i]['queryEnd'] +K2]
    else:
        seq2m = str(Seq(gene_seqs[contig][BR_s.iloc[i]['queryStart'] - 1 -K2: BR_s.iloc[i]['queryEnd'] + K1], IUPAC.unambiguous_dna).reverse_complement())     
    return seq2m

def get_CDS_hits(BR_out, ref_gene_seqs, fasta):

    BR = pd.read_table(BR_out, names= ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore'])
    BR['alnPerc'] = [abs(BR.loc[index, 'queryStart'] - BR.loc[index, 'queryEnd'])/len(ref_gene_seqs[BR.loc[index,'subject']]) for index in BR.index]
    BR = BR.set_index('gene')
    BR = BR.loc[(BR['PID'] > 80) & (BR['eVal'] < 10**(-8)) & (BR['alnPerc'] > 0.8)]

    with open(fasta, 'r') as f:
        string = f.read()
    gene_seqs = {x.split('\n')[0].replace('>','').split(' ')[0]:''.join(x.split('\n')[1:]).upper() for x in string.split('\n>')}

    cds_hits_2 = defaultdict(dict)
    for subject in set(BR['subject']):

        BR_s = BR.loc[BR['subject'] == subject]
        for i in range(len(BR_s)):
            contig = BR_s.index.tolist()[i]
            queryEnd = BR_s.iloc[i]['queryEnd']
            queryStart = BR_s.iloc[i]['queryStart']
            locus_details = '%s-s%d:e%d'%(contig,queryStart,queryEnd)

            if queryStart > 20 and len(gene_seqs[contig]) - queryEnd > 20:
                strand = '+' if BR_s.iloc[i]['subjectStart'] < BR_s.iloc[i]['subjectEnd'] else '-' 
                seq2 = gene_seqs[contig][BR_s.iloc[i]['queryStart'] - 1: BR_s.iloc[i]['queryEnd'] ]
                seq2 = seq2 if strand == '+' else str(Seq(seq2, IUPAC.unambiguous_dna).reverse_complement())

                K1, K2 = 0,0
                # check start codon:
                if str(Seq(seq2).translate()[:-1]).find('*') != -1:
                    upstream = [K for K in range(30) if str(Seq(check_seq(strand, gene_seqs, contig, BR_s, i, K,K2)).translate())[:-1].find('*') == -1]
                    downstream = [K for K in range(30) if str(Seq(check_seq(strand, gene_seqs, contig, BR_s, i, -K,K2)).translate())[:-1].find('*') == -1]

                    if len(upstream)> 0 and len(downstream) > 0:
                        K1 = min(upstream) if min(upstream) < min(downstream) else -min(downstream)
                    elif len(upstream) > 0:
                        K1 = min(upstream)
                    elif len(downstream) > 0:
                        K1 = -min(downstream)


                # check end codon:
                if str(Seq(seq2).translate())[-1] != '*':
                    upstream = [K for K in range(30) if str(Seq(check_seq(strand, gene_seqs, contig, BR_s, i, K1,K)).translate())[-1] == '*']
                    downstream = [K for K in range(30) if str(Seq(check_seq(strand, gene_seqs, contig, BR_s, i, K2,-K)).translate())[-1] == '*']

                    if len(upstream)> 0 and len(downstream) > 0:
                        K2 = min(upstream) if min(upstream) < min(downstream) else -min(downstream)
                    elif len(upstream) > 0:
                        K2 = min(upstream)
                    elif len(downstream) > 0:
                        K2 = -min(downstream)

                seq2 = check_seq(strand, gene_seqs, contig, BR_s, i, K1, K2)  

                cds_hits_2[subject].update({locus_details+'(%s)'%str(strand):seq2.upper()})

    return cds_hits_2

def get_all_mutations(ref_seq, seq_id2, sequences):
    mutations = []
    count = 0
    if len(sequences[ref_seq]) <= 1500:
        aln = pairwise2.align.globalxx(sequences[ref_seq],sequences[seq_id2])
    else:
        with open('/home/yara/Downloads/temp.ffn', 'w')as y:
            y.write('>'+ref_seq+'\n'+str(sequences[ref_seq])+'\n')
            y.write('>'+seq_id2+'\n'+str(sequences[seq_id2])+'\n')

        muscle_in = '/home/yara/Downloads/temp.ffn'
        muscle_out = '/home/yara/Downloads/temp_muscle.afa'
        muscle_cmd = 'muscle -in %s -out %s'%(muscle_in, muscle_out)
        os.system(muscle_cmd)

        with open(muscle_out, 'r') as f:
            string = f.read()
        aln = {}
        aln[0] = [''.join(x.split('\n')[1:]) for x in string.split('\n>')]

        
    for char in aln[0][0]:
        if aln[0][1][count] != char:
            pos = len(aln[0][0][:count].replace('-',''))
            mutations.append('%s%d%s'%(char,pos, aln[0][1][count]))
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
    muts_left = OrderedDict(sorted({int(re.findall('[0-9]+', mut)[0]):mut[0] for mut in mutations if mut.endswith('-')}.items(), key = lambda a: a[0]) )
    muts_right = OrderedDict(sorted({int(re.findall('[0-9]+', mut)[0]):mut[-1] for mut in mutations if mut.startswith('-')}.items(), key = lambda a: a[0]) )
    insertions_left = get_insertions(muts_left)
    insertions_right = get_insertions(muts_right)
    r = defaultdict(str)
    l = defaultdict(str)
    for mut in mutations:
        if mut.startswith('-'):
            r[int(re.findall('[0-9]+', mut)[0])] += re.findall('[A-Z]+|\*', mut)[0]
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

def get_alignment_summary(tch_seq, seq_id, sequences):
    mutations = get_all_mutations(tch_seq, seq_id, sequences)
    mutations_2 = get_all_insertions(seq_id, mutations)    
    return get_point_mutations(mutations_2)

def align_sequences(cds_hits, ref_gene_seqs):
    rows = []
    rows_nucl = []
    count = 0

    to_break = False
    for tch_seq, d in cds_hits.items():

        seq1 = ref_gene_seqs[tch_seq]
        for seq_id, seq2 in d.items():

            if set(seq2) != {'A', 'C', 'G', 'T'}:
                continue
            count = 0

            AA1 = Seq(seq1, IUPAC.unambiguous_dna).translate(table = 'Bacterial')
            AA2 = Seq(seq2, IUPAC.unambiguous_dna).translate(table = 'Bacterial')

            if AA1[-1] == '*':
                AA1 = AA1[:-1]
            if AA2[-1] == '*':
                AA2 = AA2[:-1]

            sequences = {tch_seq:AA1, seq_id:AA2}
            point_mutations, deletions, insertions = get_alignment_summary(tch_seq, seq_id, sequences)

            if 0 in insertions.keys():
                seq2 = seq2[len(insertions[0])*3:]
                del insertions[0]
                AA2 = Seq(seq2, IUPAC.unambiguous_dna).translate(table = 'Bacterial')

            ins = {x:y for x,y in insertions.items() if len(y) > 10}
            dels = {x:y for x,y in deletions.items() if len(y) > 10}
            if len(ins) > 0:
                rows.append({'TCH - locus-tag': tch_seq, 'TCH - sequence':seq1, 'strain - locus-details':seq_id, 'strain - sequence':seq2,
                            'Mutation type':'Large insertion','Mutation details':' AND '.join(['%s%s'%(x,y) for x,y in ins.items()])
                            })
                count += 1
            if len(dels) > 0:
                rows.append({'TCH - locus-tag': tch_seq, 'TCH - sequence':seq1, 'strain - locus-details':seq_id, 'strain - sequence':seq2,
                            'Mutation type':'Large deletion','Mutation details':' AND '.join(['%s%s'%(x,y) for x,y in dels.items()])
                            })

                count += 1

            if '*' in str(AA2):
                stop_codons = {x.start() for x in re.finditer('\*', str(AA2))}
                if len(stop_codons) > 0 and min(stop_codons) < 0.80*len(AA2):
                    rows.append({'TCH - locus-tag': tch_seq, 'TCH - sequence':seq1, 'strain - locus-details':seq_id, 'strain - sequence':seq2,
                            'Mutation type':'Early stop codons','Mutation details':' AND '.join([str(x) for x in stop_codons])
                            })

                count += 1

            if AA2[0] not in ['M', 'V', 'L']:
                rows.append({'TCH - locus-tag': tch_seq, 'TCH - sequence':seq1, 'strain - locus-details':seq_id, 'strain - sequence':seq2,
                            'Mutation type':'Mutation in the start codon','Mutation details': AA2[0]
                            })

                count += 1
            if count > 0:
                sequences = {tch_seq:ref_gene_seqs[tch_seq], seq_id:seq2}
                na_point_mutations, na_deletions, na_insertions = get_alignment_summary(tch_seq, seq_id, sequences)
                rows_nucl.append({'Point mutations (nucl)': ' AND '.join(na_point_mutations), 'Insertions (nucl)':' AND '.join(['%s%s'%(x,y) for x,y in na_insertions.items()]),
                 'Deletions (nucl)':' AND '.join(['%s%s'%(x,y) for x,y in na_deletions.items()]), 'strain - locus-details': seq_id, 'TCH - locus-tag':tch_seq})

    na_res = pd.DataFrame(rows_nucl)
    aa_res = pd.DataFrame(rows)
    if len(aa_res) > 0:
        res = aa_res.merge(na_res, on = ['strain - locus-details','TCH - locus-tag'], how = 'left')
    else:
        res = pd.DataFrame()
    return res


def get_bbh(BR1_out, BR2_out, PID_threshold, output_file_name = '', ref_genes = {}, aln_threshold = 0 ):
    BR2 = pd.read_table(BR2_out, names= ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore'])
    BR2 = BR2.loc[(BR2['PID'] > PID_threshold) & (BR2['eVal'] < 0.000001)]     
    if ref_genes != {}:
        BR2['aln_perc'] = [100*BR2.loc[index, 'alnLength']/len(ref_genes[BR2.loc[index, 'subject']]) for index in BR2.index]
        BR2 = BR2.loc[(BR2['aln_perc'] > aln_threshold)] 

    BR2 = BR2.set_index('gene')

    BR1 = pd.read_table(BR1_out, names= ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore'])
    BR1 = BR1.loc[(BR1['PID'] > PID_threshold) & (BR1['eVal'] < 0.000001)]     
    if ref_genes != {}:
        BR1['aln_perc'] = [100*BR1.loc[index, 'alnLength']/len(ref_genes[BR1.loc[index, 'gene']]) for index in BR1.index]
        BR1 = BR1.loc[(BR1['aln_perc'] > aln_threshold)] 

    BR1 = BR1.set_index('subject')

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

def get_syntenic_homology(sabir, locus_of_interest, ft):
    
    
    to_exclude = False
    i = ft.loc[ft['locus_tag'] == locus_of_interest].index.tolist()[0]
    
    for j in range(1,100):
        if i+j > len(ft):
            i = len(ft)
        if ft.loc[i+j]['locus_tag'] in gid2matches[sabir].keys():
            m = gid2matches[sabir][ft.loc[i+j]['locus_tag']][0]         # keep in mind that there could be two matches as opposed to just one. For later development
            s1 = int(m.split('s')[-1].split(':')[0])
            e1 = int(m.split('e')[-1])
            contig1 = m.split('-s')[0] 
            
            if int(s1) < 50 or int(contig1.split('_')[3])-int(e1) < 50:
                to_exclude = True

            right_syntenic_homolog = m  
            right_syntenic_homolog_ref = ft.loc[i+j]['locus_tag']
            S1 = ft.loc[i+j]['start']
            E1 = ft.loc[i+j]['end']
            k1 = i -j
            break

    i = ft.loc[ft['locus_tag'] == locus_of_interest].index.tolist()[0]
    for j in range(1,100):
        if i-j < 0:
            i0 = i
            i = len(ft) - i0
        if ft.loc[i-j]['locus_tag'] in gid2matches[sabir].keys():
            m = gid2matches[sabir][ft.loc[i-j]['locus_tag']][0]         # keep in mind that there could be two matches as opposed to just one. For later development
            s2 = int(m.split('s')[-1].split(':')[0])
            e2 = int(m.split('e')[-1])
            contig2 = m.split('-s')[0] 

            if int(s2) < 50 or int(contig2.split('_')[3])-int(e2) < 50:
                to_exclude = True

            left_syntenic_homolog = m  
            left_syntenic_homolog_ref = ft.loc[i-j]['locus_tag']
            S2 = ft.loc[i-j]['start']
            E2 = ft.loc[i-j]['end']
            k2 = i -j
            break

    if contig1 != contig2:
        to_exclude = True
        
    distance_in_strain = min([abs(a-b) for a in [s1,e1] for b in [s2,e2]])
    distance_in_ref = min([abs(a-b) for a in [S1,E1] for b in [S2,E2]])

    row = {'Genome ID':sabir, 'Missing gene':locus_of_interest, 'right syntenic homolog':right_syntenic_homolog, 'left syntenic homolog':left_syntenic_homolog,
     'distance in ref':distance_in_ref, 'distance in strain':distance_in_strain, 'no. genes in between in ref': abs(k1 - k2 - 1), 'Excluded':to_exclude}
    
    return row