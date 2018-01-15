import re
import os
import argparse
from Bio.Blast import NCBIXML
from Bio import SeqIO
import itertools
from Bio.pairwise2 import align
from Bio.SubsMat import MatrixInfo as matlist


def read_proteins(dataset):
    """ Reads protein dataset used for blasting as a
    dictionary.
    Input: protein sequences in fasta format.
    Returns: {prot_name:seq} """
    dataset_dict = {}
    for protein in SeqIO.parse(dataset, 'fasta'):
        dataset_dict[protein.description] = str(protein.seq)
    return dataset_dict


def read_genome(genome):
    """ Reads genome as a dictionary.
    Input: nucleotide sequences in fasta format
    Returns: {contig_name:seq} """
    genome_dict = {}
    for contig in SeqIO.parse(genome, 'fasta'):
        genome_dict[contig.name] = contig.seq
    return genome_dict


def get_translation(sequence):
    """ Simple translation. X for ambiguous codons.
    Input: nucleotide sequence
    Returns: translated nucleotide sequence (str) """
    translation_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        }
    cutting = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
    amino_acids = []
    for codon in cutting:
        try:
            amino_acids.append(translation_table[codon])
        except KeyError:
            amino_acids.append('X')
    return ''.join(amino_acids)


class Gene():
    """ Gene class. Mainly for best blast hit extraction. """

    def __init__(self, name):
        self.name = name
        self.best_blast_hits = []

    def __repr__(self):
        return self.name

    def add_blast(self, blast_record):
        """ Adds blast record to the self.best_blast_hits list """
        self.best_blast_hits.append(blast_record)

    def get_best_blast(self):
        """ Return best blast hit for given gene from different organisms.
        Input: self
        Returns: object of  Blast.Record class (Biopython)"""
        gene = self.name
        best_blast = None
        best_blast_e = 0
        for blast in gene_dict[gene].best_blast_hits:
            evalue = float(blast.alignments[0].hsps[0].expect)
            if not best_blast:
                best_blast = blast
                best_blast_e = evalue
            else:
                if evalue < best_blast_e:
                    best_blast_e = evalue
                    best_blast = blast
        return best_blast


def get_correct_orientation(contig_name, hit_frame):
    """ Gives you contig in proper orientation.
    Input: contig name, frame of contig sequence according to hit
    Returns: contig sequence in proper orientation (str) """
    if hit_frame < 0:
        contig_seq = genome_dict[contig_name].reverse_complement().upper()
    else:
        contig_seq = genome_dict[contig_name].upper()
    return contig_seq


def get_insertions(q_sequence, h_sequence, h_start):
    """Gets insertions of hit sequence
    Input: query sequence (str), hit sequence (str), hit start (int)
    Returns: list of insertions """
    insertions = []
    for m in re.finditer(r"-{3,100}", q_sequence):
        # looking for 'deletions' in query sequence in alignment
        start, end = [m.start(), m.end()]
        real_start = (start - h_sequence[:start].count('-'))*3
        insertions.append([real_start+h_start,
                           real_start+h_start+(end-start)*3])
    return insertions


def find_inframe_intron(intron_position, contig_seq):
    """Gives you potentional in-frame intron sequence.
    Input: potentional intron position (tuple of two int, start and stop),
    contig seq (str)
    Returns: list of potential introns (nucleotide sequence; str)"""
    base_start = intron_position[0]
    base_end = intron_position[1]
    # base start and base end are set up to make bigger sequence to
    # look for GT and AG (standard intron stard and end, respectively)
    if base_start - 10 >= 0:
        start = base_start - 10
    else:
        start = base_start
    if base_end + 10 <= len(contig):
        end = base_end + 10
    else:
        end = base_end
    target = str(contig_seq[start:end])
    GT = [m.start() for m in re.finditer('GT', target[:int((len(target))/2)])]
    AG = [m.end() for m in re.finditer('AG', target)]
    potential_introns = []
    for r in itertools.product(GT, AG):
        # combinations of all possible GT....AG sequences
        if r[0] < r[1]:
            intron = target[r[0]:r[1]]
            if (len(intron) > 4 and len(intron) % 3 == 0):
                potential_introns.append(str(intron))
    potential_introns.append('')
    return potential_introns


def in_frame_introns(hsp, contig_name, h_len):
    """Uses find_inframe_intron to gives you all in-frame introns from one hsp."""
    q_frame, h_frame = hsp.frame
    h_sequence = str(hsp.sbjct)
    h_start = hsp.sbjct_start
    h_end = hsp.sbjct_end
    if h_frame < 0:
        fake_start = h_start
        h_start = h_len - h_end + 1
        h_end = h_len - fake_start + 1
    q_sequence = str(hsp.query)
    intron_positions = get_insertions(q_sequence, h_sequence, h_start)

    potential_introns = []
    intron_sequences = []
    if intron_positions:
        contig_seq = get_correct_orientation(contig_name, h_frame)
        intron_sequences = []
        for intron_position in intron_positions:
            intron_seq = find_inframe_intron(intron_position, contig_seq)
            potential_introns.append(intron_seq)

        lower_distance = 500  # arbitrary used random high number
        for combination in list(itertools.product(*potential_introns)):
                test_seq = str(contig_seq[h_start-1:h_end])
                for i in combination:
                    test_seq = test_seq.replace(i, '')
                protein = get_translation(test_seq)
                if '*' not in protein:
                    distance = abs((len(str(hsp.query).replace('-', ''))
                                    - len(protein)))
                    if distance < lower_distance:
                        intron_sequences = list(combination)
    return intron_sequences


def get_blast_iteration(num, database):
    """Gives you only one blast iteration specified by number.
    Only for testing."""
    n = 0
    while n < num:
        sample = next(database)
        n += 1
    return sample


def get_introns_all_hsps(sample, hit_num):
    """gives you all in-frame introns for all hsps."""
    contig_name = sample.alignments[hit_num].hit_id
    h_len = sample.alignments[hit_num].length
    all_introns = []
    for hsp in sample.alignments[hit_num].hsps:
        intron_seqs = in_frame_introns(hsp, contig_name, h_len)
        if intron_seqs:
            all_introns += intron_seqs

    all_introns_no_none = []
    for i in all_introns:
        if i != 'None':
            all_introns_no_none.append(i)

    return all_introns_no_none


def similarity(seq):
    space = seq.count(" ")
    plus = seq.count("+")
    result = len(seq)-(space+plus)
    return result


def get_hsps_coordinates(sample, hit_num):
    """Gives you sorted (increasing) coordinates off all hsps"""
    global_start = 500000
    global_end = 0
    pseudo_coordinates = {}
    coordinates = []
    query_coordinates_unsorted = []
    for hsp in sample.alignments[hit_num].hsps:
        q_frame, h_frame = hsp.frame
        h_start = hsp.sbjct_start
        h_end = hsp.sbjct_end
        h_len = sample.alignments[hit_num].length
        middline = hsp.match
        q_start = hsp.query_start
        q_end = hsp.query_end
        q_frame, h_frame = hsp.frame
        if h_frame < 0:
            fake_start = h_start
            h_start = h_len - h_end + 1
            h_end = h_len - fake_start + 1
        if h_start < global_start:
            global_start = h_start
        if h_end > global_end:
            global_end = h_end
        pseudo_coordinates[q_start] = [int(h_start), int(h_end), str(middline)]
        query_coordinates_unsorted.append([q_start, q_end])

    query_coordinates = sorted(query_coordinates_unsorted)
    pseudo_coordinates_sorted = sorted(pseudo_coordinates)

    iter_dict = {}
    n = 0
    for i in pseudo_coordinates_sorted:
        iter_dict[n] = i
        n += 1

    for i in range(len(pseudo_coordinates)-1):
        if i == 0:

            coordinates.append(pseudo_coordinates[iter_dict[i]][0])

        x = set(range(query_coordinates[i][0], query_coordinates[i][1]))
        y = set(range(query_coordinates[i+1][0], query_coordinates[i+1][1]))

        if x.intersection(y):
            # minus first coordinate !
            xseq = (str(pseudo_coordinates[iter_dict[i]][2])
                    [-len(x.intersection(y)):])
            yseq = (str(pseudo_coordinates[iter_dict[i+1]][2])
                    [:len(x.intersection(y))])

            simx = int(similarity(xseq))
            simy = int(similarity(yseq))

            if simx > simy:
                coordinates.append(pseudo_coordinates[iter_dict[i]][1])
                coordinates.append(pseudo_coordinates[iter_dict[i+1]][0]+(len(xseq)*3))
            else:
                coordinates.append(pseudo_coordinates[iter_dict[i]][1]-(len(xseq)*3))
                coordinates.append(pseudo_coordinates[iter_dict[i+1]][0])
        else:
            coordinates.append(pseudo_coordinates[iter_dict[i]][1])
            coordinates.append(pseudo_coordinates[iter_dict[i+1]][0])
    coordinates.append(pseudo_coordinates[iter_dict[n-1]][1])
    coordinates.sort()
    return coordinates, global_start, global_end

def give_inter_intron(coordinates,index1, index2, contig_seq):
    """Gives you set of potential introns between two hsps"""
    
    try:
        start = coordinates[index1] - 10
    except:
        start = 0
    try:
        end = coordinates[index2] + 10
    except:
        end = coordinates[index2]
    target = str(contig_seq[start:end])
    GT = [m.start() for m in re.finditer('GT', str(contig_seq[start:start+20]))]
    AG = [m.end() for m in re.finditer('AG', str(contig_seq[end-20:end]))]
    for i in range(len(AG)):
        AG[i] += len(contig_seq[start:end-20])
    potential_introns = []
    for r in itertools.product(GT, AG):
        if r[0] < r[1]:
            intron = target[r[0]:r[1]]
            if len(intron) > 4:
                potential_introns.append(intron)
    return potential_introns

def get_all_inter_introns(coordinates, contig_seq):
    """Return list of all possible inter introns."""
    all_inter_introns_list = []
    inter_introns = int((len(coordinates) - 2)/2)
    index1 = 1
    index2 = 2
    for i in range(inter_introns):
        all_inter_introns_list.append(give_inter_intron(coordinates,index1, 
                                                        index2,contig_seq))
        index1 += 2
        index2 += 2
    return all_inter_introns_list

def check_best_prediction(prot_sequences, query_name):
    """Chooses best inter intron, and return CDS"""
    best_prediction = ""
    less_dashes = 10000
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -0.5
    for sequence in prot_sequences:
        for a in align.localds(query_dataset[query_name], sequence, 
            matrix, gap_open, gap_extend):
            q = (a[0].count('-'))
            h = (a[1].count('-'))
            dashes = q+h
            if dashes < less_dashes:
                less_dashes = dashes
                best_prediction = a[1]
    return best_prediction.replace('-', '')


def get_protein_predictions(genome_dict, sample, hit_num):
    query_name = sample.query
    """Return best protein predictions."""
    contig_name = sample.alignments[hit_num].hit_id
    q_frame, h_frame = sample.alignments[hit_num].hsps[0].frame
    contig_seq = get_correct_orientation(contig_name, h_frame)

    coordinates, global_start, global_end = get_hsps_coordinates(sample, hit_num)

    result_sequence = str(contig_seq[global_start-1:global_end+1])
    all_introns_no_none = get_introns_all_hsps(sample, hit_num)
    for sequence in all_introns_no_none:
        result_sequence = result_sequence.replace(sequence, '')
#        print(sequence)

    if len(sample.alignments[0].hsps) > 1:
        all_list = get_all_inter_introns(coordinates, contig_seq)
        all_list_no_none = [l for l in all_list if len(l) > 0]
        best_candidates = []
        for combination in list(itertools.product(*all_list_no_none)):
            test_seq = result_sequence[:]
            for i in combination:
                test_seq = test_seq.replace(i, '')
                protein = get_translation(test_seq)
                if '*' not in protein:
                    best_candidates.append(protein)
#                    for i in combination:
#                        print(i)
                
        best = check_best_prediction(best_candidates, query_name)
        if best:
            return best
        else:
             return(get_translation(result_sequence))
    else:
        return(get_translation(result_sequence))

    
query_dataset = read_proteins('TEST_ABCE')
genome_dict = read_genome('Metopus_genome_R1R21.fasta') 
org_name = "Metopus"
gene_dict = {}

blastout = open('Metopus_TEST.xml')
for blast_record in NCBIXML.parse(blastout):
    gene_name = blast_record.query.split('_')[-1]
    try:
        test = blast_record.alignments[0].hsps[0]
        if gene_name not in gene_dict:
            gene_dict[gene_name] = Gene(gene_name)
        gene_dict[gene_name].add_blast(blast_record)
    except:
        pass



with open('Metopus_predictions_TEST.fasta', 'w') as res:
    for gene in gene_dict:
        samples = gene_dict[gene].best_blast_hits
        contig_dict = {}
        for sample in samples:
            contig = sample.alignments[0].hit_id
            seq = get_protein_predictions(genome_dict, sample, 0)
            evalue = sample.alignments[0].hsps[0].expect
            if contig not in contig_dict and '*' not in seq:
                contig_dict[contig] = seq
                min_eval = evalue
            elif '*' in seq:
                pass
            else:
                if len(seq) > len(contig_dict[contig]) and '*' not in seq:
                    contig_dict[contig] = seq
        for contig, protein in contig_dict.items():
            res.write('>{}_{}@{}\t{}\n{}\n'.format(org_name, gene, 
                      sample.query, contig, protein))