#!/usr/bin/env python

"""
Script for raw (only blast based) protein predictions.
recommended  blast options: -outfmt 5 -gapopen 11 -gapextend 2
-dust no -soft_masking
"""

import re
import argparse
from multiprocessing import Pool
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
import itertools
from Bio.pairwise2 import align
from Bio.SubsMat import MatrixInfo as matlist
from tqdm import tqdm
from time import time


def timer(original_function):
    def wrapper_func(*args, **kwargs):
        start = time()
        func = original_function(*args, **kwargs)
        end = time()
        elapsed = end - start
        print(elapsed)
        return func
    return wrapper_func


trans_table = {
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
    genomic_dict = {}
    for contig in SeqIO.parse(genome, 'fasta'):
        genomic_dict[contig.name] = Seq(str(contig.seq).upper())
    return genomic_dict


def translation(sequence, codons=None, table1=None):
    """ Simple translation. X for ambiguous codons.
    Input: nucleotide sequence
    Returns: translated nucleotide sequence (str) """
    tr_table = trans_table
    if table1:
        tr_table = trans_table.copy()
        tr_table["TAA"] = "*"
        tr_table["TAG"] = "*"
        tr_table["TGA"] = "*"
    if not codons:
        codons = [sequence[i:i + 3] for i in range(0, len(sequence) - 2, 3)]
    amino_acids = []
    for codon in codons:
        if codon in tr_table:
            amino_acids.append(tr_table[codon])
        else:
            amino_acids.append('X')
    return ''.join(amino_acids)


class Gene:
    """ Gene class. Mainly for best blast hit extraction. """

    def __init__(self, name):
        self.name = name
        self.blast_hits = []

    def __repr__(self):
        return self.name

    def add_blast(self, blast_record):
        """ Adds blast record to the self.blast_hits list """
        self.blast_hits.append(blast_record)

    def get_best_blast(self):
        """ Return best blast hit for given gene from different organisms.
        Input: self
        Returns: object of  Blast.Record class (Biopython)"""
        gene = self.name
        best_blast = None
        best_blast_e = 0
        for blast in gene_dict[gene].blast_hits:
            evalue = float(blast.alignments[0].hsps[0].expect)
            if not best_blast:
                best_blast = blast
                best_blast_e = evalue
            else:
                if evalue < best_blast_e:
                    best_blast_e = evalue
                    best_blast = blast
        return best_blast


def correct_orientation(contig_name, hit_frame):
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
        real_start = (start - h_sequence[:start].count('-')) * 3
        insertions.append((real_start + h_start,
                           real_start + h_start + (end - start) * 3))
    return insertions


def potential_inframe_intron(intron_position, contig_seq):
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
    if base_end + 10 <= len(contig_seq):
        end = base_end + 10
    else:
        end = base_end
    target = str(contig_seq[start:end])
    gt = [m.start() for m in re.finditer('GT', target[:int((len(target)) / 2)])]
    ag = [m.end() for m in re.finditer('AG', target)]
    potential_introns = []
    for r in itertools.product(gt, ag):
        # combinations of all possible GT....AG sequences
        if r[0] < r[1]:
            intron = target[r[0]:r[1]]
            if len(intron) > 4 and len(intron) % 3 == 0:
                potential_introns.append(str(intron))
    potential_introns.append('')
    return potential_introns


def in_frame_introns(hsp, contig_name, h_len):
    """Uses potential_inframe_intron to gives you
    all in-frame introns from one hsp. It tries all possible intron combination
    and return the best combination accordint to distance:
    abs(query without dashes + dashes where introns was not found) - len of
    new protein.
    Input: hsp Bio.Blast object, contig name, length of hit (contig)
    Returns: list of best intron combination e.g. ['GTAGGAAG', '', 'GTAAAG]"""
    h_frame = hsp.frame[1]
    h_sequence = str(hsp.sbjct)
    h_start = hsp.sbjct_start
    h_end = hsp.sbjct_end
    if h_frame < 0:
        fake_start = h_start
        h_start = h_len - h_end + 1
        h_end = h_len - fake_start + 1
    q_sequence = str(hsp.query)
    # all insertions gives you list of tuples e.g. [(0,9), (30-90)]
    intron_positions = get_insertions(q_sequence, h_sequence, h_start)

    distance_dict = {}  # important for distance counting
    for count, position in enumerate(intron_positions):
        distance_dict[count] = int((position[1] - position[0]) / 3)

    potential_introns = []
    intron_sequences = []
    if intron_positions:
        # gives you hit sequence in correct orientation
        contig_seq = correct_orientation(contig_name, h_frame)
        intron_sequences = []
        # gives you intron sequence in nucleotides
        for intron_position in intron_positions:
            intron_seq = potential_inframe_intron(intron_position, contig_seq)
            potential_introns.append(intron_seq)

        lower_distance = 500  # arbitrary used random high number
        # test combinations of potential introns
        for combination in itertools.product(*potential_introns):
            test_seq = str(contig_seq[h_start - 1:h_end])
            for i in combination:
                test_seq = test_seq.replace(i, '')
            codons_ = [test_seq[i:i + 3] for i in range(0, len(test_seq) - 2, 3)]
            if len(stops.difference(set(codons_))) == len(stops):
                protein = translation(test_seq, codons=codons_)
                query_length = len(str(hsp.query).replace('-', ''))
                for count, value in enumerate(combination):
                    if not value:
                        query_length += distance_dict[count]
                distance = abs(query_length - len(protein))
                if distance < lower_distance:
                    lower_distance = distance
                    intron_sequences = list(combination)
    return intron_sequences


def ntrons_all_hsps(sample, hit_num):
    """gives you all in-frame introns for all hsps.
    Input: blast class from Biopython
    Returns: list of best in-frame intron candidates """
    contig_name = sample.alignments[hit_num].hit_id
    if '|' in contig_name:
        contig_name = contig_name.split('|')[1]
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


def hsps_coordinates(sample, hit_num):
    """Gives you sorted and corrected (according to similarity/identity)
    coordinates off all hsps.
    Input: sample - blast class from Biopython, hit number
    Returns: tuple - coordinates, global_start, global_end
    coordinates = coordinates of hsps used for building protein
    global_start = lowest contig coordinate
    global_end = highest contig coordinate"""
    global_start = 50000000000000000000
    global_end = 0
    pseudo_coordinates = {}
    coordinates = []
    query_coordinates_unsorted = []
    starting_set = set()
    hsps_range = set()
    for hsp in sample.alignments[hit_num].hsps:
        h_frame = hsp.frame[1]
        h_start = hsp.sbjct_start
        h_end = hsp.sbjct_end
        h_len = sample.alignments[hit_num].length
        middline = hsp.match
        q_start = hsp.query_start
        q_end = hsp.query_end

        if h_frame < 0:
            fake_start = h_start
            h_start = h_len - h_end + 1
            h_end = h_len - fake_start + 1

        close = False

        if not hsps_range:
            hsps_range.add(h_start)
            hsps_range.add(h_end)
            close = True
        else:
            for hsp in hsps_range:
                if abs(hsp - h_start) < 15000 or abs(hsp - h_end) < 15000:
                    close = True

        if close:
            if h_start < global_start:
                global_start = h_start
            if h_end > global_end:
                global_end = h_end

            # get rid of hsps which are located in the same place on CONTIG
            x = set(range(h_start, h_end))
            if len(starting_set.intersection(x)) < 15:
                hsps_range.add(h_start)
                hsps_range.add(h_end)
                starting_set.update(x)
                pseudo_coordinates[q_start] = (
                    [h_start, h_end, str(middline)])
                query_coordinates_unsorted.append([q_start, q_end])

    pseudo_coordinates_sorted = sorted(pseudo_coordinates)
    query_coordinates = sorted(query_coordinates_unsorted)
    n = 0
    iter_dict = {}
    for count, item in enumerate(pseudo_coordinates_sorted):
        iter_dict[count] = item
        n = count

    for i in range(len(pseudo_coordinates) - 1):
        # set lower coordinate in contig as a first coordinate
        if i == 0:
            coordinates.append(pseudo_coordinates[iter_dict[i]][0])

        x = set(range(query_coordinates[i][0], query_coordinates[i][1]))
        y = set(range(query_coordinates[i + 1][0], query_coordinates[i + 1][1]))
        # if there is overlap of hsps, this should connect overlap to more
        # similar hsp
        if x.intersection(y):
            # minus first coordinate !
            xseq = (str(pseudo_coordinates[iter_dict[i]][2])
            [-len(x.intersection(y)):])
            yseq = (str(pseudo_coordinates[iter_dict[i + 1]][2])
            [:len(x.intersection(y))])

            # CHECK WHAT GIVES YOU BETTER RESULTS
            # SIMILARITY X IDENTITY
            identx = int(len(xseq) - xseq.count(" ") - xseq.count("+"))
            identy = int(len(yseq) - yseq.count(" ") - yseq.count("+"))

            #            identx = int(len(xseq) - xseq.count(" "))
            #            identy = int(len(yseq) - yseq.count(" "))

            if identx > identy:
                coordinates.append(pseudo_coordinates[iter_dict[i]][1])
                coordinates.append(
                    pseudo_coordinates[iter_dict[i + 1]][0] + (len(xseq) * 3))
            else:
                coordinates.append(
                    pseudo_coordinates[iter_dict[i]][1] - (len(xseq) * 3))
                coordinates.append(pseudo_coordinates[iter_dict[i + 1]][0])
        else:
            coordinates.append(pseudo_coordinates[iter_dict[i]][1])
            coordinates.append(pseudo_coordinates[iter_dict[i + 1]][0])
    coordinates.append(pseudo_coordinates[iter_dict[n]][1])
    coordinates.sort()
    return coordinates, global_start, global_end


def give_inter_intron(starting_point, ending_point, contig_seq):
    """Gives you set of potential introns between two hsps.
    Input: starting_point - last expected position of first exon,
            ending_point - first expected position of  second exon
    Returns: list of potential intron sequences between two exons"""
    if (starting_point - 15) >= 0:
        start = starting_point - 15
    else:
        start = 0
    if (ending_point + 15) <= len(contig_seq):
        end = ending_point + 15
    else:
        end = ending_point
    target = str(contig_seq[start:end])
    gt = [m.start() for m in re.finditer('GT',
                                         str(contig_seq[start:start + 30]))]
    ag = [m.end() for m in re.finditer('AG', str(contig_seq[end - 30:end]))]
    for i in range(len(ag)):
        ag[i] += len(contig_seq[start:end - 30])
    potential_introns = []
    for r in itertools.product(gt, ag):
        if r[0] < r[1]:
            intron = target[r[0]:r[1]]
            if len(intron) > 4:
                potential_introns.append(intron)
    return potential_introns


def get_all_inter_introns(coordinates, contig_seq):
    """Return list of all possible inter introns.
    Input: list of coordinates of hsps used for building protein,
            contig sequence in nucleotides
    Returns: list of lists - every item in first level list represents one
            particular position between two exons and it is list of all
            potential intron sequences in nucleotides"""
    all_inter_introns_list = []
    inter_introns = int((len(coordinates) - 2) / 2)
    index1 = 1
    index2 = 2
    for _ in range(inter_introns):
        all_inter_introns_list.append(give_inter_intron(coordinates[index1],
                                                        coordinates[index2],
                                                        contig_seq))
        index1 += 2
        index2 += 2
    return all_inter_introns_list


def check_best_prediction(prot_sequences, query_name):
    """Chooses best protein prediction
    Input: protein candidates (list of protein sequenes)
    Returns: best protein prediction according to blast"""
    best_prediction = ""
    less_dashes = 10000
    matrix = matlist.blosum62
    gap_open = -11
    gap_extend = -1
    for sequence in prot_sequences:
        ga = align.globalds(query_dataset[query_name], sequence,
                            matrix, gap_open, gap_extend)[0]
        q = (ga[0].count('-'))
        h = (ga[1].count('-'))
        dashes = q + h
        if dashes < less_dashes:
            less_dashes = dashes
            best_prediction = ga[1]
    return best_prediction.replace('-', '')


def best_hsp_seq(sample, contig_seq):
    hsp = sample.alignments[0].hsps[0]
    h_len = sample.alignments[0].length
    h_start = hsp.sbjct_start
    h_end = hsp.sbjct_end
    if hsp.frame[1] < 0:
        fake_start = h_start
        h_start = h_len - h_end 
        h_end = h_len - fake_start + 1
    return contig_seq[h_start:h_end]

def hsps_prot_seq(sample):
    result = ""
    for hsp in sample.alignments[0].hsps:
        result += hsp.sbjct
    return result.replace('-', '')

def protein_prediction(sample, hit_num):
    """Return best protein predictions.
    Input: genome dict, sample - blast class from Biopython, hit_number
    Returns: best protein prediction (protein sequence) """
    query_name = sample.query
    contig_name = sample.alignments[hit_num].hit_id.split()[0]
    if '|' in contig_name:
        contig_name = contig_name.split('|')[1]
    _, h_frame = sample.alignments[hit_num].hsps[0].frame
    contig_seq = correct_orientation(contig_name, h_frame)

    coordinates, global_start, global_end = hsps_coordinates(sample,
                                                             hit_num)
    result_sequence = str(contig_seq[global_start - 1:global_end + 1]) #tady je to jeste dobre
    all_introns_no_none = ntrons_all_hsps(sample, hit_num)


    for sequence in all_introns_no_none:
        result_sequence = result_sequence.replace(sequence, '')


    if len(sample.alignments[0].hsps) > 1:
        all_list = get_all_inter_introns(coordinates, contig_seq)
        all_list_no_none = [l for l in all_list if len(l) > 0]
        for list_ in all_list_no_none:
            list_.append('')
        # append('') means that we are adding no intron at the place possibility
        best_candidates = []
        all_combination = 1
        for combination in itertools.product(*all_list_no_none):
            all_combination += 1
            if all_combination > 150000:
                break
            else:
                test_seq = result_sequence[:]
                for i in combination:
                    test_seq = test_seq.replace(i, '')
                codons_ = [test_seq[i:i + 3] for i in range(0, len(test_seq) - 2, 3)]
                if len(stops.difference(set(codons_))) == len(stops):
                    all_hsps_seqs = hsps_prot_seq(sample)
                    protein_t1 = translation(test_seq, codons=codons_, table1=True)
                    four_mers = [protein_t1[i:i + 4] for i in range(0, len(protein_t1) - 3, 4)]
                    mismatches = 0
                    for four_mer in four_mers:
                        if four_mer not in all_hsps_seqs:
                            mismatches += 1
                    if len(protein_t1)*0.03 >= mismatches:
                        protein = translation(test_seq, codons=codons_)
                        best_candidates.append(protein)
        best = check_best_prediction(best_candidates, query_name)
        if best:
            return best
        else:
            best_hsp = str(best_hsp_seq(sample, contig_seq))
            for sequence in all_introns_no_none:
                best_hsp = best_hsp.replace(sequence, '')
            return translation(best_hsp)
    else:
        return translation(result_sequence)


def final(gene):
    print(gene)
    samples = gene_dict[gene].blast_hits
    contig_dict = {}
    result = []
    for sample in samples:

        for hit_number in range(hits):
            if len(sample.alignments) > hit_number:
                contig = sample.alignments[hit_number].hit_id
                seq = protein_prediction(sample, hit_number)
                if contig not in contig_dict and '*' not in seq:
                    contig_dict[contig] = seq
                elif '*' in seq:
                    pass # ???
                else:
                    if len(seq) > len(contig_dict[contig]) and '*' not in seq:
                        contig_dict[contig] = seq
    for contig, protein in contig_dict.items():
        result.append('>{}_{}@{}\n{}\n'.format(org_name, gene,
                                               contig, protein))
    return result


@timer
def main():
    for blast_record in NCBIXML.parse(blastout):
        gene_name = blast_record.query.split('_')[-1]
        if len(blast_record.alignments) > 0:
            if gene_name not in gene_dict:
                gene_dict[gene_name] = Gene(gene_name)
            gene_dict[gene_name].add_blast(blast_record)

    with open(args.output, 'w') as res:
        # with Pool(processes=threads) as p:
        #     max_ = len(gene_dict)
        #     r = list(tqdm(p.imap(final, gene_dict), total=max_))


        # debugging
        r = []
        for i in gene_dict:
            r.append(final(i))
        # ########

        for record in r:
            if len(record) > 0:
                for i in record:
                    res.write(i)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Raw gene prediction tool')
    parser.add_argument("protein_dataset")
    parser.add_argument("genome")
    parser.add_argument("blast_output")
    parser.add_argument("organism_name")
    parser.add_argument("output")
    parser.add_argument("--genetic-code",
                        help='altrnative genetic code in a form: "TAA:Q;TAG:Q" ')
    parser.add_argument("--threads", type=int)
    parser.add_argument("--hits", type=int, help='number of hits for one protein')
    args = parser.parse_args()

    stops = {'TAA', 'TAG','TGA'}

    if args.genetic_code:
        if ';' in args.genetic_code:
            alt_code = args.genetic_code.split(';')
            for i in alt_code:
                splitted = i.split(':')
                trans_table[splitted[0]] = splitted[1]
                if splitted[0] in stops:
                    stops.remove(splitted[0])
        else:
            splitted = args.genetic_code.split(':')
            trans_table[splitted[0]] = splitted[1]
            if splitted[0] in stops:
                stops.remove(splitted[0])

    threads = 1
    if args.threads:
        threads = args.threads

    hits = 1
    if args.hits:
        hits = args.hits

    query_dataset = read_proteins(args.protein_dataset)
    genome_dict = read_genome(args.genome)
    org_name = args.organism_name
    gene_dict = {}
    blastout = open(args.blast_output)
    main()
