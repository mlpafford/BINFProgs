"""
Markov modules written by Mary Lenore Pafford

Includes the functions:
read_table(input file)
find_alpha(table of kmers, count-kmers format)
pseudocount(kmers, int pseudo value, int kmer order, alphabet)
make_prob(kmers_count)
get_kmers(sequence string)
counting(sequence file)
prefix_suffix(sequence string, order)
 
"""

import sys
import collections
import operator
import itertools
from math import log

import BINFmod

def read_table(in_file):
    """read_table expects a file with a table with two rows, sequence kmers followed by counts.
       It then puts them into a collect.Counter where count element can be used separately.
       This is useful for programs looking to work with count values, like coding cost."""
    kmers = collections.Counter()
    for line in in_file:
       split_line = line.split()
       kmers[split_line[0]] = int(split_line[1])
    return kmers, len(split_line[0])

def find_alpha(table_kmers):
    """This function finds the alphabet to plug into pseudocount.
       only the letters found in the file are needed, 
       or else probabilities will be created for unexpected kmers"""
    alpha = ""
    for kmer in table_kmers.items():
        for char in kmer[0]:
            if char not in alpha:
                alpha += char
    return alpha

def pseudocount(original_kmers, pseudo, order, alpha):
    """pseudocount takes the given alphabet and creates pseoudocounts for all kmers,
       including the unseen combinations that could exist in a test file"""
    counts = collections.Counter()
    prefix = 1
    suffix = 1
    kmers = itertools.product(alpha, repeat = order)
    for value in kmers:
        # tuple format of kmer strings
        kmer_tuple = str(value).translate(None,"(),'")
        kmer = ''
        for char in kmer_tuple.split():
            kmer += char
        counts[kmer] = pseudo
    if order == 0:
        # needs counts for start and stop
        counts["$"] = pseudo
        counts["*"] = pseudo

    for i in range(order):
        if order > 0:
        # add in any kmer combos with the prefix character
            for item in itertools.product(alpha, repeat = order - prefix):
                kmer_tuple = str(item).translate(None,"(),'")
                kmer = ""
                for char in kmer_tuple.split():
                    kmer += char
                prefix_string = ""
                for int in range(prefix):
                    prefix_string += "$"
                prefix_string += kmer
                counts[prefix_string] = pseudo
            prefix += 1
    suffix = 1
    for i in range(order):
        if order > 1:
        # add in any kmer combos with the suffix character
            for item in itertools.product(alpha, repeat = order - suffix):
                kmer_tuple = str(item).translate(None,"(),'")
                kmer = ""
                for char in kmer_tuple.split():
                    kmer += char
                suffix_string = ""
                for int in range(suffix):
                    suffix_string += "*"
                kmer += suffix_string
                counts[kmer] = pseudo
            suffix += 1
            for char in alpha:
                if len(kmer) == 1:
                    kmer += char
                    counts[kmer] = 0
                kmer = ''

    # Merge current and pseudocount  
    return original_kmers + counts

def make_prob(kmers_count):
    """This function takes counts of kmers, formatted as tuple of kmer and count.
       It creates probabilities for each. Requires pseudocounts already included to avoid zeros"""
    prob_dict = dict()
    num_kmers = 0
    #will need group partial kmers, (k-1)mer, by sorting. Group of those (k-1)mer
    kmers = sorted(kmers_count)
    group = None
    for itor, kmer in enumerate(kmers):
        current_k = kmers[0:-1]
        #start new (k-1)mer group
        if group != current_k:
            group = current_k
            for selected_kmer in kmers[itor-num_kmers:itor]:
                prob = float(kmers_count[selected_kmer])/float(sum_total)
                #want to store conditional log2 of probablities
                prob_dict[selected_kmer] = -log(prob,2)
            #reset for next group
            sum_total = 0
            num_kmers = 0
        num_kmers += 1
        sum_total+= kmers_count[kmer]
    itor += 1
    #must include final group
    for selected_kmer in kmers[itor-num_kmers:itor]:
        prob = float(kmers_count[selected_kmer])/float(sum_total)
        prob_dict[selected_kmer]= -log(prob,2)

    return prob_dict

def get_kmers(seqs,k):
        """A way to split a sequence in to kmers of given size
           Yields one kmer at a time.
        """
        for start in xrange(len(seqs) - k + 1):
                yield seqs[start:start+k]


def counting(counts, seq_file, kvalue, alpha):
    """This function both splits and counts kmers, given the file and an alphabet.
       It returns the counter passed filled with kmer keys and count values"""
    for fasta_seq in BINFmod.read_fasta(seq_file, alpha):
        sequence = fasta_seq.seq
        #to take care of empty case
        if len(sequence)== 0: continue
        #sequence needs prefix '^' and suffix '$'
        #added using prefix_suffix()
        seq = prefix_suffix(sequence, kvalue)
        for start in range(len(seq)-kvalue +1):
            counts[seq[start:start+kvalue]] += 1
    return(counts)


def prefix_suffix(sequence, order):
    """This function exists to put prefixs and suffixs on the sequence.
       This allows the program to count (and later calculate probability of)
       a letter followed by stop, or a start character with a following letter.
       Important if used for coding-cost program"""
    prefix = "^"
    suffix = "$"
    if order == 1:
    #Order 0 does not need the prefix character
    #Probability of start character alone is one
        prefix = ""
    for i in range(order-2): prefix = prefix + "^"
    for i in range(order-2): suffix = suffix + "$"
    sequence = prefix + sequence + suffix
    #for symmetry there needs to be same number of start and stops
    return(sequence)
