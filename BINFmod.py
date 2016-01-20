#!/usr/bin/env python2.7

from __future__ import print_function

"""
    Bioinformatics Modules by Mary Lenore Pafford, BME205 asg2

    This file contains 3 modules to be used in various bioinformatics programs.
    read_fasta takes in a FASTA file and returns the idenifier, comments and sequence as strings.
    read_fastq takes a FASTQ file and returns the identifier and sequence as strings.
    read_fasta_with_quality takes in a quality file as well as its quality file.
    If using quality, you must run read_fasta for the actual sequence.

    This module will help eliminate white space and unusable characters
    so that there will not be errors when running the sequences through a pipeline.

    The fast class:
        iden = string identifier, required
        comment = string, optional
        seq = string of bases
        quality = [integer scores for sequence bases]

   This class is designed to be used for both FASTA and FASTQ parsing.
   It is important to note that FASTQ scores are not originally in integer scores, therefore must be unpacked using "ord()" and the  
score
"""

import sys
import argparse
import operator
import collections
import fileinput
import re

class fast:
    """ This class contains the string values for identity, comment and bases of an individual sequence.
        Quality scores of the sequence are kept in a list, maintaining order of corresponding bases"""
    iden = None
    comment = None
    seq = None
    qual = []


def read_fasta(file, alphabet):
    """ Reads FASTA file and splits based on iden, comment and sequence
        Yeilds all as individenual string elements in the fast class """
    fsta = fast()
    for line in file:
        #first step takes any iden line and keeps it out of the sequence string
        if line.startswith('>'):
            if fsta.iden is not None: yield(fsta)
            sline = re.search("[\s]", line)
            fsta.iden = line[1: sline.start(0)]
            fsta.comment = line[sline.start(0): -1]
            fsta.seq_list = []
            fsta.seq = ""
            continue
        else:
            #capitalize and check alphabet
            line = line.rstrip()
            line = line.upper()
            for letter in line:
                if letter in alphabet:
                    fsta.seq += letter
    yield fsta


def read_fastq(in_file, number):
    """ Reads in FASTQ file and splits based on idenifying lines and sequences.
        Yields both as individual string elements in a fast class."""
    fstq = fast()
    # FASTQ files have lines of sequence, and lines of quality
    # the quality tag indicates which type the line is, set to true if quality score
    quality = False
    for line in in_file:
        if line.startswith('@'):
            if fstq.iden is not None:
               yield(fstq)
               #reset to find new sequence
               quality = False
            fstq.seq = ""
            fstq.qual = []
            splitLine = re.search("[\s,]", line)
            fstq.iden = line[1:splitLine.start(0)]
            fstq.comment = line[splitLine.start(0):-1]
            continue

        if not quality and fstq.iden is not None:
            # The "+" line is used to distinguish when quality scores start
            if line.startswith("+"):
                quality = True
                continue
            line = line.rstrip()
            line = line.upper()
            fstq.seq += line
            continue

        if quality:
            for score in line:
                fstq.qual.append(ord(score) - number)

    yield(fstq)


def read_fasta_quality(quality_file):
    """Reads in the file containing the corresponding quality values of a FASTA file
       Stored in a fast object with iden, but only quality is yielded."""
    qual_fast = fast()
    for line in quality_file:
       if line.startswith('>'):
       #quality files also have identifiers, matching the sequences they apply to
       #this can be used to seperate scores
          if qual_fast.iden == None:
             qual_fast.iden = line[1:]
             continue
          else:
             yield(qual_fast.qual)
             qual_fast.iden = line[1:]
             qual_fast.qual = []
             continue
       else:
          for score in line.split():
             #FASTA quality files have white spaces, must use split or will store as strings
             qual_fast.qual+=[int(score)]
             continue

    if qual_fast.iden is not None:
       yield(qual_fast.qual)


