over_string = "ATGGGCTATGGGCT"
over_string_rc = "AGCCCATAGCCCAT"

import random
import gzip
import sys

infile = sys.argv[1]

def gc_content(seq1, seq2):
    cur_seq = seq1[24:].strip() + seq2.strip()
    gc_count = cur_seq.count("G") + cur_seq.count("C")
    return gc_count / 126.0

def overrepresent(seq1, seq2):
    cur_seq = seq1[24:].strip() + seq2.strip()
    if over_string in cur_seq:
        return True
    if over_string_rc in cur_seq:
        return True

reader = gzip.open("./second_trim_data/%s" % infile, 'rb')
outfile = gzip.open("./second_trim_data_gc_filter_strict/%s" % infile, 'wb')


while True:
    id1 = reader.readline()
    if not id1:
        break
    seq1 = reader.readline()
    reader.readline()
    qual1 = reader.readline()
    id2 = reader.readline()
    seq2 = reader.readline()
    reader.readline()
    qual2 = reader.readline()
    if overrepresent(seq1, seq2):
        if random.random() > 0.01:
            continue
    outfile.write("%s%s+\n%s%s%s+\n%s" % (id1, seq1, qual1, id2, seq2, qual2))
outfile.close()
