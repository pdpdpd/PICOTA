
# ---------------------------------import---------------------------------------

import sys, csv, operator
import argparse

# ---------------------------------import---------------------------------------


# ---------------------------------functions------------------------------------
def read_by_line(url): # read by line
    lines = url.splitlines()
    return lines


def my_sq_range(start, end, step):
    while start < end:
        yield start
        start += step


# find the flanking sequences for primer design
def flanking_seq(chr, start_pos):
    try:
        start_ = int(snp_) - 30 + 250
        end_ = int(snp_) + 30 + 250  # now the position is not at the begining of where the SNP starts
    except:
        print
        "Wrong pos"
    else:
        url_ = "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=" + chr_ + "%3A" + str(start_) + "," + str(end_)

        # print url_
        req_ = urllib2.Request(url_)
        resp_ = urllib2.urlopen(req_)
        respHtml_ = resp_.read()
        DNA_seq = respHtml_[respHtml_.find('</DNA>') - 63:respHtml_.find('</DNA>')]
        DNA_ = DNA_seq.replace('\n', '')
        # print "DNA_"+DNA_
    return DNA_

# ---------------------------------functions------------------------------------


# get parsed information
parser = argparse.ArgumentParser(description='Automatic primer design')
parser.add_argument('--ref', default='hg38', type=str,
                    help='Reference genome. Default: hg38.')
parser.add_argument('--min', default=280, type=int,
                    help='Minimum amplicon size.')
parser.add_argument('--max', default=320, type=int,
                    help='Maximum amplicon size.')
args = parser.parse_args()