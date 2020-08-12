# python2.7
# ref: hg19
# up to 3mm and 1 base DNA/RNA bulge
# currently cant do DNA and RNA bulge together
# better local alignment requires

import sys
import math
import os
import numpy as np
import csv
import pandas as pd
import pickle
import argparse
import re
from sklearn import datasets
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE


ALL_folder_dir = '/home/yp11/Desktop/genomes/1000genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/SNPMNPINDEL'
# /hap1_HG01556_all
SNP_folder_dir = '/home/yp11/Desktop/genomes/1000genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502'
cosmid_dir = '~/Desktop/useful_tools/tagger-1.3/fetchGWI'

individual_arr = ['HG02810','HG03833', 'HG00109', 'HG01885', 'HG01556']
gRNA_arr = ['TAACGGCAGACTTCTCCAC', 'CCATTAAAGAAAATATCAT']

crista_dir = '/home/yp11/PycharmProjects/CRISTA'

mm = 3
# set indel = 1 base
# 19 bases, R66S and CFTR


def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]
# reference https://stackoverflow.com/questions/11122291/how-to-find-char-in-string-and-get-all-the-indexes


def tagscanMM(cosmid_dir, input_dir, output_dir, ind, gRNA, mm):
    # prepare gRNA sequence
    # ind = HG01556_1
    # mm only

    # ~/Desktop/useful_tools/tagger-1.3/fetchGWI -a HG01556_1_idx.txt -i HG01556_1_idx.bin -m 1 CCTGTCGTCCTTGAAGAA
    index_dir = input_dir + '/' + ind + '_idx.txt'
    target_dir = input_dir + '/' + ind + '_idx.bin'
    gRNA_scan = gRNA + 'NGG'
    qry = cosmid_dir + ' -a ' + index_dir + ' -i ' + target_dir + ' -m ' + str(mm) + ' ' +gRNA_scan
    print(qry + ' >> ' + output_dir)
    os.system(qry + ' >> ' + output_dir)
    return


def tagscanBULGE(cosmid_dir, input_dir, output_dir, ind, gRNA, mm):
    # Assume gRNA length = 19
    assert (len(gRNA) == 19)
    index_dir = input_dir + '/' + ind + '_idx.txt'
    target_dir = input_dir + '/' + ind + '_idx.bin'

    # 18 bases guides
    for i in range(0, len(gRNA)):
        j = i+1
        gRNA_18 = gRNA[:i] + gRNA[j:]
        gRNA_scan = gRNA_18 + 'NGG'
        qry = cosmid_dir + ' -a ' + index_dir + ' -i ' + target_dir + ' -m ' + str(mm) + ' ' + gRNA_scan
        os.system(qry + ' >> ' + output_dir)
        print(qry + ' >> ' + output_dir)

    # 20 bases guides
    for i in range(0, len(gRNA)):
        gRNA_20 = gRNA[:i] + 'N' + gRNA[i:]
        gRNA_scan = gRNA_20 + 'NGG'
        qry = cosmid_dir + ' -a ' + index_dir + ' -i ' + target_dir + ' -m ' + str(mm) + ' ' + gRNA_scan
        os.system(qry + ' >> ' + output_dir)
        print(qry + ' >> ' + output_dir)

        # template output
# CCTGTCGTCCTTGAANRG	CCTTTCGTCCTTGAAGAG	chr18:NC_000000.18[68885542..68885559]	-
# CCTGTCGTCCTTGAANRG	CCTGTAGTCCTTGAATAG	chr20:NC_000000.20[22955150..22955167]	-


# if PAM is different, there will be +20 for the score...
def COSMIDoff(MMpattern,NRG_):

# if there is indel, calculate with the same way
    scorelist = [0.12, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.27, 0.35, 0.5, 0.7, 0.8, 1.1, 1.3, 1.9, 2.3, 3.0,
                 4.0, 5.0, 6.0]
    COSMID_score = 0

# if indel, truncate the left side.
# del:+0.51
# ins:+0.7

    MM_index = find(MMpattern, "*")  # get a list of MM index
    # print (MM_index)

    for i in MM_index:
        COSMID_score += scorelist[i]

    if (NRG_ != 1):
        COSMID_score += 20.00

    # print -COSMID_score
    return 48.40-COSMID_score
    # larger score, higher off-target. 48.4 means perfect match.


####generate matching pattern####

def match_pattern(wt, off):
    if len(wt) == len(off) == 23:
        wt = wt[:20]
        off = off[:20]
    else:
        print('errrrrrrrrrrrrrrrrrrrrrrrrror', wt, off, len(wt), len(off))
    assert (len(wt) == len(off) == 20)

    pattern = ''
    for i in range(0, 20):
        if wt[i] != off[i]:
            pattern = pattern + '*'
        else:
            pattern = pattern + '.'
    return pattern
####generate matching pattern####


def main():
    ############################
    # #hg19
    # hg19_dir = '/home/yp11/Desktop/genomes/1000genome/chromFa_numbers'
    # ind_hap1 = 'hg19_1'
    # for gRNA in gRNA_arr:
    #     output_dir = hg19_dir + '/' + ind_hap1 + gRNA + '_MM.txt'
    #     tagscanMM(cosmid_dir, hg19_dir, output_dir, ind_hap1, gRNA, 3)
    #     output_dir_bul = hg19_dir + '/' + ind_hap1 + gRNA + '_BUL.txt'
    #     tagscanBULGE(cosmid_dir, hg19_dir, output_dir_bul, ind_hap1, gRNA, 3)
    ############################

    ############################
    # SNP only first:

    # for ind_name in individual_arr:
    #     # hap1
    #     input_dir = SNP_folder_dir + '/' + 'hap1_' + ind_name
    #     ind_hap1 = ind_name + '_1'
    #     for gRNA in gRNA_arr:
    #         output_dir = SNP_folder_dir + '/' + ind_hap1 + gRNA + '.txt'
    #         tagscanMM(cosmid_dir, input_dir, output_dir, ind_hap1, gRNA, 3)
    #         # tagscanBULGE(cosmid_dir, input_dir, output_dir, ind_hap1, gRNA, 3)
    #
    #     # hap2
    #     input_dir = SNP_folder_dir + '/' + 'hap2_' + ind_name
    #     ind_hap2 = ind_name + '_2'
    #     for gRNA in gRNA_arr:
    #         output_dir = SNP_folder_dir + '/' + ind_hap2 + gRNA + '.txt'
    #         tagscanMM(cosmid_dir, input_dir, output_dir, ind_hap2, gRNA, 3)
    #         # tagscanBULGE(cosmid_dir, input_dir, output_dir, ind_hap2, gRNA, 3)
    ############################

    ############################
    # Everything(SNP, MNP, Indel):

    # for ind_name in individual_arr:
    #     # hap1
    #     input_dir = ALL_folder_dir + '/' + 'hap1_' + ind_name +'_all'
    #     ind_hap1 = ind_name + '_1'
    #     for gRNA in gRNA_arr:
    #         output_dir = ALL_folder_dir + '/' + ind_hap1 + gRNA + '.txt'
    #         # tagscanMM(cosmid_dir, input_dir, output_dir, ind_hap1, gRNA, 3)
    #         tagscanBULGE(cosmid_dir, input_dir, output_dir, ind_hap1, gRNA, 3)
    #
    #     # hap2
    #     input_dir = ALL_folder_dir + '/' + 'hap2_' + ind_name+'_all'
    #     ind_hap2 = ind_name + '_2'
    #     for gRNA in gRNA_arr:
    #         output_dir = ALL_folder_dir + '/' + ind_hap2 + gRNA + '.txt'
    #         # tagscanMM(cosmid_dir, input_dir, output_dir, ind_hap2, gRNA, 3)
    #         tagscanBULGE(cosmid_dir, input_dir, output_dir, ind_hap2, gRNA, 3)
    ############################

    ############################
    # remove Non-NGA NAG PAMs
    # because of limitation of time, use NGG instead of NRG
    # CCTGTCGTCCTTGAANRG	CCTTTCGTCCTTGAAGAG	chr18:NC_000000.18[68885542..68885559]	-

    #for testing
    screen_dir = '/home/yp11/Desktop/genomes/1000genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/SNPMNPINDEL/' \
                 'COSMID_BULGE'
    screen_output = screen_dir + '/HG00109_1CCATTAAAGAAAATATCAT_both.txt'
    filtered_dir = screen_dir
    filtered_output = filtered_dir + '/HG00109_1CCATTAAAGAAAATATCAT_filtered.csv'
    final_output = filtered_dir + '/HG00109_1CCATTAAAGAAAATATCAT_final.csv'
    df_1 = pd.read_csv(screen_output, sep='\t', lineterminator='\n', names=['qry', 'hit', 'loc', 'strand'])

    df_gg = df_1[df_1['hit'].str[-2:] == 'GG']
    df_ga = df_1[df_1['hit'].str[-2:] == 'GA']
    df_ag = df_1[df_1['hit'].str[-2:] == 'AG']

    frames = [df_gg, df_ga, df_ag]
    df_2 = pd.concat(frames)
    df_2 = df_2.sort_values(by=['loc'])
    df_2.drop_duplicates(keep='first', inplace=True)
    # initial filtering, remove obvious duplicate
    # now still have duplicates with exact same mapping result but different alignment manner

    ############################
    # calculate COSMID score and only save the best alignment
    for row in df_2.itertuples():
        # qry and hit are aligned
        # qry	hit	loc	strand
        # CCATTAAAGAAAATATCANGG	CCATTAAACAAGATGTCAAGG	chr10:NC_000000.10[10110045..10110065]	+
        # CCATTAAAGAAAATATATNGG	CCATTAAAAACAATGTATCGG	chr10:NC_000000.10[101555033..101555053]	-
        # CNCATTAAAGAAAATATCATNGG	CACATTAAACAGAATATTATGGG	chr10:NC_000000.10[101663038..101663060]	+

        # actually we know what is wt
        # add a G to meet our design
        # since the first base shouldn't affect the efficiency, add a 'G' in ot
        wt = 'G' + 'CCATTAAAGAAAATATCAT' + 'GGG'
        # to real processing :
        # wt = 'G' + some element in that array/input + 'GGG'
        ot = 'G' + str(row.hit).upper()
        # print(wt, ot)

        if ot[-2:] == "GG" or ot[-2:] == "AG" or ot[-2:] == "GA":
            NRG_ = 1
        else:
            NRG_ = 0
        # print(NRG_)

        #cosmid calculation
        # only keep the site with highest cosmid score (lowest OT efficiency)
        if len(ot) == 23:
            # print (ot + ' len ot = 23')
            # print (wt)
            mismatchpos = match_pattern(wt, ot)
            # print (mismatchpos)
            df_2.at[row.Index, 'MM'] = mismatchpos
            df_2.at[row.Index, 'COSMID'] = COSMIDoff(mismatchpos, NRG_)


        elif len(ot) == 22:
            # need to find out which base is missing
            g_qry = 'G' + row.qry
            for i in range(0, len(wt)):
                j = i + 1
                wt_22 = wt[:i] + wt[j:]
                if (wt_22[:19] == g_qry[:19]):
                    # got the correct alignment!
                    # can use it for cosmid calculation
                    ot = ot[:i] + '-' + ot[i:]
                    # print (ot + ' len ot = 22')
                    # print (wt)
                    mismatchpos = match_pattern(wt, ot)
                    # print (mismatchpos)
                    df_2.at[row.Index, 'MM'] = mismatchpos
                    df_2.at[row.Index, 'COSMID'] = COSMIDoff(mismatchpos, NRG_) - 0.51

                    break

        elif len(ot) == 24:  # 1 base extra in ot
            wt = 'G' + str(row.qry).upper()
            wt = wt.replace('NGG', 'GGG')
            wt = wt.replace('N', '-')
            # print (ot)
            # print (wt)
            wt = wt[1:24]
            ot = ot[1:24]
            # print (ot + ' len ot = 24')
            # print (wt)
            mismatchpos = match_pattern(wt, ot)
            # print (mismatchpos)
            df_2.at[row.Index, 'MM'] = mismatchpos
            df_2.at[row.Index, 'COSMID'] = COSMIDoff(mismatchpos, NRG_) - 0.7
        df_2.at[row.Index, 'wt'] = wt
        df_2.at[row.Index, 'ot'] = ot

    # Future update required
    # using local alignment to remove MM/Bulge duplicate and keep the best alignment
    ############################

    ############################
    # then check the index
    df_2 = df_2.sort_values(by=['COSMID'], ascending=False)
    df_2.drop_duplicates(subset=['loc'], keep='first', inplace=True)

    # check point
    #####check if the table drops what we want
    ############################

    ############################
    # get CRISTA score

    # direct_output = subprocess.check_output('ls', shell=True)
    for row in df_2.itertuples():
        # get position information
        # chr18:NC_000000.18[68885542..68885559]
        # diff: the difference between the start of original and extended loc
        location_str = str(row.loc)
        chr = location_str.split(':')[0]
        start = int(location_str.split(':')[1].split('[')[1].split('.')[0])
        # end = int(location_str.split('.')[-1].replace(']', ''))

        if len(row.hit) == 23:  # 19+N+NGG
            extended_qry = 'NNN' + row.hit + 'NNN'
            diff = 3
        elif len(row.hit) == 22:  # 19+NGG
            extended_qry = 'NNNN' + row.hit + 'NNN'
            diff = 4
        elif len(row.hit) == 21:  # 19 - 1 base + NGG
            # considering they are aligned to right(PAM)
            extended_qry = 'NNNNN' + row.hit + 'NNN'
            diff = 5
        else:
            print('error: length of OT doesnt match requirements ', len(row.hit))
            break

        # search and match the target sequence
        # ~/Desktop/useful_tools/tagger-1.3/fetchGWI -a HG01556_1_idx.txt -i HG01556_1_idx.bin -m 1 CCTGTCGTCCTTGAAGAA
        # now we assume the dir directly
        # use filtered_dir
        index_dir = '/home/yp11/Desktop/genomes/1000genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/SNPMNPINDEL/hap1_HG00109_all/HG00109_1_idx.txt'
        target_dir = '/home/yp11/Desktop/genomes/1000genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/SNPMNPINDEL/hap1_HG00109_all/HG00109_1_idx.bin'
        temp_dir = '/home/yp11/Desktop/genomes/1000genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/SNPMNPINDEL/hap1_HG00109_all/temp.txt'
        qry = cosmid_dir + ' -a ' + index_dir + ' -i ' + target_dir + ' -m 0 ' + extended_qry
        # print(row.hit, extended_qry)
        # print(qry + ' > ' + temp_dir)
        os.system(qry + ' > ' + temp_dir)
        df_3 = pd.read_csv(temp_dir, sep='\t', lineterminator='\n', names=['qry', 'hit', 'loc', 'strand'])
        hit_check = 0
        for row_1 in df_3.itertuples():
            temp_location_str = str(row_1.loc)
            temp_chr = temp_location_str.split(':')[0]
            temp_origin_start = int(temp_location_str.split(':')[1].split('[')[1].split('.')[0])
            temp_start = temp_origin_start + diff
            df_3.at[row_1.Index, 'temp_start'] = temp_start

            # print(row.loc)
            # print(row_1.loc)
            # print(chr, start)
            # print(temp_chr, temp_start)

            # temp_end = int(temp_location_str.split('.')[-1].replace(']', ''))
            # print(temp_location_str, temp_chr, temp_start)
            # print(location_str,chr,start)
            # check if the seq can be matched
            if chr == temp_chr and -5 < start - temp_start < 5:  ####### A very un-recommended way to solve the missing site issue.
                df_2.at[row.Index, 'extended_seq'] = row_1.hit
                df_2.at[row.Index, 'temp_chr'] = temp_chr
                df_2.at[row.Index, 'temp_origin_start'] = temp_origin_start
                hit_check = 1
                break

        if hit_check == 0:
            print(row)
            print('no hit')
            print(df_3)
            print(qry + ' > ' + temp_dir)
            #sys.exit()

    # 2nd filtering
    # remove SNP+indel overlaps, get the best alignment(COSMID based)
    df_2 = df_2.sort_values(by=['COSMID'], ascending=False)
    df_2.drop_duplicates(subset=['temp_chr','temp_origin_start','extended_seq'], keep='first', inplace=True)
    df_2.to_csv(filtered_output, index=False, header=True)

    # then CRISTA
    for row in df_2.itertuples():
        # need to be in CRISTA repository
        # Calculate CRISTA score
        crista_qry = 'python ' + crista_dir + '/CRISTA_yp.py -s ' + 'G' + 'CCATTAAAGAAAATATCAT' + ' -d ' \
                     + row.extended_seq
        print(crista_qry)
        os.system(crista_qry + '> tmp')
        crista_score = open('tmp', 'r').read()
        df_2.at[row.Index, 'CRISTA'] = crista_score
                # save output
                # python CRISTA.py   -s  CTCAGCTGAGGTTGCTGCTG    -d  GGCCTCAGCTGATTGCTGCTGTGGAAGT


    df_2.to_csv(final_output, index=False, header=True)
    # combine repeated
    # two situations. 1: overlapping 2: one base away(only have 1 base bulge)
    # have one column calculate all (updated COSMID score)
    # and another column calculate CRISTA score
    # CRISTA requires python 3











    #CRISTA
    # Add a G before gRNA sequence
    # straight forward for SNP only query

if __name__ == '__main__':
    main()
