#!/usr/bin/env python

import pandas as pd
import argparse
import os
import math

# we don't want to get warnings. It works.
pd.options.mode.chained_assignment = None 

# Function to calculate -log10(x) without handling x as a float:
def get_log(pv):
    if 'E' in pv.upper():
        mantissa, exponent = pv.upper().split('E')
        if float(mantissa) == 0:
            return(0)
        else:
            return(math.log10(float(mantissa))+float(exponent))
    else:
        return(math.log10(float(pv)))

# Function to find and annotate top association:
def find_top_association(chrom):    
    global input_df
    global window
    global threshold

    # The dataframe is filtered by chromosome:
    test_df = input_df.loc[input_df.chromosome == chrom]

    # Calculate the -log(pv) and update column type of the position:
    test_df['mLogPv'] = test_df.pvalue.apply(get_log)
    test_df.loc[:,'bp_location'] = test_df['bp_location'].astype(float)

    # sorting rows by p-value:
    for index in test_df.sort_values(['mLogPv']).index:
        
        # Check if the index has 'false' or 'review' flag:
        if test_df.loc[index]['isTopAssociation'] == 'false' or \
            test_df.loc[index]['isTopAssociation'] == 'REQUIRES REVIEW':
            continue

        # Excluding sub significant variants.
        elif test_df.loc[index]['mLogPv'] > threshold:
            test_df.loc[index, 'isTopAssociation'] = 'false'
            continue

        # We have found a top snp!  
        else:
            pos = test_df.loc[index]['bp_location']
            pval = test_df.loc[index]['mLogPv']

            # Setting false flag for ALL variants within the window.
            test_df.loc[test_df[abs( test_df.bp_location - pos ) <= window].index,'isTopAssociation'] = 'false'
            
            # If there are multiple variants with the same p-value within the window they require review:
            if len(test_df.loc[(test_df.mLogPv == pval) & (abs( test_df.bp_location - pos ) <= window)]) > 1:
                test_df.loc[ (test_df.mLogPv == pval) & (abs( test_df.bp_location - pos ) <= window), 'isTopAssociation'] = 'REQUIRES REVIEW'
            # Otherwise it's a true peak:
            else:
                test_df.loc[ index, 'isTopAssociation'] = 'true'

    # Modifying the original dataframe: 
    input_df.loc[ test_df.index, 'isTopAssociation'] = test_df.isTopAssociation

if __name__ == '__main__':

    # Parsing commandline arguments
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='This script finds the most significant association within a defined range (100kbp by default).')

    parser.add_argument('-f', '--input', help='Input file name with table of associations.')
    parser.add_argument('-o', '--output', help='Output file name.')
    parser.add_argument('-w', '--window', default=100000, help='Window size.', type = int)
    parser.add_argument('-t', '--threshold', default=1e-5, help='p-value threshold.', type = float)
    parser.add_argument('-p', '--prune', default=False, help='Prune out sub significant associations from the output.', action='store_true')

    args = parser.parse_args()

    # inputFile = args.input
    outputFile = args.output
    window = args.window
    threshold = math.log10(args.threshold)
    prune = args.prune

    if not outputFile:
        raise(Exception("[Error] A output file needs to be specified! Exiting."))

    # If inputfile is not specified or not submitted exiting:
    if args.input and os.path.isfile(args.input):
        inputFile = args.input
    else:
        raise(Exception("[Error] A valid input file is required! Exiting."))

    # Reading input file into pandas dataframe:
    try:
        input_df = pd.read_csv(inputFile, sep="\t", dtype = str)
    except:
        raise(Exception("[Error] Input file could not be read. Please provide a tab separated table. Exiting."))
    
    # Checking header (setting lowercase):
    input_df.columns = map(str.lower, input_df.columns)
    if not pd.Series(['rs_id', 'pvalue', 'chromosome', 'bp_location']).isin(input_df.columns).all():
        raise Exception('[Error] Not all required columns were found in the file header. Required columns: "rs_id", "pvalue", "chromosome", "bp_location"')

    if prune:
        input_df = input_df[ input_df.mLogPv < threshold ]

    input_df['isTopAssociation'] = ''

    # finding peaks for each chromosome
    for chromosme in input_df.chromosome.unique():
        find_top_association(chromosme)

    # Saving the modified table into a tab separated file:
    # input_df.drop(["mLogPv"], inplace=True)
    input_df.to_csv(outputFile, sep="\t", index= False, na_rep = 'NA')
