#!/usr/bin/env python

import os
import pandas as pd
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Aggregate columns from merged MACS narrow or broad annotation files and combine with featureCounts output'
Epilog = """Example usage: python 02_NR_gatherAllInTable.py <ANNOTATE_CONSENSUS_PATH> <FEATURECOUNTS_PATH> <OUTPATH>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('ANNOTATE_CONSENSUS_PATH', help="Path where the annotated consesus peak files are located (ending with boolean.annotatePeaks.txt)")
argParser.add_argument('FEATURECOUNTS_PATH', help="Path with the featureCounts raw output (ends with featureCounts.txt) and their CPM files (ends with featureCounts.CPM.txt)")
argParser.add_argument('OUTPATH', help="Path to store final merged table")

args = argParser.parse_args()


############################################
############################################
## HELPER FUNCTIONS
############################################
############################################

def createMergedTable(inpath, CPMpath, tableOut):
    '''
    Function to merge peak annotation, coverage, CPM etc. data
    :param inpath: Path with the annotated consensus peak files
        (ending with boolean.annotatePeaks.txt)
    :param CPMpath: Path with the featureCounts raw output (ends
        with featureCounts.txt) and their CPM files (ends with
        featureCounts.CPM.txt)
    :param tableOut: Path to store final merged table
    
    '''
    
    if not os.path.exists(tableOut):
        os.makedirs(tableOut)

    # here we will get all the parameters to call the motif analysis
    for fi in os.listdir(inpath):
        if not (fi.startswith('IgG_') or fi.startswith('input_')):
            if fi.endswith('boolean.annotatePeaks.txt'):
                with open(f'{inpath}/{fi}', 'r') as f:
                    header = f.readline().strip().split('\t')
                chip = fi.split('_')[0]
                print(chip)
                print(fi)

                ## File 1
                # load annotation file
                df = pd.read_csv(f'{inpath}/{fi}', sep='\t')

                # get info about number of cells
                chipCell = []
                for li in list(df.columns):
                    if li.endswith('.bool') and li != '.bool':
                        chipCell += [li.split('.')[0]]
                    
                ## File 2
                # get equivalent files with CPM and counts
                featureFiles = [fi2 for fi2 in os.listdir(CPMpath) if fi2.startswith(fi.split('.')[0])]

                # load CPM file
                fileCPM = [fi2 for fi2 in featureFiles if fi2.endswith('featureCounts.CPM.txt')][0]
                print(fileCPM)
                df2 = pd.read_csv(f'{CPMpath}/{fileCPM}', sep='\t')
                # update column names
                cellNchange = {}
                for cell in chipCell:
                    cellNchange[cell] = f'{cell}.cpm'
                # if I added controls, update also their names
                for cell in df2.columns: 
                    if ('_IgG' in cell) or ('input_' in cell):
                        cellNchange[cell] = f'{cell}.cpm'
                df2 = df2.rename(columns=cellNchange)

                ## File 3
                # Load raw reads file
                fileCounts = [fi2 for fi2 in featureFiles if fi2.endswith('featureCounts.txt')][0]
                print(fileCounts)
                df3 = pd.read_csv(f'{CPMpath}/{fileCounts}', sep='\t', skiprows=1)
                df3 = df3.rename(columns={"Geneid": "interval_id"})

                # keep only the column with the range ID and the counts and change names
                keepCol = [df3.columns[0]] + list(df3.columns[6:])
                df3 = df3.loc[:,keepCol]
                cellNchange = dict([(k, f"{'_'.join(k.split('/')[-1].split('.')[0].split('_')[:3])}.rcount") for k in keepCol[1:]])
                df3 = df3.rename(columns=cellNchange)


                ## Merge all
                df_new = pd.merge(df, df3, on="interval_id")
                df_new = pd.merge(df_new, df2, on="interval_id")

                # store table
                fiout = fi[:-3] + 'extended.txt'
                df_new.to_csv(path_or_buf=f'{tableOut}/{fiout}', sep='\t', na_rep='', 
                             float_format=None, header=True, index=False, mode='w')


############################################
############################################
## MAIN FUNCTION
############################################
############################################

createMergedTable(inpath=args.ANNOTATE_CONSENSUS_PATH, 
                    CPMpath=args.FEATURECOUNTS_PATH,
                    tableOut=args.OUTPATH)

