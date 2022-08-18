import os
import pandas as pd
from collections import defaultdict
import itertools
from matplotlib import pyplot as plt
import seaborn as sns
from bioinfokit import analys, visuz
import numpy as np
import matplotlib.backends.backend_pdf
from matplotlib.colors import LinearSegmentedColormap


import sys
import pickle

# Function to filter the cases when we have zeros
def filterZeros(df, cells1_cpm, cells2_cpm):
    # get positions with zero in one and bellow average in other
    pos1 = df[cells1_cpm] == 0
    nonZero = df[cells2_cpm] != 0
    pos1 = np.array(pos1) & (np.array(df[cells2_cpm] < np.median(df.loc[nonZero, cells2_cpm])))

    # same here
    pos12 = df[cells2_cpm] == 0
    nonZero = df[cells1_cpm] != 0
    pos12 = np.array(pos12) & (np.array(df[cells1_cpm] < np.median(df.loc[nonZero, cells1_cpm])))
    
    # merge and invert
    pos1_invert = (np.array(pos1) | np.array(pos12)) == False
    
    # now keep only the ones that are in called peaks for this comparison
    cells1_bool = cells1_cpm[:-3] + 'bool'
    cells2_bool = cells2_cpm[:-3] + 'bool'
    posb = df[cells1_bool] == True
    posb = posb | (df[cells2_bool] == True)
    
    # return only positon of comparisons with at least one called peak and the zero criteria
    return pos1_invert & posb


# Function to get fold change (special care about what to do with zeros)
def getFC(df, pos, cells1_cpm, cells2_cpm):
    df2 = df.loc[pos,:].copy()
    
    # update zero values to minimum in cell 1
    minVal = min(df2.loc[df2[cells1_cpm] != 0, cells1_cpm])
    df2.loc[df2[cells1_cpm] == 0, cells1_cpm] = minVal
    
    # update zero values to minimum in cell 2
    minVal = min(df2.loc[df2[cells2_cpm] != 0, cells2_cpm])
    df2.loc[df2[cells2_cpm] == 0, cells2_cpm] = minVal
    

    df2.loc[:,'log2FC'] = np.log2(df2.loc[:, cells1_cpm]  / df2.loc[:, cells2_cpm] )
    return df2.loc[:,:]


def filterAndFc(df, cells1_cpm, cells2_cpm, minCPM, minRead):
    # filter the cases with one zero
    pos1 = filterZeros(df, cells1_cpm, cells2_cpm)

    # now we discard the ones bellow CPM thresshold in both cells
    pos2 = df[cells1_cpm] >= minCPM
    pos2 = np.array(pos2) | (np.array(df[cells2_cpm] >= minCPM))
    # merge
    pos = pos1 & pos2
    
    # we filter also the peaks with less than 10 reads in both cells
    cells1_rcount = cells1_cpm[:-3] + 'rcount'
    cells2_rcount = cells2_cpm[:-3] + 'rcount'

    posR = df[cells1_rcount] >= minRead
    posR = np.array(posR) | (np.array(df[cells2_rcount] >= minRead))

    pos = pos & posR
    
    # get fold change and work with zero cases as defined in function
    return getFC(df, pos, cells1_cpm, cells2_cpm)


def filterByControlAndCPM(df, cells1_cpm, cells2_cpm, foldChanges, controlFC=1,
    CPMratioFilter=None):
    ### Function to filter by control and by asigned CPM ratio value
    
    if controlFC != None:
        ## Filter out positions with too small fold change between
        # cell and its control
        ## Get controls
        controls_ = [d for d in df.columns if d.endswith('cpm') 
                                            and ('IgG' in d or 'input' in d)]
        # get cell controls
        control1 = [c for c in controls_ if c.startswith(cells1_cpm.split('_')[0])]
        if len(control1) == 0:
            control1 = [c for c in controls_ if c.startswith('allCell')][0]
        else:
            control1 = control1[0]

        control2 = [c for c in controls_ if c.startswith(cells2_cpm.split('_')[0])]
        if len(control2) == 0:
            control2 = [c for c in controls_ if c.startswith('allCell')][0]
        else:
            control2 = control2[0]

        controlFold1 = np.log2(df.loc[:, cells1_cpm]  / 
                            df.loc[:, control1] )

        posC1 = (controlFold1 <= -controlFC) & (df.loc[:,'log2FC'] >= foldChanges[0])

        controlFold2 = np.log2(df.loc[:, cells2_cpm]  / 
                                            df.loc[:, control2] )

        posC2 = (controlFold2 <= -controlFC) & (df.loc[:,'log2FC'] <= -foldChanges[1])

        # get reverse
        pos = (posC1 | posC2) == False   
    
    else:
        pos = df.loc[:,'.bool'] == False
    
    
    ## Now filter by ratio if stated
    if CPMratioFilter != None:
        compareId = (cells1_cpm[:-4], cells2_cpm[:-4])
        if compareId in CPMratioFilter:
            keepPos = ((np.log2(df[cells1_cpm]) + 
                        np.log2(df[cells2_cpm])) / 2) >= CPMratioFilter[compareId]
            pos = pos & keepPos

    df3 = df.loc[pos,:].copy()
    return df3

# Path with the ChIP tables
tableAll = "/PATH/TO/ChIP/TABLES"
# Path to output table with FC info
outFC = '/PATH/FOR/OUTPUT/TABLES'
# Path to DESeq differential analysis results
deseqP = '/PATH/TO/DESeq/FILES'



cellOrder = ['DM', 'Mye', 'LSK', 'GMP', 'Mono', 'MEP', 'Ery', 
             'Bcell']

allowedCompare = [['DM', 'Mye', 'Mono', 'MEP', 'GMP'], 
                  ['LSK', 'MEP', 'Ery', 'GMP', 
                  'Mono', 'Bcell']]

FCthressholds = {'DM':1, 'Mye':1, 'LSK':1,
                'MEP':1, 'Ery':1, 'Mono':1, 'GMP':1, 'Bcell':1}
# define thresshold to store fold hcange
# minimum CPM value needed in at lesat one of the cells
minCPM = 2

minRead = 10

noiseFilter = False

# this is to filter by the ratio of CPM (as you see in x axis of MA plots)
CPMratioFilter = None

files = [f for f in os.listdir(tableAll) if f.endswith('boolean.annotatePeaks.extended.txt')]

# Make FC tables
for peakType in ['narrowPeak']:
    print('#'*30, peakType, '#'*30)
    for fi in sorted(files):
        if peakType in fi:
            chip = fi.split('_')[0]
            print('-'*30, chip, '-'*30)
            with open(f'{tableAll}/{fi}', 'r') as f:
                header = f.readline().strip().split('\t')

            df = pd.read_csv(f'{tableAll}/{fi}', sep='\t')
            df['copy_index'] = df.index
            df_out = df.copy()
            reAsignVals = {}
            
            # get cell ids 
            chipCell = []
            allchip = []
            for li in list(df.columns):
                if li.endswith('.bool') and li != '.bool':
                    chipCell += [li.split('.')[0]]

            # get cels in desired order
            currentOrder = []
            for cell in cellOrder:
                for cc in chipCell:
                    if cc.startswith(f"{cell}_"):
                        currentOrder += [cc]
                                
            
            for cells in itertools.combinations(currentOrder, 2):
                
                # only allow certain comparisons
                allow = False
                for co in allowedCompare:
                    if (cells[0].split('_')[0] in co) and (cells[1].split('_')[0] in co):
                        allow = True
                        
                if allow == True:
                    cells1_cpm = f'{cells[0]}.cpm'
                    cells2_cpm = f'{cells[1]}.cpm'
                    
                    foldChanges = (FCthressholds[cells[0].split('_')[0]],
                                   FCthressholds[cells[1].split('_')[0]])

                    # filter values and get FC
                    df2 = filterAndFc(df, cells1_cpm, cells2_cpm,
                                                         minCPM, minRead)
                    # filter by input or IgG
                    df3 = filterByControlAndCPM(df2, cells1_cpm, 
                                                             cells2_cpm, foldChanges, 
                                                             controlFC=1,
                                                            CPMratioFilter=CPMratioFilter)
                    
                    # check if we have these cells analysed with replicates
                    cellsShort = [c.split('_')[0].split('-')[0] for c in cells]
                    deseqFile1 = f"{chip}_{peakType}_{cellsShort[0]}-vs-{cellsShort[1]}.deseq2.results.txt"
                    deseqFile2 = f"{chip}_{peakType}_{cellsShort[1]}-vs-{cellsShort[0]}.deseq2.results.txt"
                    deseqFile = None
                    if os.path.exists(f"{deseqP}/{deseqFile1}"):
                        deseqFile =  deseqFile1
                    if os.path.exists(f"{deseqP}/{deseqFile2}"):
                        deseqFile =  deseqFile2
                    if deseqFile != None:
                        # open table and locate adjusted p values and points to mark
                        df_re = pd.read_csv(f'{deseqP}/{deseqFile}', sep='\t')
                        pcol = [d for d in df_re.columns if d.endswith('padj')][0]
                        fcol = [f for f in df_re.columns if f.endswith('log2FoldChange')][0]

                        # rename and merge
                        df_re = df_re.rename(columns={'Geneid':'interval_id', 'Chr':'chr',
                                                     'Start':'start', 'End':'end'})
                        df_new = pd.merge(df3, df_re.loc[:,['chr', 'start', 'end', 'interval_id', pcol, fcol]], 
                                                 on=['interval_id', 'chr', 'start', 'end'])
                        if len(df_new.index) == 0:
                            print('No match between coordinates of DESeq and peak files (diff versions)')
                            ERROR
                            
                        # Correct FC values if order was different than what we look in here
                        if deseqFile2 == deseqFile:
                            df_new[fcol] = df_new[fcol] * -1
                        
                        # Since we have replicates, reduce log2(FC) thresshold
                        foldChanges = [0.75, 0.75]
                        df_new.index = df_new['copy_index']
                        
                        pos = (df_new[pcol] <= 0.01)
                        posu = df_new.loc[:, fcol] >= foldChanges[0]
                        posd = df_new.loc[:, fcol] <= -foldChanges[1]
                        posu = pos & posu
                        posd = pos & posd
                        
                        # rename some columns
                        df_new.drop(columns=['log2FC'], inplace=True)
                        df_new.rename(columns={fcol:'log2FC'}, inplace=True)
                        
                        
                        
                    else:
                        df_new = df3.copy()
                        df_new.index = df_new['copy_index']
                        posu = df_new.loc[:, 'log2FC'] >= foldChanges[0]
                        posd = df_new.loc[:, 'log2FC'] <= -foldChanges[1]
                        
                        
                        
                    
                    df_new['FCpass'] = 0
                    df_new.loc[posu, 'FCpass'] = 1
                    df_new.loc[posd, 'FCpass'] = -1

                    
                    # Store Fold change information in main table
                    if deseqFile != None:
                        df_out = df_out.merge(
                                    df_new[['interval_id', pcol, 'log2FC', 'FCpass']].rename(columns={
                                                            pcol: '-VS-'.join(cells) + '.DESeqPval',
                                                            "log2FC": '-VS-'.join(cells) + '.log2FC',
                                                             "FCpass": '-VS-'.join(cells) + '.FCpass'},            ),
                                    on="interval_id", how='left')
                    else:
                        df_out = df_out.merge(
                                    df_new[['interval_id', 'log2FC', 'FCpass']].rename(columns={
                                                            "log2FC": '-VS-'.join(cells) + '.log2FC',
                                                             "FCpass": '-VS-'.join(cells) + '.FCpass'},            ),
                                    on="interval_id", how='left')
                    df_out.index = df_out['copy_index']
                    
                    # The main issue now are the CPM values that were zero and we reasigned
                    # to a minimum value, here we will check them, and if at any comparison 
                    # the same CPM is computed differently will store the mean
                    ## store all values to be reasigned (because we had an
                    # above average CPM value in one cell and zero in the other)
                    # store postions that are zero valued in cpm1
                    pos = (df_out[cells1_cpm] == 0) & (df_out['interval_id'].isin(df_new['interval_id']))
                    posIndex = df_out[pos].index
                    if cells1_cpm not in reAsignVals:
                        reAsignVals[cells1_cpm] = defaultdict(list)
                    for index in posIndex:
                        if ((float(df_out.loc[index, cells1_cpm]) == 0) and
                            float(df_new.loc[index, cells1_cpm]) != 0):
                            reAsignVals[cells1_cpm][index] += [df_new.loc[index, cells1_cpm]]

                    # same for CPM 2
                    pos = (df_out[cells2_cpm] == 0) & (df_out['interval_id'].isin(df_new['interval_id']))
                    posIndex = df_out[pos].index
                    if cells2_cpm not in reAsignVals:
                        reAsignVals[cells2_cpm] = defaultdict(list)
                    for index in posIndex:
                        if ((float(df_out.loc[index, cells2_cpm]) == 0) and
                            float(df_new.loc[index, cells2_cpm]) != 0):
                            reAsignVals[cells2_cpm][index] += [df_new.loc[index, cells2_cpm]]

            
            
            # asign values and throw warning if there was more than one
            print('Listing positions with more than non-equal reasigned CPM')
            print("line 1 --> CPMcol\tindex")
            print("line 2 --> all values")
            for cpm in reAsignVals:
                for index in reAsignVals[cpm]:
                    if len(reAsignVals[cpm][index]) > 1:
                        if len(set(reAsignVals[cpm][index])) > 1:
                            print(cpm, index)
                            print(str(reAsignVals[cpm][index])[1:-1])

                        df_out.loc[index, cpm] = np.mean(reAsignVals[cpm][index])

                    else:
                        df_out.loc[index, cpm] = reAsignVals[cpm][index][0]
                        
            df_out.drop(columns=['copy_index'], inplace=True)
            df_out.to_csv(f"{outFC}/{fi[:-4]}.fc.txt", sep='\t', index=False, na_rep='NA')
            

# Check everything is ok
for peakType in ['narrowPeak']:
    print('#'*30, peakType, '#'*30)
    for fi in sorted(files):
        if peakType in fi:
            chip = fi.split('_')[0]
            print('-'*30, chip, '-'*30)
            with open(f'{tableAll}/{fi}', 'r') as f:
                header = f.readline().strip().split('\t')

            df = pd.read_csv(f'{tableAll}/{fi}', sep='\t')
            df['copy_index'] = df.index
            df_out = pd.read_csv(f"{outFC}/{fi[:-4]}.fc.txt", sep='\t')
            
            # get cell ids 
            chipCell = []
            allchip = []
            for li in list(df.columns):
                if li.endswith('.bool') and li != '.bool':
                    chipCell += [li.split('.')[0]]

            # get cels in desired order
            currentOrder = []
            for cell in cellOrder:
                for cc in chipCell:
                    if cc.startswith(f"{cell}_"):
                        currentOrder += [cc]
                                
            
            for cells in itertools.combinations(currentOrder, 2):
                
                # only allow certain comparisons
                allow = False
                for co in allowedCompare:
                    if (cells[0].split('_')[0] in co) and (cells[1].split('_')[0] in co):
                        allow = True
                        
                if allow == True:
                    cells1_cpm = f'{cells[0]}.cpm'
                    cells2_cpm = f'{cells[1]}.cpm'
                    
                    foldChanges = (FCthressholds[cells[0].split('_')[0]],
                                   FCthressholds[cells[1].split('_')[0]])

                    # filter values and get FC
                    df2 = filterAndFc(df, cells1_cpm, cells2_cpm,
                                                         minCPM, minRead)
                    # filter by input or IgG
                    df3 = filterByControlAndCPM(df2, cells1_cpm, 
                                                             cells2_cpm, foldChanges, 
                                                             controlFC=1,
                                                            CPMratioFilter=CPMratioFilter)
                    
                    # check if we have these cells analysed with replicates
                    cellsShort = [c.split('_')[0].split('-')[0] for c in cells]
                    deseqFile1 = f"{chip}_{peakType}_{cellsShort[0]}-vs-{cellsShort[1]}.deseq2.results.txt"
                    deseqFile2 = f"{chip}_{peakType}_{cellsShort[0]}-vs-{cellsShort[1]}.deseq2.results.txt"
                    deseqFile = None
                    if os.path.exists(f"{deseqP}/{deseqFile1}"):
                        deseqFile =  deseqFile1
                    if os.path.exists(f"{deseqP}/{deseqFile2}"):
                        deseqFile =  deseqFile2
                    if deseqFile != None:
                        # open table and locate adjusted p values and points to mark
                        df_re = pd.read_csv(f'{deseqP}/{deseqFile}', sep='\t')
                        pcol = [d for d in df_re.columns if d.endswith('padj')][0]
                        fcol = [f for f in df_re.columns if f.endswith('log2FoldChange')][0]

                        # rename and merge
                        df_re = df_re.rename(columns={'Geneid':'interval_id', 'Chr':'chr',
                                                     'Start':'start', 'End':'end'})
                        df_new = pd.merge(df3, df_re.loc[:,['chr', 'start', 'end', 'interval_id', pcol, fcol]], 
                                                 on=['interval_id', 'chr', 'start', 'end'])
                        
                        if len(df_new.index) == 0:
                            print('No match between coordinates of DESeq and peak files (diff versions)')
                            ERROR
                            
                        # Correct FC values if order was different than what we look in here
                        if deseqFile2 == deseqFile:
                            df_new[fcol] = df_new[fcol] * -1
                            
                        # Since we have replicates, reduce log2(FC) thresshold
                        foldChanges = [0.75, 0.75]
                        df_new.index = df_new['copy_index']
                        
                        pos = (df_new[pcol] <= 0.01)
                        posu = df_new.loc[:, fcol] >= foldChanges[0]
                        posd = df_new.loc[:, fcol] <= -foldChanges[1]
                        posu = pos & posu
                        posd = pos & posd
                        
                        # rename some columns
                        df_new.drop(columns=['log2FC'], inplace=True)
                        df_new.rename(columns={fcol:'log2FC'}, inplace=True)
                        
                        
                        
                    else:
                        df_new = df3.copy()
                        df_new.index = df_new['copy_index']
                        posu = df_new.loc[:, 'log2FC'] >= foldChanges[0]
                        posd = df_new.loc[:, 'log2FC'] <= -foldChanges[1]
                        
                        
                        
                    
                    df_new['FCpass'] = 0
                    df_new.loc[posu, 'FCpass'] = 1
                    df_new.loc[posd, 'FCpass'] = -1

                    checkIndex = df_new.index
                    pos1 = (round(df_out.loc[checkIndex, cells1_cpm], 6) != 
                           round(df_new.loc[checkIndex, cells1_cpm], 6))
                    
                    pos2 = (round(df_out.loc[checkIndex, cells2_cpm], 6) != 
                           round(df_new.loc[checkIndex, cells2_cpm], 6))

                    if sum(pos1) != 0 or sum(pos2) != 0:
                        lll
                        
                    