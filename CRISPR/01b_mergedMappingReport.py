import os
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-ip','--inPath', help='Path to the base run processed folder',required=True)

args = parser.parse_args()
inpath_ = args.inPath
inpath = inpath_ + '/QC'

# define output path
#inpath='/scratch/julen/pruebas/CRISPR/QC'
mergeStats = inpath_ + f"/RSession/{inpath_.split('/')[-1]}_sequencingStats.tsv"

# define variables with the location of the interest headers
positions = [ "N_reads\tN_guide\tPercentaje\tN_duplex_guides\n",
            "ALIGNMENT\n",
            "ALIGNMENT of previously unmapped\n",
            "Number of reads that map more than once, and their top target\n",
            "SAMTOOLS FLAGSTAT\n",
            ]
posIndex = [0] * len(positions)

# get list of valid files
symmaryF = [o for o in os.listdir(inpath) if o.startswith('summary')]

# get and format the interest text
allText = ""
for fi in symmaryF:
    with open(f'{inpath}/{fi}', 'r') as f:
        data=f.readlines()
        
    # get position of each location
    for npo, po in enumerate(positions):
        for nda, da in enumerate(data):
            if po == da:
                posIndex[npo] = nda
                
    # get text to write
    header = f'SampleID\t'
    content = f'{fi[8:-4]}\t'

    header += data[posIndex[0]][:-1] + '\t'
    content += '\t'.join(data[posIndex[0] + 1][:-1].split()) + '\t'

    header += 'UnAligned\tAlignedOnce\tAligned>1\t'
    content += '\t'.join([f'{int(d.split()[0]):,}' for d in data[posIndex[1] + 4: posIndex[1] + 7]]) + '\t'

    header += 'UnAlignedNonV\tAlignedOnceNonV\tAligned>1NonV\t'
    if posIndex[2] != 0:
        if data[posIndex[2] + 2] ==  '0 reads\n':
            content += '0\t0\t0\t'
        else:
            content += '\t'.join([f'{int(d.split()[0]):,}' for d in data[posIndex[2] + 4: posIndex[2] + 7]]) + '\t'
    else:
        content += 'NotChecked\tNotChecked\tNotChecked\t'

    header += 'Aligned>1Info'
    tmptext = ''
    for dd in [d.split() for d in data[posIndex[3] + 2: posIndex[4] - 1]]:
        tmptext += f'{dd[1]}({dd[0]});'
    tmptext = tmptext[:-1]
    content += tmptext
    header += '\n'

    allText += content + '\n'
    
# write all
with open(mergeStats, 'w') as fout:
    fout.write(header)
    fout.write(allText)
