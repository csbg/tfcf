import pyfastx
import re
import io
import gzip
import argparse

def writefastq(name1, seq1, uselessl, qual1, fout1):
    fout1.write(name1 + '\n')
    fout1.write(seq1 + '\n')
    fout1.write(uselessl)
    fout1.write(qual1 + '\n')
    
def writefastqGz(name1, seq1, uselessl, qual1, fout1):
    name1 = name1 + '\n'
    seq1 = seq1 + '\n'
    qual1 = qual1 + '\n'
    fout1.write(name1.encode())
    fout1.write(seq1.encode())
    fout1.write(uselessl.encode())
    fout1.write(qual1.encode())
    

def filterFastqPattern(fi1, flanquingSeq='CACCG(.{20})GTTTTAGAGC', 
                       writeGz=False):
    '''
    :param fi1: path to the fastq file to be processed. It must
        be gziped and finish in .gz
    :param 'CACCG(.{20})GTTTTAGAGC' flanquingSeq: regexp pattern
        to use at the time to search for the sequence flanquing 
        our 20bp sequence of interest
    :param False writeGz: set to True if you want to save output
        extrated fastq file in gzipped format
    '''
    # first define the function to use if we want to write 
    # output copressed or readable
    if writeGz == True:
        writeF = writefastqGz
        outname = fi1.split('.gz')[0] + '_extracted.fastq.gz'
        fout1 = io.BufferedWriter( gzip.open(outname,'wb'))
    else:
        writeF = writefastq
        outname = fi1.split('.fastq')[0] + '_extracted.fastq'
        fout1 = open(outname,'w')

    # Some variables to extrac the flanked match from the whole match
    # get the number of base pairs flanquing our sequence of interest
    regxLeft = len(flanquingSeq.split('(.{')[0])
    # get the length of our sequence of interest
    regxIn = int(flanquingSeq.split('(.{')[-1].split('})')[0])

    # useless line of fastq
    uselessl = '+\n'

    fileInfo = {'nRead':0, 'nGuides':0, 'nDuplex':0}
    fq1 = pyfastx.Fastx(fi1)
    for data1 in fq1:
        fileInfo['nRead'] += 1
        name1,seq1,qual1,comment1 = data1
        readId = f'@{name1} {comment1}'
        # look for regexp pattern in the sequence
        patternLoc=list(re.finditer(flanquingSeq, seq1))
        # one match is what we look for
        if len(patternLoc) == 1:
            fileInfo['nGuides'] += 1
            indice = [m.start(0) for m in patternLoc][0]

            trimSeq1 = seq1[indice+regxLeft:indice+regxLeft+regxIn] 
            trimQual = qual1[indice+regxLeft:indice+regxLeft+regxIn] 

            writeF(readId, trimSeq1, uselessl, trimQual, fout1)

        # more than one match is weird
        elif len(patternLoc) > 1:
            fileInfo['nDuplex'] += 1
            print(f'read {readId} was found to match the pattern more than once')
            print(seq1)
            print('Skipped')

    print('N_reads\tN_guide\tPercentaje\tN_duplex_guides')
    percent = round(fileInfo["nGuides"]/fileInfo["nRead"] * 100)
    print(f'{fileInfo["nRead"]:,}', f'{fileInfo["nGuides"]:,}', 
          f'{percent}%',  f'{fileInfo["nDuplex"]:,}')

############

parser = argparse.ArgumentParser(description='')
parser.add_argument('-fq','--fastq_path', help='Path to GZIP-ed FASTQ file',required=True)
parser.add_argument('-ptrn','--pattern', help='Regexp patern to retrieve the sequences of interest',required=False)
parser.add_argument('-outgz','--outgz', help='Set to "True" to write output fastq file already GZIP-ed',required=False)

## Load script input variables
args = parser.parse_args()
fi1 = args.fastq_path
flanquingSeq = args.pattern
outgz = args.outgz
if outgz == 'True':
    outgz = True
else:
    outgz = False
## Run
# if we didnt provide a pattern, use default
if flanquingSeq == None:
    filterFastqPattern(fi1, 
                       writeGz=outgz)
else:
    filterFastqPattern(fi1, flanquingSeq=flanquingSeq,
                       writeGz=outgz)

