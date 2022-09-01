
# This script will convert HOMER motif files to a format and naming 
#suitable for TOBIAS


# Path to R 
R=/PATH/TO/R/bin/R
# base output dir
outdir=/PATH/TO/OUTDIR
# Path to HOMER motifs to be check
motifsCheck=/PATH/TO/HOMER/data/knownTFs/vertebrates/known.motifs


## First of all we convert homer motifs to a format readable by TOBIAS
mkdir -p ${outdir}/motifs
# Convert homer motifs to JASPAR  matrix
Rcommand="
library(universalmotif)
library('stringr')
motifsCheck='${motifsCheck}'
fileout='${outdir}/motifs/known_jaspar'
# First we get all homer Ids
rl <- readLines(motifsCheck)
homerIds <- rl[grep('>.*', rl)]
homerIds <- unlist(lapply(homerIds, 
                function (x) {strsplit(x, '\t')[[1]][2]}))
homerIds <- unlist(lapply(homerIds,
                function (x) {strsplit(x, '/')[[1]][1]}))

# Change weird characters (mandatory to use TOBIAS)
homerIds <- unlist(lapply(homerIds,
                function (x) {gsub(' ', '', x)[[1]]}))
homerIds <- unlist(lapply(homerIds,
                function (x) {gsub(':', '-', x)[[1]]}))
homerIds <- unlist(lapply(homerIds,
                function (x) {gsub('\\\?', '', x)[[1]]}))
homerIds <- unlist(lapply(homerIds,
                function (x) {gsub('\\\|', '--', x)[[1]]}))
homerIds <- unlist(lapply(homerIds,
                function (x) {gsub('\\\+', '--', x)[[1]]}))
homerIds <- unlist(lapply(homerIds,
                function (x) {gsub('\\\(', '_', x)[[1]]}))
homerIds <- unlist(lapply(homerIds,
                function (x) {gsub('\\\)', '', x)[[1]]}))
homerIds <- unlist(lapply(homerIds,
                function (x) {gsub(',', '-', x)[[1]]}))
homerIds <- unlist(lapply(homerIds,
                function (x) {gsub('_-', '_', x)[[1]]}))


# Look for duplicates
n_occur <- data.frame(table(homerIds))
for (id1 in n_occur[n_occur\$Freq > 1,'homerIds']) {
    pos <- homerIds == id1
    homerIds[pos] <- paste0(id1, paste0('_', 1:sum(pos)))
}

homerIds <- paste0('>', homerIds)
# Then we open file
motifs = read_homer(motifsCheck, skip = 0)
# rename motifs as they are in homer
for (i in 1:length(homerIds)) {
    motifs[[i]]['name'] = homerIds[i]
}

write_matrix(motifs, paste0(fileout, '.motifs'), positions = 'columns', 
            rownames = FALSE,
        'jaspar', sep = '', headers = TRUE, overwrite = TRUE,
        append = FALSE)
"
echo "$Rcommand" | $R --no-save
# Multiply matrix values by 100 so that motifs are visible in summary
# Python
finalMotifs=${outdir}/motifs/known_jaspar_t100.motifs
PyCommand="
with open('${finalMotifs}', 'w') as fout:
    with open('${outdir}/motifs/known_jaspar.motifs', 'r') as fin:
        for line in fin:
            if line[0] != '>' and line != '\n':
                line = [str(round(float(l) * 100, 4)) for l in line.split()]
                line = '\t'.join(line) + '\n'
            
            fout.write(f'{line}')
"
echo "$PyCommand" | python
