

require(monocle3)
require(Seurat)

x <- Read10X_h5("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SHARE/LINES_ECCITE4_Cas9_LINES.h5", use.names = TRUE)
x <- x$`Gene Expression`
str(row.names(x))
tail(row.names(x))
grep("line", row.names(x), value = TRUE, ignore.case = TRUE)
