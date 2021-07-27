source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
baseDir <- "PPI_01_CollectData/"
out <- dirout(baseDir)

require(readxl)

inDir <- dirout_load("PPI_00_getData/")

list.files(inDir(""))

# HUMAN - MOUSE mapping ---------------------------------------------------
# Run on July 27, 2021
# Dataset
# - Mouse genes (GRCm39)
# Filters
# - [None selected]
# Attributes
# - Gene stable ID
# - Human gene stable ID
# - Human homology type
# - Gene name
# - Human gene name
hm <- fread(inDir("BioMart_Human_Mouse_2021_07_27.txt"), check.names = T)
hm <- setNames(hm[Human.homology.type == "ortholog_one2one",c("Gene.name", "Human.gene.name"),with=F], c("Mouse", "Human"))
hm <- unique(hm)
hm <- hm[-unique(c(which(duplicated(hm$Mouse)), which(duplicated(hm$Human)))),]
stopifnot(!any(duplicated(hm$Mouse)))
stopifnot(!any(duplicated(hm$Human)))



# UNIPROT ID MAP ----------------------------------------------------------
uniprot.file <- out("Uniprot.tsv")
if(!file.exists(uniprot.file)){
  uniprot.id.map <- fread("https://www.uniprot.org/uniprot/?query=reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22&format=tab&columns=id,entry%20name,genes&sort=score")
  names(uniprot.id.map) <- make.names(names(uniprot.id.map))
  uniprot.id.map[,Gene := gsub(" .+$", "", Gene.names)]
  uniprot.id.map <- unique(setNames(uniprot.id.map[,c("Entry.name", "Gene", "Entry")], c("ID", "Gene", "ACC")))
  write.tsv(uniprot.id.map, uniprot.file)
} else {
  uniprot.id.map <- fread(uniprot.file)
}


# HIPPIE ------------------------------------------------------------------
hippie <- fread(inDir("hippie.txt"))
names(hippie.id.map) <- make.names(names(hippie.id.map))
hippie <- merge(hippie, uniprot.id.map, by.x="V1", by.y="ID")
hippie <- merge(hippie, uniprot.id.map, by.x="V3", by.y="ID")
hippie <- setNames(hippie[,c("Gene.x", "Gene.y", "V5"), with=F], c("A","B","Score"))
hippie <- unique(hippie)
hippie <- data.table(hippie, db="HIPPIE", dataset="HIPPIE", organism="Human")


# HUMAP2 ------------------------------------------------------------------
humap2 <- fread(inDir("hu.MAP2.Drew.2021.MolSysBio"))
humap2 <- humap2[V3 > 0.5]
message("MISSING SOME INFORMATION HERE - UNIPROT ACCESSIONS")
#humap2[grepl("\\|", V1) | ]
humap2 <- setNames(unique(humap2), c("A","B","Score"))
humap2 <- data.table(humap2, db="hu.Map2", dataset="hu.Map2", organism="Human")


# CFMS --------------------------------------------------------------------
cfms <- fread(inDir("CFMS.Skinnider.2021.NatureMethods.tsv"))
ggplot(cfms, aes(x=score, y=precision)) + geom_hex()
ggsave(out("CFMS_ScorevsPrecision.pdf"),w=5,h=5)
cfms <- setNames(cfms[,c("protein_A", "protein_B", "score"),with=F], c("A","B","Score"))
cfms <- data.table(cfms, db="CF-MS", dataset="CF-MS", organism="Human")


# PCP SILAM ---------------------------------------------------------------
pcp <- data.table(readxl::read_xlsx(inDir("PCP.SILAM.Skinnider.2021.Cell.xlsx"), sheet = 1))
pcp.list <- split(setNames(pcp[,-"Tissue"], c("A","B", "Score")), pcp$Tissue)
pcp.list <- lapply(names(pcp.list), function(nam) data.table(pcp.list[[nam]], db="PCP SILAM", dataset=nam, organism="Mouse"))


# CORUM -------------------------------------------------------------------
corumCore <- fread(cmd = paste0('unzip -cq ', inDir("CorumCore.txt.zip")), check.names = T)
corumAll <- fread(cmd = paste0('unzip -cq ', inDir("CorumAll.txt.zip")), check.names = T)
interactions.from.corum <- function(x){
  x[,subunits.simple := gsub(" ", ";", gsub(",", ";", subunits.Gene.name.syn.))]
  gg <- unique(do.call(c, strsplit(x$subunits.simple, ";")))
  adjm <- matrix(0, length(gg), length(gg))
  row.names(adjm) <- gg
  colnames(adjm) <- gg
  i <- 1
  for(i in 1:nrow(x)){
    gs <- unique(strsplit(x[i, subunits.simple], ";")[[1]])
    gs <- gs[gs != ""]
    adjm[gs, gs] <- 1
  }
  res <- unique(melt(data.table(data.frame(adjm),keep.rownames = TRUE), id.vars ="rn")[value == 1])
  res <- setNames(res, c("A","B", "Score"))
  return(res)
}
corum.list <- list(
  data.table(interactions.from.corum(corumCore[Organism == "Human"]), db="Corum", dataset="Core", organism="Human"),
  data.table(interactions.from.corum(corumCore[Organism == "Mouse"]), db="Corum", dataset="Core", organism="Mouse"),
  data.table(interactions.from.corum(corumAll[Organism == "Human"]), db="Corum", dataset="All", organism="Human"),
  data.table(interactions.from.corum(corumAll[Organism == "Mouse"]), db="Corum", dataset="All", organism="Mouse")
)


# COMBINE -------------------------------------------------------------------------
all.interactions <- data.table()
all.interactions <- rbind(all.interactions, hippie)
all.interactions <- rbind(all.interactions, humap2)
all.interactions <- rbind(all.interactions, cfms)
all.interactions <- rbind(all.interactions, do.call(rbind, pcp.list))
all.interactions <- rbind(all.interactions, do.call(rbind, corum.list))
all.interactions <- rbind(all.interactions, hippie)
all.interactions.orig <- copy(all.interactions)

# Get human names
all.interactions <- merge(all.interactions, hm, by.x="A", by.y="Mouse", all.x=TRUE)
stopifnot(nrow(all.interactions) == nrow(all.interactions.orig))
all.interactions <- merge(all.interactions, hm, by.x="B", by.y="Mouse", all.x=TRUE, suffixes = c("_A", "_B"))
stopifnot(nrow(all.interactions) == nrow(all.interactions.orig))

# Get mouse names
all.interactions <- merge(all.interactions, hm, by.x="A", by.y="Human", all.x=TRUE)
stopifnot(nrow(all.interactions) == nrow(all.interactions.orig))
all.interactions <- merge(all.interactions, hm, by.x="B", by.y="Human", all.x=TRUE, suffixes = c("_A", "_B"))
stopifnot(nrow(all.interactions) == nrow(all.interactions.orig))

# Finalize names
all.interactions[organism == "Mouse", Mouse_A := A]
all.interactions[organism == "Mouse", Mouse_B := B]
all.interactions[organism == "Human", Human_A := A]
all.interactions[organism == "Human", Human_B := B]
#write.tsv(all.interactions, out("PPI.tsv"))


# Identify interesting pairs ----------------------------------------------
mds <- fread(dirout_load("POOLED_10_IndividualAnalysis")("Correlation_","hits","_MDS.tsv"))
mds.dist <- as.matrix(mds[,2:3])
row.names(mds.dist) <- mds$rn
mds.dist <- as.matrix(dist(mds.dist))
diag(mds.dist) <- max(mds.dist)
ii <- unique(data.table(t(apply(which(mds.dist <= sort(mds.dist[upper.tri(mds.dist)], decreasing = FALSE)[200], arr.ind = TRUE), 1, sort))))
pairs <- data.table(A=mds$rn[ii$V1], B=mds$rn[ii$V2])
gg <- unique(c(pairs$A, pairs$B))

all.interactions.small <- all.interactions[Mouse_A %in% gg & Mouse_B %in% gg]

# Plots pairs -------------------------------------------------------------
res <- data.table()
for(pi in 1:nrow(pairs)){
  pa=pairs[pi]$A
  pb=pairs[pi]$B
  res <- rbind(res, data.table(all.interactions.small[(Mouse_A == pa & Mouse_B == pb) | (Mouse_A == pb & Mouse_B == pa)], pair=paste(pa, pb)))
}
res <- res[!is.na(dataset)]

p.ppi <- ggplot(res, aes(y=pair, x=dataset, color=Score)) +
  geom_point() +
  theme_bw(12) + 
  facet_grid(. ~ dataset, switch = "x", space = "free", scales = "free") +
  scale_color_gradient(low="blue", high="red", limits=c(0,1))
ggsave(out("PPI.result.pdf"), w=6,  h=15, plot=p.ppi + xRot())

# ADD DEPMAP ------------------------------------------------------------------
depmap <- read.table(inDir("depmap.CRISPR.csv"), row.names = 1, sep=",", header = TRUE)
str(depmap)
stopifnot(all(sapply(depmap, is.numeric)))
depmap <- as.matrix(depmap)
depmap.ann <- fread(inDir("depmap.ann.csv"))
depmap <- depmap[row.names(depmap) %in% depmap.ann$DepMap_ID,]
depmap.ann <- depmap.ann[DepMap_ID %in% row.names(depmap)]
depmap.clean <- depmap
colnames(depmap.clean) <- gsub("\\..+", "", colnames(depmap.clean))

depmap.groups <- with(depmap.ann, split(DepMap_ID, lineage))
depmap.groups <- depmap.groups[sapply(depmap.groups, length) > 20]
depmap.groups <- depmap.groups[names(depmap.groups) != "unknown"]

# Calculate specific correlations

dm.res <- get.dm.cor(depmap.clean, pairs, "All")
for(xnam in names(depmap.groups)){
  dm.res <- rbind(dm.res, get.dm.cor(depmap.clean[depmap.groups[[xnam]],], pairs, xnam))
}

get.dm.cor <- function(dm, pairs, name){
  gg <- unique(c(pairs$A, pairs$B))
  cMT <- cor(dm[,hm[Mouse %in% gg & Human %in% colnames(dm)]$Human], use="pairwise.complete.obs")
  cDT <- melt(data.table(cMT, keep.rownames = TRUE), id.vars="rn")
  cDT <- merge(cDT, hm, by.x="rn", by.y="Human")
  cDT <- merge(cDT, hm, by.x="variable", by.y="Human", suffixes=c("_A", "_B"))
  pairs[,pair := paste(A, B)]
  return(data.table(merge(cDT, pairs, by.x=c("Mouse_A", "Mouse_B"), by.y=c("A","B")), db="DepMap", dataset=name, organism="Human"))
}

p.dm <- ggplot(dm.res[pair %in% res$pair], aes(y=pair, x=dataset, color=value)) +
  geom_point() +
  theme_bw(12) + 
  facet_grid(. ~ db, switch = "x", space = "free", scales = "free") +
  #xRot() +
  scale_color_gradient2(low="blue", mid="white", high="red", limits=c(-1,1))
ggsave(out("DepMap.result.pdf"), w=6,  h=15, plot=p.dm+xRot())


p <- gridExtra::grid.arrange(
  p.ppi,
  p.dm, 
  nrow=1, ncol=2, widths=c(3,4))
ggsave(out("Results_combined.pdf"), h=10, w=9, plot=p)



