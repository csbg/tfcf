source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
baseDir <- "PPI_01_CollectData/"
out <- dirout(baseDir)

require(readxl)
require(igraph)

inDir <- dirout_load("PPI_00_getData/")

list.files(inDir(""))


# Genes of interest -------------------------------------------------------
m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
(GENES.OF.INTEREST <- unique(gsub("_\\d+$", "", row.names(m))))

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



# Walktrap function -------------------------------------------------------
walktrap.probabilities <- function(g, steps=5){
  adj <- get.adjacency(g)
  diag(adj) <- 0
  degree <- rowSums(adj)
  degree[degree == 0] <- 1
  prob <- adj/degree
  stopifnot(all(round(rowSums(prob),4)[rowSums(adj) > 0]==1))
  stopifnot(all(round(rowSums(prob),4)[rowSums(adj) == 0]==0))
  #table(Matrix::rowSums(prob))
  probs <- list()
  probs[[1]] <- prob
  for(i in 2:steps){
    print(paste("Step", i))
    probs[[i]] <- probs[[i-1]] %*% prob
  }
  for(i in 1:steps){
    print(paste("Checking", i))
    stopifnot(all(round(rowSums(probs[[i]]),4)[rowSums(adj) > 0]==1))
    stopifnot(all(round(rowSums(probs[[i]]),4)[rowSums(adj) == 0]==0))
  }
  probs
}
# Making sure this works
# https://towardsdatascience.com/similarity-of-documents-with-random-walks-98f94fd2c76c
# https://arxiv.org/pdf/physics/0512106.pdf
g <- graph.edgelist(as.matrix(data.table(A=paste0("D", c(1,1,2,3,3)), B=paste0("T", c(1,2,1,1,2)))), directed = FALSE)
plot(g)
adj <- get.adjacency(g)
diag(adj) <- 0
degree <- rowSums(adj)
degree[degree == 0] <- 1
prob <- adj/degree
t(prob)
prob %*% prob


# HIPPIE ------------------------------------------------------------------
hippie <- fread(inDir("hippie.txt"))
hippie <- merge(hippie, uniprot.id.map, by.x="V1", by.y="ID")
hippie <- merge(hippie, uniprot.id.map, by.x="V3", by.y="ID")
hippie <- setNames(hippie[,c("Gene.x", "Gene.y", "V5"), with=F], c("A","B","Score"))
hippie <- unique(hippie)
hippie <- data.table(hippie, db="HIPPIE", dataset="HIPPIE", organism="Human")


# HIPPIE DISTANCE ---------------------------------------------------------
g <- graph.edgelist(as.matrix(hippie[Score > 0.7][,c("A","B"),with=F]), directed = FALSE)
g <- simplify(g)
goi <- intersect(hm[Mouse %in% GENES.OF.INTEREST]$Human, V(g)$name)

# distances (not very helpful)
dd <- distances(g, goi, goi)
table(dd)

# calculate walktrap probabilities
steps=8
probs <- walktrap.probabilities(g, steps = steps)

# Look at differences
i <- 1
j <- 2
for(i in 1:(steps-1)){
  for(j in (i+1):steps){
    message(i, " vs ", j)
    print(sum(abs(probs[[i]][goi,] - probs[[j]][goi,])))
  }
}

# Get similarity
sim <- cor(t(as.matrix(probs[[4]][goi,])))
diag(sim) <- NA
cleanDev(); pdf(out("HIPPIE_DISTANCES.pdf"), w=30,h=30)
pheatmap(sim, fontsize_row = 5)
dev.off()

# similarity into a table
sum(sim > 0.7, na.rm = TRUE)
adj2 <- sim
adj2[adj2 < 0.7] <- 0
adj2[is.na(adj2)] <- 0
sum(adj2 > 0.7, na.rm = TRUE)
g <- graph.adjacency(adj2, weighted = "max")
hippie_distance <- data.table(get.edgelist(g), E(g)$max)
for(i in 1:10){
  stopifnot(sim[hippie_distance[i]$V1, hippie_distance[i]$V2] == hippie_distance[i]$V3)
}
hippie_distance <- data.table(setNames(hippie_distance, c("A","B", "Score")), db="HIPPIE", dataset="HIPPIE_Distances", organism="Human")
# hippie[A == "KLF9" & B == "JDP2"]
# hippie[A == "JDP2" & B == "KLF9"]
# hippie[A %in% goi & B %in% goi]
# hippie_distance[grepl("SMARC", A) & grepl("SMARC", B)]
# hippie[grepl("SMARC", A) & grepl("SMARC", B)]



# HUMAP2 ------------------------------------------------------------------
humap2 <- fread(inDir("hu.MAP2.Drew.2021.MolSysBio"))
quantile(humap2$V3)
humap2 <- humap2[V3 > 0.01]

# Clean up protein ids
humap2.clean <- humap2[!grepl("^sp", V2) & !grepl("^sp", V1)]
humap2.corr <- humap2[grepl("^sp", V2) | grepl("^sp", V1)]
i <- 1
for(i in 1:nrow(humap2.corr)){
  if(grepl("^sp", humap2.corr[i]$V2)){
    acc <- gsub("^sp\\|(.+)\\|(.+)?", "\\1", humap2.corr[i]$V2)
    name <- uniprot.id.map[ACC == acc]$Gene
    if(length(name) == 0) name <- NA
    humap2.corr[i, V2 := name]
  }
  if(grepl("^sp", humap2.corr[i]$V1)){
    acc <- gsub("^sp\\|(.+)\\|(.+)?", "\\1", humap2.corr[i]$V1)
    name <- uniprot.id.map[ACC == acc]$Gene
    if(length(name) == 0) name <- NA
    humap2.corr[i, V1 := name]
  }
}
stopifnot(nrow(humap2.corr) + nrow(humap2.clean) == nrow(humap2))
humap2 <- rbind(humap2.clean, humap2.corr[!is.na(V1) & !is.na(V2) & V1 != "" & V2 != ""])
humap2 <- setNames(unique(humap2), c("A","B","Score"))
humap2 <- data.table(humap2, db="hu.Map2", dataset="hu.Map2", organism="Human")


# Hu.Map DISTANCE ---------------------------------------------------------
g <- graph.edgelist(as.matrix(humap2[Score > 0.7][,c("A","B"),with=F]), directed = FALSE)
g <- simplify(g)
goi <- intersect(hm[Mouse %in% GENES.OF.INTEREST]$Human, V(g)$name)

# distances (not very helpful)
dd <- distances(g, goi, goi)
table(dd)

# calculate walktrap probabilities
steps=8
probs <- walktrap.probabilities(g, steps = steps)

# Look at differences
i <- 1
j <- 2
for(i in 1:(steps-1)){
  for(j in (i+1):steps){
    message(i, " vs ", j)
    print(sum(abs(probs[[i]][goi,] - probs[[j]][goi,])))
  }
}

# Get similarity
sim <- cor(t(as.matrix(probs[[4]][goi,])))
diag(sim) <- NA
cleanDev(); pdf(out("HUMAP2_DISTANCES.pdf"), w=30,h=30)
pheatmap(sim, fontsize_row = 5)
dev.off()

# similarity into a table
sum(sim > 0.7, na.rm = TRUE)
adj2 <- sim
adj2[adj2 < 0.7] <- 0
adj2[is.na(adj2)] <- 0
sum(adj2 > 0.7, na.rm = TRUE)
g <- graph.adjacency(adj2, weighted = "max")
humap2_distance <- data.table(get.edgelist(g), E(g)$max)
for(i in 1:10){
  stopifnot(sim[humap2_distance[i]$V1, humap2_distance[i]$V2] == humap2_distance[i]$V3)
}
humap2_distance <- data.table(setNames(humap2_distance, c("A","B", "Score")), db="hu.Map2", dataset="hu.Map2_Distances", organism="Human")



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
  x[,subunits.simple := gsub(" ", ";", gsub(",", ";", subunits.Gene.name.))]
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



# Manual Complex annotation -----------------------------------------------
compexes.David <- fread(dirout_load("EXT_01_GetAnnotationsDavid")("ManualComplexes.tsv"))

x <- compexes.David
interactions.from.david <- function(x){
  gg <- unique(x$Gene)
  adjm <- matrix(0, length(gg), length(gg))
  row.names(adjm) <- gg
  colnames(adjm) <- gg
  cx <- x$Complex[1]
  for(cx in unique(x$Complex)){
    gs <- x[Complex == cx]$Gene
    adjm[gs, gs] <- 1
  }
  res <- unique(melt(data.table(data.frame(adjm),keep.rownames = TRUE), id.vars ="rn")[value == 1])
  res <- setNames(res, c("A","B", "Score"))
  return(res)
}
complexes.manual.list <- list(
  data.table(interactions.from.david(compexes.David), db="Manual", dataset="Manual", organism="Mouse")
)



# DEPMAP ------------------------------------------------------------------
depmap <- read.table(inDir("depmap.CRISPR.csv"), row.names = 1, sep=",", header = TRUE)
str(depmap)
stopifnot(all(sapply(depmap, is.numeric)))
depmap <- as.matrix(depmap)
depmap.ann <- fread(inDir("depmap.ann.csv"))
depmap <- depmap[row.names(depmap) %in% depmap.ann$DepMap_ID,]
depmap.ann <- depmap.ann[DepMap_ID %in% row.names(depmap)]
depmap.clean <- depmap
colnames(depmap.clean) <- gsub("\\..+", "", colnames(depmap.clean))
depmap.genes <- hm[Mouse%in% GENES.OF.INTEREST & Human %in% colnames(depmap.clean)]$Human
depmap.clean <- depmap.clean[,depmap.genes]
str(depmap.clean)

m <- depmap.clean
get.depmap.cor <- function(m, name){
  cMT <- cor(m, use="pairwise.complete.obs")
  cMT[upper.tri(cMT, diag = TRUE)] <- NA
  cDT <- melt(data.table(cMT, keep.rownames = TRUE), id.vars="rn")[!is.na(value)]
  names(cDT) <- c("A", "B", "Score")
  return(data.table(cDT, db="DepMap", dataset=name, organism="Human"))
}


depmap.list <- list(All = get.depmap.cor(m = depmap.clean, name = "All"))
for(tx in unique(depmap.ann[,.N, by="lineage"][N>20]$lineage)){
  depmap.list[[tx]] <- get.depmap.cor(m = depmap.clean[depmap.ann[lineage == tx]$DepMap_ID,], name = tx)
}
# lapply(depmap.list, function(dt) dt[A == "ACTL6A" & B == "AHCTF1"])


# 
# depmap.groups <- with(depmap.ann, split(DepMap_ID, lineage))
# depmap.groups <- depmap.groups[sapply(depmap.groups, length) > 20]
# depmap.groups <- depmap.groups[names(depmap.groups) != "unknown"]

# # Calculate specific correlations
# get.dm.cor <- function(dm, pairs, name){
#   gg <- unique(c(pairs$A, pairs$B))
#   cMT <- cor(dm[,hm[Mouse %in% gg & Human %in% colnames(dm)]$Human], use="pairwise.complete.obs")
#   cDT <- melt(data.table(cMT, keep.rownames = TRUE), id.vars="rn")
#   cDT <- merge(cDT, hm, by.x="rn", by.y="Human")
#   cDT <- merge(cDT, hm, by.x="variable", by.y="Human", suffixes=c("_A", "_B"))
#   pairs[,pair := paste(A, B)]
#   return(data.table(merge(cDT, pairs, by.x=c("Mouse_A", "Mouse_B"), by.y=c("A","B")), db="DepMap", dataset=name, organism="Human"))
# }
# 
# dm.res <- get.dm.cor(depmap.clean, pairs, "All")
# for(xnam in names(depmap.groups)){
#   dm.res <- rbind(dm.res, get.dm.cor(depmap.clean[depmap.groups[[xnam]],], pairs, xnam))
# }
# 
# p.dm <- ggplot(dm.res[pair %in% res$pair], aes(y=pair, x=dataset, color=value)) +
#   geom_point() +
#   theme_bw(12) + 
#   facet_grid(. ~ db, switch = "x", space = "free", scales = "free") +
#   #xRot() +
#   scale_color_gradient2(low="blue", mid="white", high="red", limits=c(-1,1))
# ggsave(out("DepMap.result.pdf"), w=6,  h=15, plot=p.dm+xRot())



# COMBINE -------------------------------------------------------------------------
all.interactions <- data.table()
all.interactions <- rbind(all.interactions, hippie)
all.interactions <- rbind(all.interactions, hippie_distance)
all.interactions <- rbind(all.interactions, humap2)
all.interactions <- rbind(all.interactions, humap2_distance)
all.interactions <- rbind(all.interactions, cfms)
all.interactions <- rbind(all.interactions, do.call(rbind, pcp.list))
all.interactions <- rbind(all.interactions, do.call(rbind, corum.list))
all.interactions <- rbind(all.interactions, do.call(rbind, complexes.manual.list))
all.interactions <- rbind(all.interactions, do.call(rbind, depmap.list))
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

all.interactions[,.N, by=c("dataset", "db", "organism")]

# Finalize names
all.interactions[organism == "Mouse", Mouse_A := A]
all.interactions[organism == "Mouse", Mouse_B := B]
all.interactions[organism == "Human", Human_A := A]
all.interactions[organism == "Human", Human_B := B]
all.interactions$A <- NULL
all.interactions$B <- NULL
save(all.interactions, file=out("All.Interactions.RData"))


interactions <- all.interactions[Mouse_A %in% GENES.OF.INTEREST & Mouse_B %in% GENES.OF.INTEREST]
interactions[,A := Mouse_A]
interactions[,B := Mouse_B]
interactions$Human_A <- NULL
interactions$Human_B <- NULL
interactions$Mouse_A <- NULL
interactions$Mouse_B <- NULL
interactions <- interactions[A != B]
pm <- as.matrix(interactions[,c("A", "B"),with=F])
pms <- t(apply(pm, 1, sort))
pm <- apply(pms, 1, paste, collapse= " - ")
interactions$pair <- pm
res <- data.table()
for(dbx in unique(interactions$db)){
  res <- rbind(res, unique(interactions[db == dbx]))
}
interactions <- res
save(interactions, file=out("GOI.Interactions.RData"))


