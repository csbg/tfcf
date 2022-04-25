require(nebula)
require(monocle3)

neb.file <- out("Nebular.RDS")

# DE for GUIDES -----------------------------------------
obj.de <- monocle.obj#[,monocle.obj$clusterDE %in% c("1", "21")]

# Keep only KO, NPs, and NTCs
obj.de <- obj.de[,obj.de$mixscape_class.global %in% c("KO", "NTC", "NP")]
obj.de <- obj.de[,obj.de$clusterDE %in% names(which(table(obj.de$clusterDE) > 30))]

# Remove lowly expressed genes
obj.de <- obj.de[Matrix::rowSums(counts(obj.de)) > 20,]

# Order by sample
obj.de <- obj.de[,order(obj.de$sample_broad)]
# ggplot(data.table(sample=obj.de$sample, i=1:ncol(obj.de)), aes(x=i, y=sample)) + geom_point()

# prepare annotation for DE
obj.de.ann <- data.frame(
  row.names=colnames(obj.de),
  GuideDE=gsub("_.+$", "", obj.de$guide),
  ClusterDE=obj.de$clusterDE
)

obj.de.ann <- mutate(obj.de.ann, GuideDE = gsub("^Men$", "Men1", GuideDE))
obj.de.ann <- mutate(obj.de.ann, GuideDE = gsub("^Pu.1$", "Spi1", GuideDE))
#obj.de.ann <- filter(obj.de.ann, !GuideDE %in% c("Pu.1", "Spi1"))
table(obj.de.ann$GuideDE)
table(obj.de.ann$ClusterDE)
obj.de.ann$GuideDE <- relevel(factor(obj.de.ann$GuideDE), ref = "NTC")
obj.de <- obj.de[, row.names(obj.de.ann)]
stopifnot(row.names(obj.de.ann) == colnames(obj.de))
str(obj.de.ann)
write.tsv(data.table(obj.de.ann, keep.rownames = TRUE), out("DEG_Annnotation.tsv"))



# Skip if no cells left ---------------------------------------------------
if(ncol(obj.de) < 100) next


# Calculate non-zero values -----------------------------------------------
guide.barcodes <- split(row.names(obj.de.ann), factor(obj.de.ann$GuideDE))
m <- counts(obj.de)
stopifnot(class(m) == "dgCMatrix")
guide.fnz <- lapply(guide.barcodes, function(bcx) tabulate(m[,bcx,drop=F]@i + 1))
guide.fnz <- do.call(cbind, guide.fnz)
row.names(guide.fnz) <- row.names(m)
total.count <- data.table(obj.de.ann)[,.(totalCount = .N), by="GuideDE"]

# Run nebula --------------------------------------------------------------

# Model matrix
mm <- model.matrix(data=obj.de.ann, ~ GuideDE)
#colnames(mm) <- gsub("^(GuideDE.+)$", paste0("\\1_"), colnames(mm))
mm <- mm[,colSums(mm) != 0]

# run nebula
nebRes <- nebula(
  count = counts(obj.de),
  id = obj.de$sample_broad,
  pred = mm
)

# export results
res <- data.table()
for(cx in colnames(mm)){
  res <- rbind(res, data.table(
    term=cx,
    p_value=nebRes$summary[,paste("p", cx, sep="_")],
    se=nebRes$summary[,paste("se", cx, sep="_")],
    estimate=nebRes$summary[,paste("logFC", cx, sep="_")],
    gene_id=nebRes$summary$gene,
    convergence=nebRes$convergence
  ))
}

# finalize results
res[,q_value := p.adjust(p_value, method="BH")]
res[, estimate := estimate / log(2)] # convert to log2 FC
saveRDS(res, neb.file)


#  . Export results -------------------------------------------------------

# Use only reliable results
resGuides <- res[grepl("GuideDE", term)][convergence >= -15]

# Clean up names
resGuides[, guide := gsub("^GuideDE", "", term)]
resGuides[, guide := gsub("_.+$", "", guide)]
resGuides[, guide := gsub("^GuideDE", "", guide)]

# Add counts
resGuides <- merge(
  resGuides, 
  setNames(melt(data.table(guide.fnz, keep.rownames = TRUE), id.vars = "rn"), c("gene_id", "guide", "nonzeroCount")),
  by=c("gene_id", "guide"))
resGuides <- merge(resGuides, total.count, by.x="guide", by.y="GuideDE")

# Clean estimates
resGuides <- resGuides[!is.na(estimate)]
resGuides[, estimate_raw := estimate]
resGuides[, estimate := ifelse(p_value > 0.9, 0, estimate)]
write.tsv(resGuides[q_value < 1][,-"term",with=F], file=out("DEG_Results_Export.tsv"))
write.tsv(resGuides, file=out("DEG_Results_all.tsv"))

#  . Vulcano / p-val distribution -----------------------------------------
ggplot(resGuides, aes(x=estimate, y=-log10(p_value))) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  facet_wrap(~gsub("_", "\n", gsub("\\:tissueDE", "\nInteraction: ", term)), scales = "free", ncol = 5)
ggsave(out("DEG_Vulcano.pdf"), w=15,h=15)

ggplot(resGuides, aes(x=p_value)) + 
  theme_bw(12) +
  geom_histogram() +
  facet_wrap(~gsub("_", "\n", gsub("\\:tissueDE", "\nInteraction: ", term)), scales = "free", ncol = 5)
ggsave(out("DEG_PVal_histogram.pdf"), w=15,h=15)