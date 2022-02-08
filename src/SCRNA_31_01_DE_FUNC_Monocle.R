require(monocle3)

neb.file <- out("Nebular.RDS")

# DE for GUIDES -----------------------------------------
obj.de <- monocle.obj#[,monocle.obj$clusterDE %in% c("1", "21")]

# Keep only KO and NTCs
obj.de <- obj.de[,obj.de$mixscape_class.global %in% c("KO", "NTC")]
obj.de <- obj.de[,obj.de$clusterDE %in% names(which(table(obj.de$clusterDE) > 30))]

# Remove lowly expressed genes
obj.de <- obj.de[Matrix::rowSums(counts(obj.de)) > 20,]

# Order by sample
obj.de <- obj.de[,order(obj.de$sample_broad)]
# ggplot(data.table(sample=obj.de$sample, i=1:ncol(obj.de)), aes(x=i, y=sample)) + geom_point()

# prepare object for DE
obj.de$guideDE <- relevel(factor(gsub("_.+", "", obj.de$guide)), ref = "NTC")


# Run nebula --------------------------------------------------------------

# run nebula
gene_fits <- fit_models(obj.de, model_formula_str = "~guideDE + clusterDE + sample_broad")

res <- coefficient_table(gene_fits)
res <- data.table(select(res, -c("model", "model_summary")))
resGuides <- res[grepl("guideDE", term)][status == "OK"][, -c("status", "std_err", "test_val", "model_component")]
resGuides[, guide := gsub("^guideDE", "", term)]
resGuides[, guide := gsub("_.+$", "", guide)]
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