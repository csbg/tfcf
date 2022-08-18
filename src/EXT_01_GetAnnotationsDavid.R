source("src/00_init.R")
out <- dirout("EXT_01_GetAnnotationsDavid/")



# READ DATA ---------------------------------------------------------------
x <- fread("metadata/TFCF_Annotations.tsv")
x[,GENE := gsub("^(.)", "\\U\\1", tolower(GENE), perl = TRUE)]


# TFCF MAPPING ------------------------------------------------------------
write.tsv(x[,c("GENE", "Category")], out("Manual.TFCF.mapping.tsv"))



# COMPLEXES ---------------------------------------------------------------
complexes <- x[,c("GENE", "COMPLEX"),with=F]
complexes[, COMPLEX := gsub("&", ",", COMPLEX)]
complexes[, COMPLEX := gsub("and", ",", COMPLEX)]
complexes[, COMPLEX := gsub("rna complexes:", "", COMPLEX, ignore.case = T)]
complexes[, COMPLEX := gsub(" complex", "", COMPLEX, ignore.case = T)]
complexes[, COMPLEX := gsub("\\/(\\d+)$", ".\\1", COMPLEX, ignore.case = T)]
complexes[, COMPLEX := gsub("\\/(\\d+),", ".,", COMPLEX, ignore.case = T)]
complexes[, COMPLEX := gsub("\\(.+?\\)", "", COMPLEX, ignore.case = T)]
complexes[, COMPLEX := gsub("  ", " ", COMPLEX)]
complexes[, COMPLEX := gsub(" ?, ?", ",", COMPLEX)]
complexes[, COMPLEX := gsub("^ ", "", COMPLEX)]
complexes[, COMPLEX := gsub(" $", "", COMPLEX)]
complexes[COMPLEX != ""]$COMPLEX

unique(sort(do.call(c, strsplit(complexes$COMPLEX, ","))))

splitComplexes <- function(x){
  ret <- strsplit(x, ",")[[1]]
  if(any(grepl("\\/", ret))){
    ret <- unique(c(ret, strsplit(ret[grepl("\\/", ret)], "/")[[1]]))  
  }
  return(ret)
}

res <- data.table()
for(i in 1:nrow(complexes)){
  if(complexes[i]$COMPLEX == "") next
  res <- rbind(res, data.table(Gene=complexes[i]$GENE, Complex=splitComplexes(complexes[i]$COMPLEX)))
  # print(complexes[i]$COMPLEX)
  # message(paste(splitComplexes(complexes[i]$COMPLEX), collapse = " --- "))
}

res <- unique(res)
write.tsv(res, out("ManualComplexes.tsv"))




# READ DATA ---------------------------------------------------------------
#ann2 <- fread("metadata/TFCF_Annotations_v2.tsv")
#write.tsv(res, out("ManualComplexes.tsv"))
