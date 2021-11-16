# Rscript src/FULLINT_01_01_Integration.R
Rscript src/FULLINT_05_01_SingleR.R
# Rscript src/FULLINT_06_01_CytoTRACE.R
Rscript src/FULLINT_10_01_BasicAnalysis.R

for tissue in "leukemia" "in vivo" "in vitro"; do
  echo "$tissue"
	Rscript src/FULLINT_10_01_BasicAnalysis.R "$tissue"
done