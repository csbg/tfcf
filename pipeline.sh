Rscript src/FULLINT_01_01_Integration.R
Rscript src/FULLINT_05_01_SingleR.R
Rscript src/FULLINT_06_01_CytoTRACE.R
Rscript src/FULLINT_10_01_BasicAnalysis.R

# "leukemia" "in vivo" "in vitro"
for tissue in "in vivo" "in vitro"; do
  echo "$tissue"
	Rscript src/FULLINT_10_01_BasicAnalysis.R "$tissue"
done