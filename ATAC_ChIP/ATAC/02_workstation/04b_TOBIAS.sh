#!/bin/bash
# -*- ENCODING: UTF-8 -*-

# HOW TO RUN ME
# bash /PATH/TO/SCRIPTS/00b_TOBIAS.sh

# This scripts runs TOBIAS comparing datasets with a specific
# string with others. It also discerns comparison by batches
# specified in the name after the first _
# cell-condition_batch-[extra]_...

# Path to TOBIAS bin
TOBIAS=/PATH/TO/TOBIAS/bin/TOBIAS
# Path to R bin
R=/PATH/TO/R/bin/R
# Path to Uropa bin
uropa="/PATH/TO/bin/uropa"
# best factors script
#bestFactors="/PATH/TO/get_best_bingfactors.py"
# base merged bams path (replicates are not used separated)
bambase=/PATH/TO/bams/merged
# base peaks path
peakbase=/PATH/TO/peaks
# base output dir
outdirbase=/PATH/TO/OUTPUT
# Path to motifs to be check (In format suited for TOBIAS, 00a_...)
motifsCheck=/PATH/TO/OUTDIR/motifs/known_jaspar_t100.motifs
# Path to reference genome
refGenome=/PATH/TO/genomes/mm10_reordered/mm10.reordered.fa
# Cell type of interest
cell=LSK
# Specify control conditions string in name
controlCondition=NTC5
# define peak type
peaktype=broadPeak
# path for the location of the pipeline scripts
scriptsPath="/PATH/TO/SCRIPTS"
# Number of CPU to use
nCPU=16
# Annotation file
gtfFile="/PATH/TO/annotations/Mus_musculus.GRCm38.99.chrYes.gtf"
# Clean heavy files? "yes" "no"
clean="yes"

##===============================================================================

# function to check if the given first file doesnt exist or is older than 
# the second input file
fileNotExistOrOlder () {
    # check if the file exists of it was created with a previous bam version 
    analyse="no"
    if [ ! -e $1 ]; then
        analyse="yes"
    # only proceed if the output file is older than the bam file
    # in this way if we resequenced and kept the name the analysis 
    # will be repeated
    else
        for tfile in $2; do
            if [[ $1 -ot ${tfile} ]] ; then
                analyse="yes"
                echo $1" older than "${tfile}
            fi
        done
    fi
}

##===============================================================================

outdir="${outdirbase}/${peaktype}"
if [ ! -e ${outdir}/${cell} ]; then
	mkdir -p ${outdir}/${cell}
fi

# get NTC and all KO
allbams=$(find -L ${bambase}/${cell}*bam)

# lets split by batches
batches=$(for ba in $allbams; do
    ba=$(basename $ba)
    ba=(${ba//_/ }); 
    ba=${ba[1]};
    ba=(${ba//-/ }); 
    ba=${ba[0]};
    echo $ba;
done | sort | uniq)
 
for batch in ATAC6; do

    echo "Processing batch ${batch}"

    NTC_bam=$(echo $allbams| tr ' ' '\n' | grep "${cell}-${controlCondition}_${batch}-\|${cell}-${controlCondition}_${batch}_") 
    ko_bams=$(echo $allbams| tr ' ' '\n' | grep -v "${cell}-${controlCondition}_") 
    ko_bams=$(echo $ko_bams| tr ' ' '\n' | grep "${cell}-.*_${batch}-\|${cell}-.*_${batch}_") 
    
    # Proceed only if we have both control and ko
    continue="yes"
    if [ -z "$NTC_bam" ] ; then 
        echo "No control ${controlCondition} found in ${batch}"
        continue="no"
    fi
    if [ -z "$ko_bams" ] ; then 
        echo "No KO found in ${batch}"
        continue="no"
    fi

    if [ ${continue} == "yes" ]; then

        for ko_bam in ${ko_bams}; do

            echo "Processing KO ${ko_bam}"

            prefix_ko=$(basename ${ko_bam} | sed 's/.sort.rmdup.*.bam//g')
            outKo=${outdir}/${cell}/${prefix_ko}
            prefix_ntc=$(basename ${NTC_bam} | sed 's/.sort.rmdup.*.bam//g')
            outNTC=${outdir}/${cell}/${prefix_ko}/${prefix_ntc}

            mkdir -p ${outKo}


            #########################
            # CONSENSUS PEAKS: Get consensus peaks for same cell type
            #########################

            # get all file names to use
            consensusPeakBed=${outdir}/${cell}/${prefix_ko}/${cell}.bed
            annotPeaks=${consensusPeakBed::-4}  
            consensusPeakBed_annot=${annotPeaks}_annotated.bed

            # if we did the analysis and deleted intermediate files we dont want to redoo all
            if [ ! -e "${outKo}/${prefix_ko}_footprints.bw" ]; then

                # check if the file exists of it was created with a previous bam version 
                fileNotExistOrOlder "${consensusPeakBed_annot}" "${NTC_bam} ${ko_bam}"
                # this outputs analyse as yes or no in lowercase
                if [ ! -e ${consensusPeakBed_annot} ]; then
                    echo "Processing consensus peaks for ${cell}-${controlCondition}_${batch}_${prefix_ko}"
                
                #if [[ ${analyse} == "yes" ]]; then
                    # Lets create a consensus peak coordinates file for all the conditions of same cell
                    allPeaks=$(find -L ${peakbase}/${cell}*${peaktype})
                    if [ ${peaktype} == "narrowPeak" ]; then
                        mergecols=`seq 2 10 | tr '\n' ','`
                        expandparam='--is_narrow_peak'
                    elif [ ${peaktype} == "broadPeak" ]; then
                        mergecols=`seq 2 9 | tr '\n' ','`
                        expandparam=''
                    fi
                    fileLabels=$(for f in $allPeaks; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',')

                    sort -T '.' -k1,1 -k2,2n ${allPeaks} \
                                    | mergeBed -c $mergecols -o collapse > ${outdir}/${cell}/${prefix_ko}/${cell}_consenus.txt
                    python ${scriptsPath}/subscripts/01_NR_macs2_merged_expand.py ${outdir}/${cell}/${prefix_ko}/${cell}_consenus.txt \
                                        ${fileLabels} \
                                        ${outdir}/${cell}/${prefix_ko}/${cell}_consenus.boolean.txt \
                                        $expandparam

                    
                    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' \
                        ${outdir}/${cell}/${prefix_ko}/${cell}_consenus.boolean.txt > ${consensusPeakBed}
                    rm ${outdir}/${cell}/${prefix_ko}/${cell}_consenus.txt
                    rm ${outdir}/${cell}/${prefix_ko}/${cell}_consenus.boolean.txt

                    # Annotate the peaks
                    $uropa --bed ${consensusPeakBed} --gtf ${gtfFile} \
                                --show_attributes gene_id gene_name --feature_anchor start \
                                --distance 20000 10000 --feature gene --outdir ${outdir}/${cell}/${prefix_ko}
                            
                    cut -f 1-6,16-17 ${annotPeaks}_finalhits.txt | head -n 1 > ${annotPeaks}_annotated_header.txt
                    cut -f 1-6,16-17 ${annotPeaks}_finalhits.txt | tail -n +2 > ${annotPeaks}_annotated.bed
                else
                    echo "Consensus peaks for ${cell}-${controlCondition}_${batch}_${prefix_ko} already processed"
                fi

                #########################
                # Get footprints bigwigs
                #########################

                # check if the file exists of it was created with a previous bam version 
                fileNotExistOrOlder "${outKo}/${prefix_ko}_footprints.bw" "${bambase}/${prefix_ko}.sort.rmdup.rmblackls.rmchr.Tn5.bam"
                # this outputs analyse as yes or no in lowercase
                if [ ! -e "${outKo}/${prefix_ko}_footprints.bw" ]; then
                #if [[ ${analyse} == "yes" ]]; then
                    echo "Processing footprints for ${cell}-${controlCondition}_${batch}_${prefix_ko}"

                    ### In NTC
                    bamfile=${bambase}/${prefix_ntc}.sort.rmdup.rmblackls.rmchr.Tn5.bam
                    $TOBIAS ATACorrect --bam ${bamfile} \
                            --genome ${refGenome} \
                            --peaks ${consensusPeakBed} --read_shift 0 0 --prefix ${prefix_ntc} \
                            --outdir ${outNTC} --cores ${nCPU}

                    #Calculate footprint scores per condition*/
                    correctedbw=${outNTC}/${prefix_ntc}_corrected.bw
                    $TOBIAS ScoreBigwig --signal ${correctedbw} --regions ${consensusPeakBed} \
                                        --cores ${nCPU} --output ${outNTC}/${prefix_ntc}_footprints.bw &> \
                                        ${outNTC}/${prefix_ntc}_footprints.log
                        



                    ### In KO

                    bamfile=${bambase}/${prefix_ko}.sort.rmdup.rmblackls.rmchr.Tn5.bam
                    $TOBIAS ATACorrect --bam ${bamfile} \
                            --genome ${refGenome} \
                            --peaks ${consensusPeakBed} --read_shift 0 0 --prefix ${prefix_ko} \
                            --outdir ${outKo} --cores ${nCPU}
                    #TOBIAS ATACorrect ${atacorrect} -b $allmerged -g $fasta -p $allbed --cores 99 ${atacorrect} --blacklist $blacklist --prefix $condition &> ${condition}_atacorrect.log
                    #Calculate footprint scores per condition*/
                    correctedbw=${outKo}/${prefix_ko}_corrected.bw
                    $TOBIAS ScoreBigwig --signal ${correctedbw} --regions ${consensusPeakBed} \
                                    --cores ${nCPU} --output ${outKo}/${prefix_ko}_footprints.bw &> \
                                    ${outKo}/${prefix_ko}_footprints.log
                    #TOBIAS ScoreBigwig --signal ${correctedbw} ${footprinting} --regions $allmerged ${footprinting} --cores 99 --output ${condition}_footprints.bw &> ${condition}_footprinting.log
                else
                    echo "Footprints for ${cell}-${controlCondition}_${batch}_${prefix_ko} already processed"
                fi


                #########################
                # Get differentially bound motifs 
                #########################
                #outprefix1=(${prefix_ko//_/ }); 
                #outprefix1=${outprefix1[0]}; 
                #outprefix2=(${prefix_ntc//_/ }); 
                #outprefix2=${outprefix2[0]};
                outprefix=${prefix_ko}_${prefix_ntc}

                # check if the file exists of it was created with a previous bam version 
                fileNotExistOrOlder "${outKo}/bindetect_output/bindetect-${outprefix}_figures.pdf" "${outKo}/${prefix_ko}_footprints.bw"
                # this outputs analyse as yes or no in lowercase
                if [ ! -e "${outKo}/bindetect_output/bindetect-${outprefix}_figures.pdf" ]; then
                #if [[ ${analyse} == "yes" ]]; then
                    echo "Processing motifs for ${cell}-${controlCondition}_${batch}_${prefix_ko}"
                    #Estimate bound sites from scored */
                    # here we also can do differential
                    # It works as it should when you use raw PFM motif files

                    $TOBIAS BINDetect --motifs ${motifsCheck} \
                                    --signals ${outKo}/${prefix_ko}_footprints.bw ${outNTC}/${prefix_ntc}_footprints.bw \
                                    --genome ${refGenome} --peaks ${consensusPeakBed_annot} \
                                    --cond-names ${prefix_ko} ${prefix_ntc} --cores ${nCPU} \
                                    --prefix bindetect-${outprefix} \
                                    --naming 'name' \
                                    --outdir ${outKo}/bindetect_output &> ${outKo}/${outprefix}_bindetect.log
                    #TOBIAS BINDetect --motifs ${motifsCheck} --signals $footprints ${bindetect} --genome $fasta --peaks  $annotated_headerbed --peak_header $annotated_header --cores 99 --cond_names $condition --outdir TFBS &> bindetect.log
                else
                    echo "Motifs for ${cell}-${controlCondition}_${batch}_${prefix_ko} already processed"
                fi

                #########################
                # Networks analysis
                # very heavy files as output
                #########################

                # # check if the file exists of it was created with a previous bam version 
                # fileNotExistOrOlder "${outKo}/bindetect_output/network_${prefix_ko}/edges.txt" "${outKo}/bindetect_output/bindetect-${outprefix}_figures.pdf"
                # # this outputs analyse as yes or no in lowercase
                # if [ ! -e "${outKo}/bindetect_output/network_${prefix_ko}/edges.txt" ]; then
                # #if [[ ${analyse} == "yes" ]]; then
                #     echo "Processing networks for ${cell}-${controlCondition}_${batch}_${prefix_ko}"
                #     ### Network analysis
                #     # You can get the TF IDs from the motif file and add them here to get the gene symbol
                #     # https://www.genenames.org/tools/multi-symbol-checker/
                #     origin="/home/julen/programas/annotations/motifTFlink/homer_vertebrates_known_Motif-hgnc_TOBIAS_valid.csv"
                #     # In ko ATAC
                #     mkdir -p ${outKo}/bindetect_output/network_${prefix_ko}
                #     $TOBIAS CreateNetwork --TFBS ${outKo}/bindetect_output/*/beds/*${prefix_ko}_bound.bed \
                #         --origin ${origin} --outdir ${outKo}/bindetect_output/network_${prefix_ko}
                #     # in control ATAC
                #     mkdir -p ${outKo}/bindetect_output/network_${prefix_ntc}
                #     $TOBIAS CreateNetwork --TFBS ${outKo}/bindetect_output/*/beds/*${prefix_ntc}_bound.bed \
                #         --origin ${origin} --outdir ${outKo}/bindetect_output/network_${prefix_ntc}
                # else
                #     echo "Networks for ${cell}-${controlCondition}_${batch}_${prefix_ko} already processed"
                # fi

                if [ $clean == "yes" ]; then
                    mv ${outKo}/bindetect_output/*txt ${outKo}/; 
                    mv ${outKo}/bindetect_output/*pdf ${outKo}/; 
                    mv ${outKo}/bindetect_output/*xlsx ${outKo}/;
                    mv ${outKo}/bindetect_output/*html ${outKo}/; 
                    #mv ${outKo}/bindetect_output/network* ${outKo}/; 
                    rm -r ${outKo}/bindetect_output/; 
                    rm -r ${outNTC}/; 
                    rm ${outKo}/*uncorrected.bw
                    rm ${outKo}/*expected.bw
                    rm ${outKo}/*corrected.bw
                    rm ${outKo}/*bias.bw
                    #rm *bw
                fi
            fi
        done
    fi
done

## Appart from all
#$TOBIAS ClusterMotifs --motifs ${motifsCheck} --threshold 0.3 --type png --outdir /scratch/julen/ATAC/allData/02_firstATAC/TOBIAS/motifClustering

# get a selection of the best TF
# ${bestFactors} -in ${outKo}/bindetect_output/bindetect-${outprefix}_results.txt \
#         -filter null -o ${outKo}/bindetect_output/bindetect-${outprefix}_results_best_hits.xlsx

#Join subset of bound TFs*/
#cat $bed | cut -f1-4 | sort -k1,1 -k2,2n > all_${condition}_bound.bed; igvtools index all_${condition}_bound.bed 2>&1;

# delPath=/PATH/TO/OUTDIR/TOBIAS/broadPeak
# for cell in */; do 
#     for i in ${cell}/${cell}-*/; do 
#         cd ${delPath}/${i}; 
#         mv bindetect_output/*txt .; 
#         mv bindetect_output/*pdf .; 
#         mv bindetect_output/*xlsx .;
#         mv bindetect_output/*html .; 
#         mv bindetect_output/network* . ; 
#         rm -r bindetect_output/; 
#         rm -r ${cell}*/; 
#         rm *bw
#     done
# done
