cd $GFS/PROJECTS/TfCf/

mkdir CRUK

for f in CITESEQ1 CITESEQ2 ECCITE1 ECCITE2 ECCITE5 ECCITE7 ECCITE7_2;do
	zip -r CRUK/$f.zip RawData/$f
done

cadfile="$HOME/tmp/cadfile.rcfile"
touch $cadfile
echo "mput *" >> $cadfile
echo "exit" >> $cadfile
cd CRUK/
cadaver -r $cadfile https://myfiles.sbg.ac.at/remote.php/dav/files/b1075184/papers/2022_TFCF/DataDrop_NF/CRUK/
rm -rf *
rm $cadfile