cd $GFS/PROJECTS/TfCf/Analysis/SHARE/
cadfile="$HOME/tmp/cadfile.rcfile"
touch $cadfile
echo "mput *" >> $cadfile
echo "exit" >> $cadfile
cadaver -r $cadfile https://myfiles.sbg.ac.at/remote.php/dav/files/b1075184/papers/TFCF/DataDrop_NF/
rm -rf *
rm $cadfile