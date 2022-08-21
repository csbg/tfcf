mkdir -p $HOME/tmp/WeTransfer/

cd $RAWDATA


rm $HOME/tmp/files.txt
find CITESEQ* -type f >> $HOME/tmp/files.txt
find ECCITE1* -type f >> $HOME/tmp/files.txt
find ECCITE2* -type f >> $HOME/tmp/files.txt
find ECCITE5* -type f >> $HOME/tmp/files.txt
find ECCITE7* -type f >> $HOME/tmp/files.txt

#grep "md5sums" $HOME/tmp/files1.txt > $HOME/tmp/files.txt

echo "URL,FILE" > $GFS/PROJECTS/TfCf/Nisha.files.csv

while read file; do
  echo "$file"
  
  file2=$(echo $file | sed 's|\/|___|g')
  
  WT_API_PLUS_USER_TOKEN=bqxkm8o694-pcQiDBB4nYn9watVqxeqhvu5caRbNhELjjjvvUhSnH9pRdm8 $CODEBASE/wtclient//wtclient-linux-386 upload $file &> $HOME/tmp/WeTransfer/$file2.tmp
  
  url=$(grep "https" $HOME/tmp/WeTransfer/$file2.tmp | grep -v "Thank you")
  
  echo "${url},${file}" >> $GFS/PROJECTS/TfCf/Nisha.files.csv
  
done < $HOME/tmp/files.txt