sed -i 's/dnadiff-results\///g' $1
sed -i 's/_vs_/  /' $1
sed -i 's/\.report:AvgIdentity//g' $1
sed -i 's/\.report:AlignedBases//g' $1
sed -i -E 's/ +/	/g' $1
