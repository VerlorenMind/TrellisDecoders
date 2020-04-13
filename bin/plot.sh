touch output
for SNR in `seq $1 0.25 $2`
do
	echo $SNR
	./beast 1000 1000000 $SNR $3 $4 >> $4-out 
done
