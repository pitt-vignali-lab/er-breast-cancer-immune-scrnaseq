mkdir Outcomes
mkdir Logs

for i in $(ls -d *GEX*.gz | rev | cut -c 13- | rev | uniq)
do 
	mkdir $i 
	mv $i*.gz $i 
done

NJOBS=($(ls -d *GEX | wc -l))
sbatch --array=1-$NJOBS 2.array.sbatch