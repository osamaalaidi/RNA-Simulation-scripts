
dirList0=' 
		145_complete_AT
		146_complete_AT
		147_complete_AT
		148_complete_AT
		149_complete_AT
		150_complete_AT
		151_complete_AT
		152_complete_AT '


dirList1='
		145_incomplete_AT
		146_incomplete_AT
		147_incomplete_AT
		148_incomplete_AT
		149_incomplete_AT
		150_incomplete_AT
		151_incomplete_AT
		152_incomplete_AT '

for i in $dirList0; do
	cd $i
	echo 'Getting suboptimal structures ..'
	python RiboSwitch_v15.py
	echo 'Starting.. ' $i
	python Substruct_v15_mismatches_0_shortRun.py
	echo $i, ' Complete ..'
	cd ..
done


for i in $dirList1; do
	cd $i
	echo 'Getting suboptimal structures ..'
	python RiboSwitch_v15.py
	echo 'Starting.. ' $i
	python Substruct_v15_mismatches_1_shortRun.py
	echo $i, ' Complete ..'
	cd ..
done
