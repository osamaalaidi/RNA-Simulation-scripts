
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
	mkdir $i
	cd $i
	cp '../RiboSwitch_v15.py' .
	cp '../substructs_energies_main_plus_Aptmer.txt' .
	cp '../substructs_seqs_main_plus_Aptmer.txt' .
	cp '../Substruct_v15_mismatches_0_shortRun.py' .
	cp /media/osama/Data/Riboswitch/scriptsriboswitch/steady_state/version_13/e_range_11_310/$i/SAMI.fasta .
	cd ..
done


for i in $dirList1; do
	mkdir $i
	cd $i
	cp '../RiboSwitch_v15.py' .
	cp '../substructs_energies_main_plus_Aptmer.txt' .
	cp '../substructs_seqs_main_plus_Aptmer.txt' .
	cp '../Substruct_v15_mismatches_1_shortRun.py' .
	cp /media/osama/Data/Riboswitch/scriptsriboswitch/steady_state/version_13/e_range_11_310/$i/SAMI.fasta .
	cd ..
done
