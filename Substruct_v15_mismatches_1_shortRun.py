from RiboSwitch_v15 import Vienna
import sys
import os
import matplotlib.pyplot as mplt
import cdecimal
import numpy as np
import re
import time
time_start = time.time()

# in version 13:
	# Aptmer AT is not 5AT
	# Plots were made to be 100 points only
	# an alternative AT search option was added e.g. ))).)) or )))).)

# Setting numbers percisions upto 100 decimal places
myContext = cdecimal.Context(prec=100, rounding='ROUND_HALF_DOWN')
cdecimal.setcontext(myContext)
D= cdecimal.Decimal


def read_tables():
	substructs_energies_file= open('substructs_energies_main_plus_Aptmer.txt', 'r')
	substructs_energies		= substructs_energies_file.readlines()
	substructs_energies_file.close()

	energies_header_list	=	[]
	construct_B_over_A_list	=	[]
	std_dev_B_over_A_list	=	[]
	complex_conc_AB_list	=	[]
	RNA_conc_B_list			=	[]
	SAM_conc_A_list			=	[]
	ka_list					=	[]
	kd_list					=	[]
	dG0_ka_list				=	[]
	dG0_kd_list				=	[]
	dG0_non_specific_list	=	[]
	dG0_specific_list		=	[]


	for i in range(len(substructs_energies)):
		substruct_energies= substructs_energies[i].strip('\n')
		if i==0:
			pass
		else:
			name,construct_B_over_A,std_dev_B_over_A,complex_conc_AB,RNA_conc_B,SAM_conc_A,ka,kd,dG0_ka,dG0_kd,dG0_non_specific,dG0_specific	=	substruct_energies.split('\t')

			energies_header_list.append(str(name))
			construct_B_over_A_list.append(float(construct_B_over_A))
			std_dev_B_over_A_list.append(float(std_dev_B_over_A))
			complex_conc_AB_list.append(float(complex_conc_AB))
			RNA_conc_B_list.append(float(RNA_conc_B))
			SAM_conc_A_list.append(float(SAM_conc_A))
			ka_list.append(float(ka))
			kd_list.append(float(kd))
			dG0_ka_list.append(float(dG0_ka))
			dG0_kd_list.append(float(dG0_kd))
			dG0_non_specific_list.append(float(dG0_non_specific))
			dG0_specific_list.append(float(dG0_specific))


	substructs_seqs_file= open('substructs_seqs_main_plus_Aptmer.txt', 'r')
	substructs_seqs		= substructs_seqs_file.readlines()
	substructs_seqs_file.close()

	seqs_header_list		=	[]
	construct_list			=	[]
	P1_helix_indix_1_list	=	[]
	P1_helix_indix_2_list	=	[]
	P1_sequence_list		=	[]
	P1_structure_list		=	[]
	AT_helix_indix_1_list	=	[]
	AT_helix_indix_2_list	=	[]
	AT_sequence_list		=	[]
	AT_structure_list		=	[]
	Mg_list					=	[]
	AT_structure_list_slip	=	[]


	for i in range(len(substructs_seqs)):
		substruct_seqs= substructs_seqs[i].strip('\n')
		if i==0:
			pass
		else:
			construct,P1_helix_indix_1,P1_helix_indix_2,P1_sequence,P1_structure,AT_helix_indix_1,AT_helix_indix_2,AT_sequence,AT_structure,Mg,AT_structure_slip	=	substruct_seqs.split('\t')
		
			construct_list.append(str(construct))
			P1_helix_indix_1_list.append(str(P1_helix_indix_1))
			P1_helix_indix_2_list.append(str(P1_helix_indix_2))
			P1_sequence_list.append(str(P1_sequence))
			P1_structure_list.append(str(P1_structure))
			AT_helix_indix_1_list.append(str(AT_helix_indix_1))
			AT_helix_indix_2_list.append(str(AT_helix_indix_2))
			AT_sequence_list.append(str(AT_sequence))
			AT_structure_list.append(str(AT_structure))
			Mg_list.append(str(Mg))
			AT_structure_list_slip.append(str(AT_structure_slip))
	
	# Generate a dictionary for all data
	substrucs_data	= {
						'energies_header_list':energies_header_list,
						'construct_B_over_A_list':construct_B_over_A_list,
						'std_dev_B_over_A_list':std_dev_B_over_A_list,
						'complex_conc_AB_list':complex_conc_AB_list,
						'RNA_conc_B_list':RNA_conc_B_list,
						'SAM_conc_A_list':SAM_conc_A_list,
						'ka_list':ka_list,
						'kd_list':kd_list,	
						'dG0_ka_list':dG0_ka_list,
						'dG0_kd_list':dG0_kd_list,
						'dG0_non_specific_list':dG0_non_specific_list,
						'dG0_specific_list':dG0_specific_list,

						'seqs_header_list':seqs_header_list,
						'construct_list':construct_list,
						'P1_helix_indix_1_list':P1_helix_indix_1_list,
						'P1_helix_indix_2_list':P1_helix_indix_2_list,
						'P1_sequence_list':P1_sequence_list,
						'P1_structure_list':P1_structure_list,
						'AT_helix_indix_1_list':AT_helix_indix_1_list,
						'AT_helix_indix_2_list':AT_helix_indix_2_list,
						'AT_sequence_list':AT_sequence_list,
						'AT_structure_list':AT_structure_list,
						'Mg_list':Mg_list,
						'AT_structure_list_slip': AT_structure_list_slip
					}
	return substrucs_data


def multi_mplot(x, y, label_x= 'x', label_y='y', outFile='FreeEnergy_versus_P.png', col1='r.', col2='b-'):
		'''A plotting fuction for suing matPlotLib'''
		print 'Creating plots of free energies and probabilities ..'
		mplt.plot(x,y, col1)
		mplt.plot(x,y, col2)		
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.savefig(outFile)
		#mplt.cla()
		#mplt.clf()		
		return

def mplot(x, y, label_x= 'x', label_y='y', outFile='FreeEnergy_versus_P.png'):
		'''A plotting fuction for suing matPlotLib'''
		print 'Creating plots of free energies and probabilities ..'
		mplt.plot(x,y, 'r.')
		#mplt.plot(x,y, 'b-')		
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.savefig(outFile)
		mplt.cla()
		mplt.clf()		
		return


if __name__=='__main__':
	commands			=	sys.argv[1:]
	main_dir			=	os.getcwd()
	mainPath			=	os.getcwd()+'/'#'/' instead for linux					# *** changed in linux version
	Vienna_path			=	'/home/osama/Programs/viennarna_2.2.7-1_amd64/usr/bin/'	# *** changed in linux version
	RNAFold_path		=	Vienna_path+'RNAFold'									# *** changed in linux version
	RNASubOpt_path		=	Vienna_path+'RNAsubopt'									# *** changed in linux version
	RNAplot_path		=	Vienna_path+'RNAplot'									# *** changed in linux version

	#seqFilePath			=	'2gis.fasta'
	seqFilePath			=	'SAMI.fasta'
	logFileName			=	'report.log'
	logFile				= 	open(logFileName,'w')
	mismatches			= 1

	################################################################################

	# Get dictionary of sub-structures and thier binding energies
	substrucs_data		= read_tables()

	# Read subopt output
	outPath				= mainPath+'RNASubOpt.out'
	p					= Vienna(mainPath, RNAFold_path, RNASubOpt_path, RNAplot_path, seqFilePath, commands)
	strucs, freeEnergies= p.readOutFile(outPath)
	seq_id, seq			= p.readFastaSeqs(seqFilePath)

	# Partition function
	Q					=	p.PartitionFunction(freeEnergies)
	# Porbabilities
	probArray			=	p.Probability(freeEnergies, Q, fileName='probabilities_before_binding')


	substrucs_Names			= substrucs_data['construct_list']
	binding_energies		= substrucs_data['dG0_specific_list']
	P1_start_end_index_1	= substrucs_data['P1_helix_indix_1_list']
	P1_start_end_index_2	= substrucs_data['P1_helix_indix_2_list']	
	P1_substrucs			= substrucs_data['P1_structure_list']
	P1_sequence_list		= substrucs_data['P1_sequence_list']

	AT_start_end_index_1	= substrucs_data['AT_helix_indix_1_list']
	AT_start_end_index_2	= substrucs_data['AT_helix_indix_2_list']
	AT_substrucs			= substrucs_data['AT_structure_list']
	AT_substrucs_slip		= substrucs_data['AT_structure_list_slip']
	AT_sequence_list		= substrucs_data['AT_sequence_list']

	all_indecies=[]
	all_Qs=[]
	all_probArrays=[]
	all_freeEnergies=[]
	all_diff=[]
	print 'Number of substructures scanned: ',len(P1_substrucs)

	cwd=os.getcwd()

	## Writing table - input simulation		
	tableFile =	open('substructs_energies_results.txt','w')
	print >> tableFile, 'Construct','\t','Mean B/A','\t','Standard deviation B/A','\t','Complex conc. [AB]','\t','RNA conc. [B]','\t','SAM conc. [A]','\t','Ka= [AB]/[A]x[B]','\t','Kd = 1/Ka','\t','dG0= -RTxln(Ka)','\t','dG0= RTxln(Kd)','\t','dG0 (non-specific)','\t','dG0 (specific)','\t','P','\t','Q','\t','P_afterBinding','\t','Q_afterBindng'


	# Search for sub-structures through all suboptimal structures
	for i in range(len(substrucs_Names)):

		os.mkdir(substrucs_Names[i])
		os.chdir(substrucs_Names[i])
		substrucs_Name =substrucs_Names[i]
		#print substrucs_data
		logFileName			= 'sub_struct_'+str(i)+'.log'
		logFile				= open(logFileName,'w')
		
		logFileName2			= 'sub_struct_2_'+str(i)+'.log'
		logFile2				= open(logFileName2,'w')

		print >> logFile, 'The Search for substructures'
		#print >> logFile, 'Molecule Name: ',substrucs_Names[i]
		print >> logFile, 'Binding Energy: ', binding_energies[i]
		print >> logFile, 'Substructure:', P1_substrucs[i], AT_substrucs[i]
		
		# P1 search
		P1_substrucs[i]					= P1_substrucs[i].strip(' ')
		P1_1_start, P1_2_end 			= P1_substrucs[i].split('NNN')
		P1_start_index_1, P1_end_index_1 = P1_start_end_index_1[i].split(':')
		P1_start_index_2, P1_end_index_2 = P1_start_end_index_2[i].split(':')
#		P1_pattern = re.compile(P1_1_start+'.*'+P1_2_end)
		## need to add a check for the sequence
		P1_seq_1,P1_seq_2  = P1_sequence_list[i].split('NNN')
		P1_seq_pattern = re.compile(P1_seq_1+'.*'+P1_seq_2)
		## need to add the check for bracket blanace

		# AT search
		AT_substrucs[i]						= AT_substrucs[i].strip(' ')
		AT_1_start, AT_2_end				= AT_substrucs[i].split('AAAAUCACUGACAAA')
		AT_start_index_1, AT_end_index_1	= AT_start_end_index_1[i].split(':')
		AT_start_index_2, AT_end_index_2	= AT_start_end_index_2[i].split(':')

		# Accounting for slippage
		AT_1_start_slip, AT_2_end_slip		= AT_substrucs_slip[i].split('AAAAUCACUGACAAA')

#		AT_pattern = re.compile(AT_1_start+r'...............'+AT_2_end)
		## need to add a check for the sequence
		AT_seq_pattern = re.compile(AT_sequence_list[i])

#		P1_seq_index	= P1_seq_pattern.match(seq, int(P1_start_index_1)-1, int(P1_end_index_2) )
#		AT_seq_index	= AT_seq_pattern.match(seq, int(AT_start_index_1)-1, int(AT_end_index_2) )

		#P1_start_index_1, P1_end_index_1 = (8-P1_end_index_1)+P1_start_index_1, (8-P1_end_index_1)+P1_end_index_1 # To make all construct same start
		P1_seq_index	= P1_seq_pattern.match(seq[8-int(P1_end_index_1): int(P1_end_index_2) ])	# ok
		AT_seq_index	= AT_seq_pattern.match(seq[ int(AT_start_index_1)-1: int(AT_end_index_2)] )# ok

		print >> logFile, 'Sequence:'
		print >> logFile, seq[8-int(P1_end_index_1): 8], ' , ', seq[int(AT_start_index_1)-1: int(AT_end_index_2)] # ok

		print 'P1 helix: ', P1_substrucs[i], P1_start_index_1, P1_end_index_1, P1_start_index_2, P1_end_index_2
		print ' AT: ', AT_substrucs[i], AT_start_index_1, AT_end_index_1, AT_start_index_2, AT_end_index_2

		ind_file= open('indicies', 'w')
		print >> ind_file, 's', 'j', 'P1_index_1', 'P1_index_2', 'P1_seq_index', 'AT_index_1', 'AT_index_2', 'AT_seq_index', 'P1_count_index', 'AT_count_index', 'AT_index_2_slip'

		indecies = []

		for j in range(len(strucs)):

			P1_index_1=		P1_1_start	== strucs[j][8-int(P1_end_index_1): 8-int(P1_end_index_1)+len(P1_1_start) ]
			P1_index_2=		P1_2_end	== strucs[j][int(P1_start_index_2)-1: int(P1_end_index_2)]

			# Bracket count condition
			P1_count_index =strucs[j][8-int(P1_end_index_1): int(P1_end_index_2)].count('(') == strucs[j][8-int(P1_end_index_1): int(P1_end_index_2)].count(')')

#			print 'P1',P1_index_1, P1_1_start,'= ?', strucs[j][8-int(P1_end_index_1): 8-int(P1_end_index_1)+len(P1_1_start) ]
#			print 'P2',P1_index_2, P1_2_end,'= ?', strucs[j][int(P1_start_index_2)-1: int(P1_end_index_2)]
			'''
			if 'Aptamer' in substrucs_Name:
				AT_index_1 		=	str(AT_1_start)		!= strucs[j][ int(AT_start_index_1)-1: int(AT_end_index_1)]
				AT_index_2 		=	str(AT_2_end)		!= strucs[j][ int(AT_start_index_2)-1: int(AT_end_index_2)]
				AT_index_2_slip =	str(AT_2_end_slip)	!= strucs[j][ int(AT_start_index_2)-1: int(AT_end_index_2)]			# AT slippage
			else:
				AT_index_1 		=	str(AT_1_start)		== strucs[j][ int(AT_start_index_1)-1: int(AT_end_index_1)]
				AT_index_2 		=	str(AT_2_end)		== strucs[j][ int(AT_start_index_2)-1: int(AT_end_index_2)]
				AT_index_2_slip =	str(AT_2_end_slip)	== strucs[j][ int(AT_start_index_2)-1: int(AT_end_index_2)]			# AT slippage

			'''

#			print AT_index
#			print str(AT_1_start+'...............'+AT_2_end), '  |  ', strucs[j][ int(AT_start_index_1)-1: int(AT_end_index_2)]
#			AT_index= re.match(AT_1_start+r'...............'+AT_2_end, strucs[j][ int(AT_start_index_1)-1: int(AT_end_index_2)] )

			# Bracket count conditions
			AT_count_index=strucs[j][ int(AT_start_index_1)-1: int(AT_end_index_2)].count('(') == strucs[j][ int(AT_start_index_1)-1: int(AT_end_index_2)].count(')')

			if AT_count_index ==True:
				#N_bp = strucs[j][ int(AT_start_index_1)-1: int(AT_end_index_2)].count('(')
				N_bp = strucs[j][ int(AT_start_index_1)-1: int(AT_end_index_1)].count('(')
				if 'Aptamer' in substrucs_Name:
					AT_len = 4
					if N_bp < AT_len:
						AT_index_1 = True
						AT_index_2 = True
						AT_index_2_slip = True
					else:
						AT_index_1 = False
						AT_index_2 = False
						AT_index_2_slip = False

				else:
					AT_len= int( substrucs_Name[4:].strip('AT') )
					if (N_bp == AT_len) or (N_bp == (AT_len - mismatches)):
						AT_index_1 = True
						AT_index_2 = True
						AT_index_2_slip = True
					else:
						AT_index_1 = False
						AT_index_2 = False
						AT_index_2_slip = False
			else:
				AT_index_1 = False
				AT_index_2 = False
				AT_index_2_slip = False
			

			'''
			if P1_index_1 and P1_index_2:
				print >> logFile, 'P1_index is True'
			else:
				print >> logFile, 'P1_index is False'			

			if P1_seq_index:
				print >> logFile, 'P1_seq_index is True'
			else:
				print >> logFile, 'P1_seq_index is False'
			
			if AT_index_1:
				print >> logFile, 'AT_index_1 is True'
			else:
				print >> logFile, 'AT_index_1 is False'

			if AT_index_2:
				print >> logFile, 'AT_index_2 is True'
			else:
				print >> logFile, 'AT_index_2 is False'

			if AT_seq_index:
				print >> logFile, 'AT_seq_index is True'
			else:
				print >> logFile, 'AT_seq_index is False'
			'''

		#	if P1_index_1 and P1_index_2 and P1_seq_index and AT_index and AT_seq_index and P1_count_index and AT_count_index: #changed for mutants
#			if AT_index and AT_count_index:
			if P1_index_1 and P1_index_2 and AT_index_1 and AT_index_2 and P1_count_index and AT_count_index and P1_seq_index and AT_seq_index:
				s=True
			elif P1_index_1 and P1_index_2 and AT_index_1 and AT_index_2_slip and P1_count_index and AT_count_index and P1_seq_index and AT_seq_index:
				s=True # slipage
			else:
				s=False
			indecies.append(s)

			print >> logFile2,'Structure:', j, strucs[j] 
			print >> logFile2, j, '\t', strucs[j][8-int(P1_end_index_1):int(P1_end_index_2)],' , ', strucs[j][int(AT_start_index_1)-1: int(AT_end_index_2)], s 

			print >> ind_file, s, j, P1_index_1, P1_index_2, P1_seq_index, AT_index_1, AT_index_2, AT_seq_index, P1_count_index, AT_count_index, AT_index_2_slip
#			print >> ind_file, s, j, AT_index, AT_seq_index, AT_count_index

			### print j, s, int(P1_start_index_1), int(P1_end_index_1), int(P1_start_index_2), int(P1_end_index_2), int(AT_start_index_1), int(AT_end_index_1), int(AT_start_index_2), int(AT_end_index_2)
			#print P1_index.string[P1_index.start(P1_index.group(0)):P1_index.end(P1_index.group(0))], P1_pattern
			#print AT_index.string[AT_index.start(AT_index.group(0)):AT_index.end(AT_index.group(0))], AT_pattern
			#indecies					= p.subStructure_search(strucs, substrucs[i], P1_start_index[i])
		logFile2.close()
		indecies = np.array(indecies, dtype='bool')
		FreeEnergies_with_Substrucst= np.where( indecies, freeEnergies+D(binding_energies[i]), freeEnergies  )
		# Partition function
		Q_with_Substrucst			= p.PartitionFunction(FreeEnergies_with_Substrucst)
		# Porbabilities
		probArray_with_Substrucst	= p.Probability(FreeEnergies_with_Substrucst, Q_with_Substrucst, fileName='probabilities_substruc_'+str(i))
		print >> logFile, '\n'
		print >> logFile, 'Struc No.','\t','Presence of substructure(0 if not found)','\t','Free_Energy_before_binding','\t','Probability_before_binding','\t','Free_Energy_after_binding','\t','Probability_after_binding'#, 'Probability it does not contain the sub-structure' 
		for j in range(len(strucs)):
			print >> logFile, j,'\t', indecies[j],'\t','%G'%freeEnergies[j],'\t','%G'%probArray[j],'\t',FreeEnergies_with_Substrucst[j],'\t','%G'%probArray_with_Substrucst[j]#, '%G'%P_with_Substrucst[i]#, '%G'%P_without_Substrucst[i]

		# Plots
		mplot(freeEnergies, probArray, label_x='Free Energy', label_y='Probability', outFile= 'FreeEnergies_versus_P1.png')
		mplot(FreeEnergies_with_Substrucst, probArray_with_Substrucst, label_x='Free Energy', label_y='Probability', outFile= 'FreeEnergies_with_Substrucst_versus_P1.png')

		multi_mplot(freeEnergies, probArray, label_x='Free Energy', label_y='Probability', outFile= 'FreeEnergies_versus_P2.png', col1='.g',col2='-m')
		mplt.cla()
		mplt.clf()
		multi_mplot(FreeEnergies_with_Substrucst, probArray_with_Substrucst, label_x='Free Energy', label_y='Probability', outFile= 'FreeEnergies_with_Substrucst_versus_P2.png', col1='.g',col2='-m')
		mplt.cla()
		mplt.clf()

#		FreeEnergies_with_Substrucst_only= np.where( indecies, FreeEnergies_with_Substrucst, np.empty_like(FreeEnergies_with_Substrucst)  )
#		probArray_with_Substrucst_only	=  np.where( indecies, probArray_with_Substrucst, np.empty_like(probArray_with_Substrucst)  )
#		mplot(FreeEnergies_with_Substrucst_only, probArray_with_Substrucst_only, label_x='Free Energy', label_y='Probability', outFile= 'FreeEnergies_versus_P_withsubstruct_only.png')

		logFile.close()
		logFile2=open('Summary.txt', 'w')
		all_indecies.append(indecies)
		all_Qs.append(Q_with_Substrucst)
		all_probArrays.append(probArray_with_Substrucst)
		all_freeEnergies.append(FreeEnergies_with_Substrucst)

		P_with_Substrucst_before		= np.where( indecies, probArray, np.zeros_like(probArray) )
		P_without_Substrucst_before		= np.where( indecies, np.zeros_like(probArray), probArray )
		P_Total_before					= np.sum(P_with_Substrucst_before) + np.sum(P_without_Substrucst_before)
		mplot(FreeEnergies_with_Substrucst, P_with_Substrucst_before, label_x='Free Energy', label_y='Probability', outFile= 'FreeEnergies_versus_P_withsubstruct_before_binding.png')

		print >> logFile2, 'Sum of Probabilities (before binding) of structures with the substructure =', '%E'%np.sum(P_with_Substrucst_before)
		print >> logFile2, 'Sum of Probabilities (before binding) of structures without the substructure =', '%E'%np.sum(P_without_Substrucst_before)
		print >> logFile2, 'Sum of all Probabilities (before binding) =', '%E'%P_Total_before
		print >> logFile2, 'Q_before=', '%E'%Q

		P_with_Substrucst_after		= np.where( indecies, probArray_with_Substrucst, np.zeros_like(probArray) )
		P_without_Substrucst_after	= np.where( indecies, np.zeros_like(probArray), probArray_with_Substrucst )
		P_Total_after				= np.sum(P_with_Substrucst_after) + np.sum(P_without_Substrucst_after)
		diff						= np.sum(P_with_Substrucst_after) - np.sum(P_with_Substrucst_before)

		mplot(FreeEnergies_with_Substrucst, P_with_Substrucst_after, label_x='Free Energy', label_y='Probability', outFile= 'FreeEnergies_versus_P_withsubstruct_after_binding.png')

		print >> logFile2, 'Sum of Probabilities (after binding) of structures with the substructure =', '%E'%np.sum(P_with_Substrucst_after)
		print >> logFile2, 'Sum of Probabilities (after binding) of structures without the substructure =', '%E'%np.sum(P_without_Substrucst_after)
		print >> logFile2, 'Sum of all Probabilities (after binding) =', '%E'%P_Total_after
		print >> logFile2, 'Q_after=', '%E'%Q_with_Substrucst
	
		print >> logFile2, 'Difference in probability between before and after dP=', '%E'%diff
		
		try:
			print >> logFile2, 'P_after /P_before=', '%E'%(np.sum(P_with_Substrucst_after)/np.sum(P_with_Substrucst_before) )
		except ZeroDivisionError:
			print >> logFile2, '0.0'

	#	try:
	#		print >> logFile2, '(DP/P_before)*100=', '%E'%( np.array(diff*D(100))/np.sum(P_with_Substrucst_before) )
	#	except ZeroDivisionError:
	#		print >> logFile2, '0.0'

		all_diff.append(diff)
		logFile2.close()
		os.chdir(cwd)

#		break
	all_indecies	= np.array(all_indecies)
	all_Qs			= np.array(all_Qs)
	all_probArrays	= np.array(all_probArrays)
	all_freeEnergies= np.array(all_freeEnergies)
	all_diff		= np.array(all_diff)
	
	'''
	# Needs iteration
	FreeEnergies_with_Substrucst_all= np.array(all_freeEnergies[0])
	for i in range(len(all_indecies)-1):
		a=FreeEnergies_with_Substrucst_all
		FreeEnergies_with_Substrucst_all	= np.where( all_indecies[i+1], all_freeEnergies[i+1], a )#all_freeEnergies[i] )
	'''
	# Getting final energies and probabilities

	FreeEnergies_with_Substrucst_all= freeEnergies
	for i in range(len(all_indecies)):
		a=FreeEnergies_with_Substrucst_all
		FreeEnergies_with_Substrucst_all	= np.where( all_indecies[i], all_freeEnergies[i], a )#all_freeEnergies[i] )


		###	#FreeEnergies_with_Substrucst_zeros=	np.where( all_indecies==1, freeEnergies, np.zeros_like(probArray)  )  #??
	# Partition function for all substructures
	Q_all			=	p.PartitionFunction(FreeEnergies_with_Substrucst_all)
	# Porbabilities for all substructures
	probArray_all	=	p.Probability(FreeEnergies_with_Substrucst_all, Q_all, fileName='probabilities_after_binding_final')

	print 'probArray_all sum = ', np.sum(probArray_all)
	print 'System Q = ', Q_all

	Total_P_with_Substrucst_before		= D(0.0)
	Total_P_with_Substrucst_after		= D(0.0)
	Total_P_without_Substrucst_before	= D(0.0)
	Total_P_without_Substrucst_after	= D(0.0)
	indicies_without_Substrucst_before_final= D(0.0)
	indicies_without_Substrucst_after_final= D(0.0)
#P_with_Substrucst_before_final, P_without_Substrucst_before_final
	for i in range(len(substrucs_Names)):

		P_with_Substrucst_before_final		= np.where( all_indecies[i], probArray, np.zeros_like(probArray) )
		P_with_Substrucst_after_final		= np.where( all_indecies[i], probArray_all, np.zeros_like(probArray_all) )

		indicies_without_Substrucst_before_final= all_indecies[i] + indicies_without_Substrucst_before_final#np.where( all_indecies[i], np.ones_like(probArray), np.zeros_like(probArray) )
		indicies_without_Substrucst_after_final	= all_indecies[i] + indicies_without_Substrucst_after_final#np.where( all_indecies[i], np.ones_like(probArray_all), np.zeros_like(probArray_all) )

		print 'substruct',str(i+1), substrucs_data['energies_header_list'][i]
		print >> tableFile, substrucs_data['energies_header_list'][i],'\t',\
							substrucs_data['construct_B_over_A_list'][i],'\t',\
							substrucs_data['std_dev_B_over_A_list'][i],'\t',\
							substrucs_data['complex_conc_AB_list'][i],'\t',\
							substrucs_data['RNA_conc_B_list'][i],'\t',\
							substrucs_data['SAM_conc_A_list'][i],'\t',\
							substrucs_data['ka_list'][i],'\t',\
							substrucs_data['kd_list'][i],'\t',\
							substrucs_data['dG0_ka_list'][i],'\t',\
							substrucs_data['dG0_kd_list'][i],'\t',\
							substrucs_data['dG0_non_specific_list'][i],'\t',\
							substrucs_data['dG0_specific_list'][i],'\t',\
							'%E'%np.sum(P_with_Substrucst_before_final),'\t',\
							'%E'%Q,'\t',\
							'%E'%np.sum(P_with_Substrucst_after_final),'\t',\
							'%E'%Q_all

		Total_P_with_Substrucst_before	= Total_P_with_Substrucst_before 	+ np.sum(P_with_Substrucst_before_final)
		Total_P_with_Substrucst_after	= Total_P_with_Substrucst_after	 	+ np.sum(P_with_Substrucst_after_final)

#		print 'substruct', i, '= ', np.sum(P_without_Substrucst_before_final), P_without_Substrucst_before_final 
#		print 'substruct', i, '= ', np.sum(P_without_Substrucst_after_final), P_without_Substrucst_after_final

	tableFile.close()

	P_without_Substrucst_before_final	= np.where( indicies_without_Substrucst_before_final, np.zeros_like(probArray), probArray )
	P_without_Substrucst_after_final	= np.where( indicies_without_Substrucst_after_final, np.zeros_like(probArray_all), probArray_all )

	Total_P_without_Substrucst_before	= Total_P_without_Substrucst_before + np.sum(P_without_Substrucst_before_final)
	Total_P_without_Substrucst_after	= Total_P_without_Substrucst_after	+ np.sum(P_without_Substrucst_after_final)

	Sum_of_all_system_P_before			= Total_P_with_Substrucst_before	+ Total_P_without_Substrucst_before
	Sum_of_all_system_P_after			= Total_P_with_Substrucst_after		+ Total_P_without_Substrucst_after


	f=open('all_substruct_summary', 'w')
	print >> f, 'Total_P_with_Substrucst_before =', Total_P_with_Substrucst_before
	print >> f, 'Total_P_with_Substrucst_after =', Total_P_with_Substrucst_after
	print >> f, 'Total_P_without_Substrucst_before =', Total_P_without_Substrucst_before
	print >> f, 'Total_P_without_Substrucst_after =', Total_P_without_Substrucst_after
	print >> f, 'Sum_of_all_system_P_before =', Sum_of_all_system_P_before
	print >> f, 'Sum_of_all_system_P_after =', Sum_of_all_system_P_after
	f.close()

	f=open('all_substruct_data', 'w')
	print >>f, 'Sum of all probablity differences between before and after binding = ', np.sum(all_diff)
	print 'all_indecies=', all_indecies
	print >> f,'Partition function =', Q_all
	print >> f,'No.','\t','Free_Energy_after_all_Substructures','\t','Probability_after_all_substructures'
	for i in range(len(all_indecies[0])):
		print >> f,i,'\t', FreeEnergies_with_Substrucst_all[i],'\t',probArray_all[i]
	f.close()

	##########################################################################
	time_end = time.time()
	print '\n'
	print 'Time of excultion: ', time_end - time_start, 'seconds'

