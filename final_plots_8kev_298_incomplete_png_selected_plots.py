import matplotlib.pyplot as mplt
import cdecimal
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

# Setting numbers percisions upto 100 decima places
myContext = cdecimal.Context(prec=100, rounding='ROUND_HALF_DOWN')
cdecimal.setcontext(myContext)
D= cdecimal.Decimal
ext='_new_selected.png'

# 26,24 or 7.5,7 or for NAR: 9,7 (full page) => size_a,size_b
# 19,12 or 7,4.4 or for NAR: 5.2,3.3 (vertical) or 3.3,2.08 (horizontal) => size_c,size_d

size_a,size_b= 7.5,7#26,24
size_c,size_d= 3.3,2.08#19,12
font_size_ratio=28/7.0 # font/A4 page width
subplotHPad=0.5
subplotVPad=0.5
figPad=0.25

font_size= font_size_ratio * size_b
markeredgewidth_ratio=1.0/12.0

mpl.rcParams['font.weight']='bold'
mpl.rcParams['axes.linewidth']='3'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['grid.linestyle']= '-'
mpl.rcParams['font.size']=font_size
mpl.rcParams['savefig.bbox']='tight' # standard default in tiff
mpl.rcParams['savefig.pad_inches']=figPad # not working in tiff


def set_globVar():
	global globVar 
	globVar= 0
	return


def multi_mplot(x, y, label_x= 'x', label_y='y', outFile='FreeEnergy_versus_P'+ext, Title='', col1='r.', col2='b.', legends=None, legend_loc=2, N=111):
		'''A plotting fuction for using matPlotLib'''
		print 'Creating plots of free energies and probabilities ..'
		#fig 	= mplt.figure()
		#ax 	= fig.add_subplot(N)
		if type(x[0])==str:
		#	ax.set_xticklabels(x, rotation=45)
			x= np.array(range(len(x)))
			mplt.xticks(x)
			#mplt.xticks(range(len(substrucs_names)),substrucs_names,)
		mplt.plot(x,y, col1, linewidth=6, markersize=font_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
		#mplt.plot(x,y, col2)		
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.title(Title)
		mplt.grid(True)
		#mplt.tight_layout(pad=subplotPad)
		if legends != None:
			#mplt.legend(legends)
			mplt.legend(legends, bbox_to_anchor=(1.01, 1), shadow=True, fancybox=True, loc=legend_loc, borderaxespad=0)
		mplt.savefig(outFile, dpi=600)
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)
		#mplt.cla()
		#mplt.clf()
		return


def multi_subplot(x, y, label_x= 'x', label_y='y', outFile='FreeEnergy_versus_P'+ext, Title='', col1='r.', col2='b.', legends=None, legend_loc=2, N=111, rows=1, col=1, n=1, ni=0):
		'''A plotting fuction for using matPlotLib'''
		print 'Creating plots of free energies and probabilities ..'

		rows	= [1,2,3,1,2,3,1,2,3]
		col		= [1,1,1,2,2,2,3,3,3]
		Subplot_N	= int(str(rows[n])+str(col[n])+str(n+1))
		print Subplot_N
		if n==0 and ni==0:
			fig 	= mplt.figure()
			ax 		= fig.add_subplot(3, 3, n+1)
		else:
			fig 	= mplt.gcf()
			ax 		= fig.add_subplot(3, 3, n+1)

		fig.subplots_adjust(hspace=subplotHPad, wspace=subplotVPad)
		if type(x[0])==str:
			ax.set_xticklabels(x, rotation=45)
			x= np.array(range(len(x)))
			#mplt.xticks(x)
			#mplt.xticks(range(len(substrucs_names)),substrucs_names,)
		mplt.plot(x, y, col1, linewidth=6, markersize=font_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
		#mplt.plot(x,y, col2)		
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.grid(True)
		#mplt.tight_layout(pad=subplotPad)
		if legends != None:
			#mplt.legend(legends)
			#mplt.legend(legends, bbox_to_anchor=(1.01, 1), shadow=True, fancybox=True, loc=legend_loc, borderaxespad=0)
			mplt.title('TL = '+legends[n]+' nucleotides')
		#mplt.savefig(outFile, dpi=600)
		#f=mplt.gcf()
		#f.set_size_inches(14,12)

		#mplt.cla()
		#mplt.clf()
		return


def subplot_hist(x, outfile='hist'+ext, label_x= 'x', label_y='y', Title= '', histtype='stepfilled', bins=500,legends=None, n=1, facecolor='b'):
		''' Histograms Subplots '''
		rows	= [1,2,3,1,2,3,1,2,3]
		col		= [1,1,1,2,2,2,3,3,3]
		Subplot_N	= int(str(rows[n])+str(col[n])+str(n+1))
		print Subplot_N
		if n==0:
			fig 	= mplt.figure()
			ax 		= fig.add_subplot(3, 3, n+1)

		else:
			fig 	= mplt.gcf()
			ax 		= fig.add_subplot(3, 3, n+1)
		fig.subplots_adjust(hspace=subplotHPad, wspace=subplotVPad)
		mplt.hist(x, bins=bins, histtype=histtype, facecolor=facecolor)
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.title(Title)
		mplt.grid(True)
		#mplt.tight_layout(pad=subplotPad)
		if legends != None:
			mplt.title('TL = '+legends[n]+' nucleotides')
		return



def multi_mplot_str(x, y, label_x= 'x', label_y='y', outFile='FreeEnergy_versus_P'+ext, Title='', col1='r.', col2='b.', legends=None, legend_loc=2, N=111, counter=0):
		'''A plotting fuction for using matPlotLib'''
		print 'Creating plots of free energies and probabilities ..'
		if counter==0:
			fig 	= mplt.figure()
			ax 		= fig.add_subplot(N)
			fig.subplots_adjust(hspace=subplotHPad, wspace=subplotVPad)
			ax.set_xticklabels(x, rotation=45)
		if type(x[0])==str:
			x= np.array(range(len(x)))
			#mplt.xticks(x)
			#mplt.xticks(range(len(substrucs_names)),substrucs_names,)
		mplt.plot(x,y, col1, linewidth=6, markersize=font_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
		#mplt.plot(x,y, col2)		
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.title(Title)
		mplt.grid(True)
		#mplt.tight_layout(pad=subplotPad)
		if legends != None:
			#mplt.legend(legends)
			mplt.legend(legends, bbox_to_anchor=(1.01, 1), shadow=True, fancybox=True, loc=legend_loc, borderaxespad=0)
		mplt.savefig(outFile, dpi=600)
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)
		#mplt.cla()
		#mplt.clf()
		return


def mplot(x, y, label_x= 'x', label_y='y', outFile='FreeEnergy_versus_P'+ext, Title='', legends=None, col='r.', N=111):
		'''A plotting fuction for using matPlotLib'''
		print 'Creating plots of free energies and probabilities ..'
		if type(x[0])==str:
			#mplt.xticks(x)
			fig 	= mplt.figure()
			ax 		= fig.add_subplot(N)
			ax.set_xticklabels(x, rotation=45)
			x= np.array(range(len(x)))
		#mplt.subplot(N)
		mplt.plot(x,y, col, linewidth=6, markersize=font_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
		#mplt.plot(x,y, 'b-')		
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.title(Title)
		#mplt.tight_layout(pad=subplotPad)
		if legends != None:
			#mplt.legend(legends)
			mplt.legend(legends, bbox_to_anchor=(1.02, 1), shadow=True, fancybox=True, loc=legend_loc, borderaxespad=3)
		mplt.grid(True)
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()	
		return



def mplot_err(x, y, label_x= 'x', label_y='y', outFile='FreeEnergy_versus_P'+ext, Title='', legends=None, col='r.', N=111, yerr=None):
		'''A plotting fuction for using matPlotLib'''
		print 'Creating plots of free energies and probabilities ..'
		if type(x[0])==str:
			#mplt.xticks(x)
			fig 	= mplt.figure()
			ax 		= fig.add_subplot(N)
			ax.set_xticklabels(x, rotation=45)
			x= np.array(range(len(x)))
		#mplt.subplot(N)
		print 'yerr', yerr
		#mplt.plot(x,y, col, linewidth=6, markersize=font_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
		mplt.errorbar(x,y, yerr=yerr, fmt='-bo', ecolor='r')
		#mplt.plot(x,y, 'b-')		
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.title(Title)
		#mplt.tight_layout(pad=subplotPad)
		if legends != None:
			#mplt.legend(legends)
			mplt.legend(legends, bbox_to_anchor=(1.02, 1), shadow=True, fancybox=True, loc=legend_loc, borderaxespad=3)
		mplt.grid(True)
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()	
		return


def plot_hist(x, outfile='hist'+ext, label_x= 'x', label_y='y', Title='', legends=None, col='r.', histtype='bar', bins=500 ):
	mplt.hist(x, bins=bins)#, histtype=histtype, bins=bins)
	mplt.xlabel(label_x)
	mplt.ylabel(label_y)
	mplt.title(Title)
	#mplt.tight_layout(pad=subplotPad)
	if legends != None:
		mplt.legend(legends)
	mplt.grid(True)
	f=mplt.gcf()
	f.set_size_inches(size_c,size_d)
	mplt.savefig(outfile, dpi=600)
	mplt.cla()
	mplt.clf()
	return





def read_probabilities(inFile='probabilities', outPath='./'):
	''' Reads probabilities file '''
	f=open(inFile,'r')
	lines=f.readlines()
	data=lines[1:]
	struc_ids, FreeEnergies, BolzFactors, Probabilities=[],[],[],[]
	for line in data:
		struc_id, FreeEnergy, BolzFactor, Probability= line.split()
		#struc_ids.append(int(struc_id))
		FreeEnergies.append(float(FreeEnergy))
		#BolzFactors.append(float(BolzFactor))
		Probabilities.append(float(Probability))
	#struc_ids		= np.array(struc_ids)
	FreeEnergies	= np.array(FreeEnergies)
	#BolzFactors		= np.array(BolzFactors)
	Probabilities	= np.array(Probabilities)
	# Sort
	dataArr= np.array([FreeEnergies, Probabilities])
#	dataArr= np.sort(dataArr, axis=1)
	FreeEnergies, Probabilities = dataArr[0], dataArr[1]

#	dataArr= np.sort(dataArr, axis=0)
	FreeEnergies, Probabilities = dataArr[0], dataArr[1]

	return FreeEnergies, Probabilities


def read_all_substruct_data(inFile='all_substruct_data', outPath='./'):
	''' Reads final probabilities file, i.e. after binding '''
	f=open(inFile,'r')
	lines=f.readlines()
	data=lines[3:]
	struc_ids, FreeEnergies, Probabilities=[],[],[]
	for line in data:
		struc_id, FreeEnergy, Probability= line.split()
		#struc_ids.append(int(struc_id))
		FreeEnergies.append(float(FreeEnergy))
		Probabilities.append(float(Probability))
	#struc_ids		= np.array(struc_ids)
	FreeEnergies	= np.array(FreeEnergies)
	Probabilities	= np.array(Probabilities)
	# Sort
	dataArr= np.array([FreeEnergies, Probabilities])
#	dataArr= np.sort(dataArr, axis=1)
	FreeEnergies, Probabilities = dataArr[0], dataArr[1]


#	dataArr= np.sort(dataArr, axis=0)
	FreeEnergies, Probabilities = dataArr[0], dataArr[1]

	return FreeEnergies, Probabilities



def read_tables(path):
	substructs_energies_file= open(path+'substructs_energies_results.txt', 'r')
	substructs_energies		= substructs_energies_file.readlines()
	substructs_energies_file.close()

	energies_header_list			= []
	construct_B_over_A_list			= []
	std_dev_B_over_A_list			= []
	complex_conc_AB_list			= []
	RNA_conc_B_list					= []
	SAM_conc_A_list					= []
	ka_list							= []
	kd_list							= []
	dG0_ka_list						= []
	dG0_kd_list						= []
	dG0_non_specific_list			= []
	dG0_specific_list				= []
	P_with_Substructs_before_list	= []
	Q_before_list					= []
	P_with_Substructs_after_list	= []
	Q_after_list					= []


	for i in range(len(substructs_energies)):
		substruct_energies= substructs_energies[i].strip('\n')
		if i==0:
			pass
		else:
			print substruct_energies.split('\t')
			name,construct_B_over_A,std_dev_B_over_A,complex_conc_AB,RNA_conc_B,SAM_conc_A,ka,kd,dG0_ka,dG0_kd,dG0_non_specific,dG0_specific, P_with_Substructs_before, Q_before, P_with_Substructs_after, Q_after	=	substruct_energies.split('\t')
			
			energies_header_list.append(str(name.strip()))
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
			P_with_Substructs_before_list.append(float(P_with_Substructs_before))
			Q_before_list.append(float(Q_before))
			P_with_Substructs_after_list.append(float(P_with_Substructs_after))
			Q_after_list.append(float(Q_after))



	substructs_seqs_file= open(path+'substructs_seqs_main_plus_Aptmer.txt', 'r')
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
						'P_with_Substructs_before_list':P_with_Substructs_before_list,
						'Q_before_list':Q_before_list,
						'P_with_Substructs_after_list':P_with_Substructs_after_list,
						'Q_after_list':Q_after_list,

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


def readSubstrucsData(path, subDirs, inFile='indicies'):
	substructFreq=[]
	for subDir in subDirs:
		substructIndicies=[]
		f= open(path+'/'+subDir+'/'+inFile, 'r')
		lines = f.readlines()
		#s, j, P1_index_1, P1_index_2, P1_seq_index, AT_index_1, AT_index_2, AT_seq_index, P1_count_index, AT_count_index, AT_index_2_slip = line.split('\t')
		lines= lines[1:] # to skip first line
		for line in lines:
			all_vars	= line.split(' ')
			s			= all_vars[0]
			if s=='False':
				s=0
			elif s=='True':
				s=1
			substructIndicies.append(bool(s))
		substructIndicies = np.array(substructIndicies, dtype='bool')
		substructFreq.append( np.sum(substructIndicies) )
	substructFreq = np.array(substructFreq)
	return substructFreq


def visTable(path, cwd):
	''' Prints and plots final results tables for a single length'''
	# get data for a single run/file i.e. same length, energy and temperature
	substrucs_data =	read_tables(path)

	substrucs_names		= substrucs_data['energies_header_list'][:-1]
	for i in range(len(substrucs_names)):
		substrucs_names[i] = substrucs_names[i].upper()
	substrucs_names.append('Aptamer_TS')
	return substrucs_data, substrucs_names


def sampleArr(arr, res=1000, flip=False):
	'Samples and flips an array'
	'''sampledArr =[]
	delta= int(len(arr)/res)
	for i in range(delta):
		sampledArr.append(arr[i])
	for i in range(res):
		sampledArr.append(arr[i*delta])
		print i*delta
	sampledArr = np.array(sampledArr)
	if flip==True:
		sampledArr=np.flipud(sampledArr)
	print sampledArr.size
	return sampledArr'''
	return arr



def wTitle(Title, outFileName):
	global globVar
	if globVar==0:
		f= open('plot_titles_mod', 'w')
		globVar = 1
	else:
		f= open('plot_titles_mod', 'a')
	print >> f, outFileName, ': ', Title
	f.close()
	return ''


#		color=('r','b','k','m','g','c','y')
#		color= color * plotNum


#savfig 
#pad_inches=0.07, bbox_inches='tight'
#plt
#dpi
#fig 	= plt.figure()

#		f=plt.gcf()
#		f.set_size_inches(self.plot_size_in_inch)

#subplt= plt.subplot(rows, cols, n+1)



markerfacecolor='green'


if __name__=='__main__':
	set_globVar()
	R	= 0.0019872041
	#T	= 298.15

	# for all main folders
	mainDirs=['e_range_8_298']
				#'e_range_11_310',
				#'e_range_13_295',
				#'e_range_8_295',
				#'e_range_8_310']
				#'e_range_4_295']#
	for mainDir in mainDirs:
		T	= float(mainDir[-3:]+'.15')
		cwd	=os.getcwd()
		os.chdir(cwd+'/'+mainDir)

		# for all subfolders
		cwd	=os.getcwd()
		dirs_complete = ['145_complete_AT',
						 '146_complete_AT',
						 '147_complete_AT',
						 '148_complete_AT',
						 '149_complete_AT',
						 '150_complete_AT',
						 '151_complete_AT',
						 '152_complete_AT']


		dirs_incomplete = [ '145_incomplete_AT',
							'146_incomplete_AT',
							'147_incomplete_AT',
							'148_incomplete_AT',
							'149_incomplete_AT',
							'150_incomplete_AT',
							'151_incomplete_AT',
							'152_incomplete_AT']

		trascLen				= []
		trascLen_name			= []

		P_trascLen_BB			= []
		P_trascLen_AB			= []
		E_trascLen_BB			= []
		E_trascLen_AB			= []

		P_substruct_trascLen_BB	= []
		P_substruct_trascLen_AB	= []

		P_substruct_sum_trascLen_BB	= []
		P_substruct_sum_trascLen_AB	= []


		#still to be added # accumulative increase by length or % increase
		#still  to be added # no of total structures
		#still  to be added # P of substructures


		dirs	= dirs_incomplete

		for dir_name in dirs:
			
			path	= cwd+'/'+dir_name+'/'
			
#			struc_ids_BB, FreeEnergies_BB, Probabilities_BB			= read_probabilities(inFile=path+'probabilities', outPath=path)
#			struc_ids_AB, FreeEnergies_AB, Probabilities_AB			= read_all_substruct_data(inFile=path+'all_substruct_data', outPath=path)

			substrucs_data, substrucs_names  = visTable(path, cwd)

			print 'substrucs_names =' ,substrucs_names 

			FreeEnergies_BB, Probabilities_BB			= read_probabilities(inFile=path+'probabilities', outPath=path)
			FreeEnergies_AB, Probabilities_AB			= read_all_substruct_data(inFile=path+'all_substruct_data', outPath=path)

			# frequencies
			substructFreq = readSubstrucsData(path, substrucs_data['energies_header_list'], inFile='indicies')


			# stats:
			# overall frequences per lengths
			trascLen.append(int(dir_name[:3]))
			trascLen_name.append(dir_name[:3])


			P_trascLen_BB.append(Probabilities_BB)
			P_trascLen_AB.append(Probabilities_AB)

			E_trascLen_BB.append(FreeEnergies_BB)
			E_trascLen_AB.append(FreeEnergies_AB)


			# nfold change per substruct
			# Accomulative probabilies: contribution of each substructure seperate and accomulative

		
		######## Lengths comparisons

		path	= cwd+'/'
		# stats arrays: 
		# Plot number of conformers
		color=('ro-','bo-','ko-','mo-','go-','co-','yo-', 'kx-')
		color_dot=('ro','bo','ko','mo','go','co','yo', 'kx')
		colors=('r','b','k','m','g','c','y', 'w')
		#substrucs_names			= substrucs_data['energies_header_list']
		#substrucs_names[-1]		= 'Aptamer_TS'
		trascLen_name			= trascLen_name
		trascLen				= np.array(trascLen)

		#################################################################################################### Probabilities


		########################## series plot
		# subplot both BB and AB
		outFile = path+'transLen_versus_P_BB_AB_subplot'+'_new.png'#ext
		for i in range(len(P_trascLen_BB)):
			multi_subplot(E_trascLen_BB[i], P_trascLen_BB[i], label_x='Free energy [kcal/mol]', label_y='P', outFile=outFile, Title=wTitle('P binding at increasing transcription lengths', 'transLen_versus_P_BB__AB_subplot'+ext), col1=color_dot[1], legends=trascLen_name, rows=len(P_trascLen_BB), col=1, n=i, ni=0)
			multi_subplot(E_trascLen_AB[i], P_trascLen_AB[i], label_x='Free energy [kcal/mol]', label_y='P', outFile=outFile, Title=wTitle('P binding at increasing transcription lengths', 'transLen_versus_P_BB__AB_subplot'+ext), col1=color_dot[0], legends=trascLen_name, rows=len(P_trascLen_AB), col=1, n=i, ni=1)
			mplt.legend(['BB','AB'])
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()


		########################### Histograms subplots

		# Histograms 500 bins
		outFile = path+'transLen_versus_P_BB_hist_subplot_500'+ext
		for i in range(len(P_trascLen_BB)):
			subplot_hist(P_trascLen_BB[i], outfile=path+'P_hist_before'+ext, label_x= 'Probability', label_y='Frequency', Title= wTitle('Probability distribution', 'transLen_versus_P_BB_hist_subplot'+ext), histtype='stepfilled', bins=500, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_P_AB_hist_subplot_500'+ext
		for i in range(len(P_trascLen_AB)):
			subplot_hist(P_trascLen_AB[i], outfile=path+'P_hist_after'+ext, label_x= 'Probability', label_y='Frequency', Title=wTitle('Probability distribution','transLen_versus_P_AB_hist_subplot'+ext), histtype='stepfilled', bins=500, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()


		##########################

		# Histograms 150 bins
		outFile = path+'transLen_versus_P_BB_hist_subplot_150'+ext
		for i in range(len(P_trascLen_BB)):
			subplot_hist(P_trascLen_BB[i], outfile=outFile, label_x= 'Probability', label_y='Frequency', Title= wTitle('Probability distribution', 'transLen_versus_P_BB_hist_subplot_150'+ext), histtype='stepfilled', bins=150, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_P_AB_hist_subplot_150'+ext
		for i in range(len(P_trascLen_AB)):
			subplot_hist(P_trascLen_AB[i], outfile=outFile, label_x= 'Probability', label_y='Frequency', Title=wTitle('Probability distribution','transLen_versus_P_AB_hist_subplot_150'+ext), histtype='stepfilled', bins=150, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()


		##########################

		# Histograms 50 bins
		outFile = path+'transLen_versus_P_BB_hist_subplot_50'+ext
		for i in range(len(P_trascLen_BB)):
			subplot_hist(P_trascLen_BB[i], outfile=outFile, label_x= 'Probability', label_y='Frequency', Title= wTitle('Probability distribution', 'transLen_versus_P_BB_hist_subplot_50'+ext), histtype='stepfilled', bins=50, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_P_AB_hist_subplot_50'+ext
		for i in range(len(P_trascLen_AB)):
			subplot_hist(P_trascLen_AB[i], outfile=outFile, label_x= 'Probability', label_y='Frequency', Title=wTitle('Probability distribution','transLen_versus_P_AB_hist_subplot_50'+ext), histtype='stepfilled', bins=50, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		####################################
		# Histograms 20 bins
		outFile = path+'transLen_versus_P_BB_hist_subplot_20'+ext
		for i in range(len(P_trascLen_BB)):
			subplot_hist(P_trascLen_BB[i], outfile=outFile, label_x= 'Probability', label_y='Frequency', Title= wTitle('Probability distribution', 'transLen_versus_P_BB_hist_subplot_50'+ext), histtype='stepfilled', bins=20, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_P_AB_hist_subplot_20'+ext
		for i in range(len(P_trascLen_AB)):
			subplot_hist(P_trascLen_AB[i], outfile=outFile, label_x= 'Probability', label_y='Frequency', Title=wTitle('Probability distribution','transLen_versus_P_AB_hist_subplot_50'+ext), histtype='stepfilled', bins=20, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()


		####################################
		# Histograms 10 bins
		outFile = path+'transLen_versus_P_BB_hist_subplot_10'+ext
		for i in range(len(P_trascLen_BB)):
			subplot_hist(P_trascLen_BB[i], outfile=outFile, label_x= 'Probability', label_y='Frequency', Title= wTitle('Probability distribution', 'transLen_versus_P_BB_hist_subplot_50'+ext), histtype='stepfilled', bins=10, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_P_AB_hist_subplot_10'+ext
		for i in range(len(P_trascLen_AB)):
			subplot_hist(P_trascLen_AB[i], outfile=outFile, label_x= 'Probability', label_y='Frequency', Title=wTitle('Probability distribution','transLen_versus_P_AB_hist_subplot_50'+ext), histtype='stepfilled', bins=10, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()


		# % change per length
		# % change over various length per substructure



