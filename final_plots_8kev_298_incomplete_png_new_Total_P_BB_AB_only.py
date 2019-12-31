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
ext='_new.pdf'

# 26,24 or 7.5,7 or for NAR: 9,7 (full page) => size_a,size_b
# 19,12 or 7,4.4 or for NAR: 5.2,3.3 (vertical) or 3.3,2.08 (horizontal) => size_c,size_d

size_a,size_b= 7.5,7#26,24
size_c,size_d= 3.3,2.08#19,12
font_size_ratio=6/7.0 # font/A4 page width
subplotHPad=0.5
subplotVPad=0.5
figPad=0.5

font_size= font_size_ratio * size_b
markeredgewidth_ratio=1.0/12.0
marker_size= font_size

mpl.rcParams['font.weight']='bold'
mpl.rcParams['axes.linewidth']='1.5'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['grid.linestyle']= '-'
mpl.rcParams['font.size']=font_size
mpl.rcParams['savefig.bbox']='tight'
mpl.rcParams['savefig.pad_inches']=figPad


mpl.rcParams['xtick.labelsize']='small'      #: medium # fontsize of the tick labels
mpl.rcParams['ytick.labelsize']='small' 

def set_globVar():
	global globVar 
	globVar= 0
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
		mplt.plot(x,y, col, linewidth=6, markersize=marker_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
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


def mplot_twinx(x, y1, label_x= 'x', label_y1='y1', outFile='FreeEnergy_versus_P'+ext, Title='', legends=None, col1='b.', N=111, y2=None, col2='r.',label_y2='y2'):
		'''A plotting fuction for using matPlotLib'''
		print 'Creating plots of free energies and probabilities ..'
		if type(x[0])==str:
			#mplt.xticks(x)
			#fig 	= mplt.figure()
			fig, ax	= mplt.subplots()
			#ax1
			#ax 		= fig.add_subplot(N)
			ax.set_xticklabels(x, rotation=45)
			ax.set_ylabel(label_y1, color=col1[0])
			#ax.tick_params('y', color=col1[0])
			#ax2
			ax_twinx = ax.twinx()
			ax_twinx.set_xticklabels(x, rotation=45)
			ax_twinx.set_ylabel(label_y2, color=col2[0])
			#ax_twinx.tick_params('y', color=col2[0])
			x= np.array(range(len(x)))

		ax.set_xlabel(label_x)
		ax.plot(x,y1, col1, linewidth=1.5, markersize=marker_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
		#mplt.ylabel(label_y1)

		ax_twinx.plot(x,y2, col2, linewidth=1.5, markersize=marker_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
		#ax_twinx.ylabel(label_y2)

		#mplt.tight_layout(pad=subplotPad)
		if legends != None:
			#mplt.legend(legends)
			mplt.legend(legends, bbox_to_anchor=(1.02, 1), shadow=True, fancybox=True, loc=legend_loc, borderaxespad=3)
		ax.grid(True)
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)
		mplt.savefig(outFile, dpi=1200)
		mplt.cla()
		mplt.clf()	
		return

def mplot_err(x, y, label_x= 'x', label_y='y', outFile='FreeEnergy_versus_P'+ext, Title='', legends=None, col1='r.', N=111, yerr=None):
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
		#mplt.plot(x,y, col, linewidth=6, markersize=marker_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
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
		mplt.plot(x,y, col1, linewidth=2, markersize=font_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
		#mplt.plot(x,y, col2)		
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.title(Title)
		mplt.grid(True)
		#mplt.tight_layout(pad=subplotPad)
		if legends != None:
			#mplt.legend(legends)
			mplt.legend(legends, bbox_to_anchor=(1.01, 1), shadow=True, fancybox=True, loc=legend_loc, borderaxespad=0)
		mplt.savefig(outFile, dpi=1200)
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)
		#mplt.cla()
		#mplt.clf()
		return


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



def visTable(path, cwd):
	''' Prints and plots final results tables for a single length'''
	# get data for a single run/file i.e. same length, energy and temperature
	substrucs_data =	read_tables(path)

	substrucs_names		= substrucs_data['energies_header_list'][:-1]
	for i in range(len(substrucs_names)):
		substrucs_names[i] = substrucs_names[i].upper()
	substrucs_names.append('Aptamer_TS')
	return substrucs_data, substrucs_names





if __name__=='__main__':
	set_globVar()
	R	= 0.0019872041
	#T	= 298.15

	# for all main folders
	mainDirs=['e_range_8_298']

	for mainDir in mainDirs:
		T	= float(mainDir[-3:]+'.15')
		cwd	=os.getcwd()
		dirs_incomplete = [ '145_incomplete_AT',
							'146_incomplete_AT',
							'147_incomplete_AT',
							'148_incomplete_AT',
							'149_incomplete_AT',
							'150_incomplete_AT',
							'151_incomplete_AT',
							'152_incomplete_AT']

		trascLen					= []
		trascLen_name				= []
		P_substruct_sum_trascLen_BB	= []
		P_substruct_sum_trascLen_AB	= []
		color						= ('ro-','bo-','ko-','mo-','go-','co-','yo-', 'kx-')
		dirs	= dirs_incomplete

		for dir_name in dirs:
			
			path	= cwd+'/'+mainDir+'/'+dir_name+'/'
			trascLen.append(int(dir_name[:3]))
			trascLen_name.append(dir_name[:3])

			substrucs_data, substrucs_names  = visTable(path, cwd)

			P_substruct_sum_trascLen_BB.append(np.sum(substrucs_data['P_with_Substructs_before_list']))
			P_substruct_sum_trascLen_AB.append(np.sum(substrucs_data['P_with_Substructs_after_list']))

		trascLen= np.array(trascLen)

		P_substruct_sum_trascLen_BB		= np.array(P_substruct_sum_trascLen_BB)
		P_substruct_sum_trascLen_AB		= np.array(P_substruct_sum_trascLen_AB)


		multi_mplot_str(trascLen, P_substruct_sum_trascLen_BB, label_x= 'Transcription length', label_y='Total Probability', 
		outFile='transLen_versus_Total_P_BB_AB_no_inset'+ext, Title='', col1='bo-', counter=0)

		multi_mplot_str(trascLen, P_substruct_sum_trascLen_AB, label_x= 'Transcription length', label_y='Total Probability', 
		outFile='transLen_versus_Total_P_BB_AB_no_inset'+ext, Title='', legends=['Before binding', 'After binding'], col1='ro-', counter=1)
		mplt.cla()
		mplt.clf()
















