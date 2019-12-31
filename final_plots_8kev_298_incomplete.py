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
ext='.pdf'

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

mpl.rcParams['font.weight']='bold'
mpl.rcParams['axes.linewidth']='1.5'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['grid.linestyle']= '-'
mpl.rcParams['font.size']=font_size
mpl.rcParams['savefig.bbox']='tight'
mpl.rcParams['savefig.pad_inches']=figPad

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
		mplt.plot(x, y, col1, linewidth=2, markersize=font_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
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
		mplt.plot(x,y, col, linewidth=2, markersize=font_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
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
		#mplt.plot(x,y, col, linewidth=2, markersize=font_size, markeredgewidth=(markeredgewidth_ratio)*font_size)
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

def plot_surf(x,y,z, xlabel='', ylabel='', zlabel='', outFile='surf'+ext, N=111):
	fig 	= mplt.figure()
	ax 		= fig.add_subplot(N, projection='3d')

	if type(x[0])==str:
		x= np.array(range(len(x)))
		mplt.xticks(x)

	if type(y[0])==str:
		x= np.array(range(len(y)))
		mplt.yticks(y)

	ax.plot_surface(x, y, z, cmap=mplt.cm.hot, rstride=1, cstride=1, linewidth=2, antialiased=False)#flag)#coolwarm)
	#ax.set_zlim(0, 1)
	ax.set_xlabel(xlabel)#r'$\phi_\mathrm{real}$')
	ax.set_ylabel(ylabel)#r'$\phi_\mathrm{im}$')
	ax.set_zlabel(zlabel)#r'$V(\phi)$')
	mplt.savefig(outFile, dpi=600)
	#mplt.show()
	mplt.cla()
	mplt.clf()
	return


def plot_wire(x,y,z, xlabel='', ylabel='', zlabel='', outFile='surf'+ext, N=111):
	fig 	= mplt.figure()
	ax 		= fig.add_subplot(N, projection='3d')

	if type(x[0])==str:
		x= np.array(range(len(x)))
		mplt.xticks(x)

	if type(y[0])==str:
		x= np.array(range(len(y)))
		mplt.yticks(y)

	ax.plot_wireframe(x, y, z, cmap=mplt.cm.hot, rstride=1, cstride=1)#, linewidth=0, antialiased=False)#flag)#coolwarm)
	#ax.set_zlim(0, 1)
	ax.set_xlabel(xlabel)#r'$\phi_\mathrm{real}$')
	ax.set_ylabel(ylabel)#r'$\phi_\mathrm{im}$')
	ax.set_zlabel(zlabel)#r'$V(\phi)$')
	mplt.savefig(outFile, dpi=600)
	#mplt.show()
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
	# Plot free energies and prbablities
	mplot(FreeEnergies, Probabilities, label_x='Free energy [kcal/mol]', label_y='Probability', outFile=outPath+'FreeEnergy_versus_P_before_binding'+ext)

#	dataArr= np.sort(dataArr, axis=0)
	FreeEnergies, Probabilities = dataArr[0], dataArr[1]

	mplot(Probabilities, FreeEnergies, label_y='Free energy [kcal/mol]', label_x='Probability', outFile=outPath+'P_versus_FreeEnergy_before_binding'+ext)
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
	# Plot free energies and prbablities
	mplot(FreeEnergies, Probabilities, label_x='Free energy [kcal/mol]', label_y='Probability', outFile= outPath+'FreeEnergy_versus_P_after_binding'+ext)

#	dataArr= np.sort(dataArr, axis=0)
	FreeEnergies, Probabilities = dataArr[0], dataArr[1]

	mplot(Probabilities, FreeEnergies, label_y='Free energy [kcal/mol]', label_x='Probability', outFile= outPath+'P_versus_FreeEnergy_after_binding'+ext)
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
	mplot(subDirs, substructFreq, label_x= 'Substructure', label_y='Frequency', outFile=path+'Substruct_versus_Frequency'+ext, Title='', legends=None, col='ro-')
	return substructFreq


def visTable(path, cwd):
	''' Prints and plots final results tables for a single length'''
	# get data for a single run/file i.e. same length, energy and temperature
	substrucs_data =	read_tables(path)

	# print P tables
	f=open(path+'p_tables', 'w')
	print >> f, 'Construct','\t','Q before binding','\t', 'P before binding','\t', 'Q after binding', '\t', 'P after binding', '\t', 'dP'

	for i in range(len(substrucs_data['energies_header_list'])):
		print >> f, substrucs_data['energies_header_list'][i],'\t', substrucs_data['Q_before_list'][i],'\t', substrucs_data['P_with_Substructs_before_list'][i],'\t', substrucs_data['Q_after_list'][i],'\t', substrucs_data['P_with_Substructs_after_list'][i],'\t', (substrucs_data['P_with_Substructs_after_list'][i] - substrucs_data['P_with_Substructs_before_list'][i])
	print >> f, 'Total Probabilty', '\t', '-', '\t', np.sum(substrucs_data['P_with_Substructs_before_list']), '\t','-', np.sum(substrucs_data['P_with_Substructs_after_list']),'\t',(np.sum(substrucs_data['P_with_Substructs_after_list']) - np.sum(substrucs_data['P_with_Substructs_before_list']) )
	f.close()

	# print All P tables
	if globVar==0:
		f=open(cwd+'/'+'All_p_tables', 'w')
	else:
		f=open(cwd+'/'+'All_p_tables', 'a')
	print >> f, path
	print >> f, 'Construct','\t','Q before binding','\t', 'P before binding','\t', 'Q after binding', '\t', 'P after binding', '\t', 'dP'

	for i in range(len(substrucs_data['energies_header_list'])):
		print >> f, substrucs_data['energies_header_list'][i],'\t', substrucs_data['Q_before_list'][i],'\t', substrucs_data['P_with_Substructs_before_list'][i],'\t', substrucs_data['Q_after_list'][i],'\t', substrucs_data['P_with_Substructs_after_list'][i],'\t', (substrucs_data['P_with_Substructs_after_list'][i] - substrucs_data['P_with_Substructs_before_list'][i])
	print >> f, 'Total Probabilty', '\t', '-', '\t', np.sum(substrucs_data['P_with_Substructs_before_list']), '\t','-', np.sum(substrucs_data['P_with_Substructs_after_list']),'\t',(np.sum(substrucs_data['P_with_Substructs_after_list']) - np.sum(substrucs_data['P_with_Substructs_before_list']) )
	print >> f, '\n'	
	f.close()


	substrucs_names		= substrucs_data['energies_header_list'][:-1]
	for i in range(len(substrucs_names)):
		substrucs_names[i] = substrucs_names[i].upper()
	substrucs_names.append('Aptamer_TS')


	# b/a 
	mplot(substrucs_names, np.array(substrucs_data['construct_B_over_A_list']), label_x='Construct', label_y='b/a', outFile=path+'Substruct_versus_b-over-a'+ext, Title='', legends=None, col='ro-')
	#b/a_with_err  np.array(substrucs_data['std_dev_B_over_A_list'])
	mplot_err(substrucs_names, np.array(substrucs_data['construct_B_over_A_list']), label_x='Construct', label_y='b/a', yerr=np.array(substrucs_data['std_dev_B_over_A_list']), outFile=path+'Substruct_versus_b-over-a_with_err'+ext, Title='', legends=None, col='ro-')


#	substrucs_names		= 
	#['2P1_11AT', '3P1_10AT', '4P1_9AT', '5P1_8AT', '6P1_7AT', '7P1_6AT', '8P1_5AT', 'Aptamer_TS']
	# Ka/Kd
	mplot(substrucs_names, np.array(substrucs_data['ka_list'])*1.0E-6, label_x='Construct', label_y='Ka x  1.0E6', outFile=path+'Substruct_versus_ka'+ext, Title='', legends=None, col='ro-')
	mplot(substrucs_names, np.array(substrucs_data['kd_list'])*1.0E6, label_x='Construct', label_y='Kd [micromolars]', outFile=path+'Substruct_versus_kd'+ext, Title='', legends=None, col='bo-')
	# Ka & Kd
	multi_mplot_str(substrucs_names, np.array(substrucs_data['ka_list'])*1.0E-6, label_x= 'Construct', label_y='Ka x  1.0E6 / Kd  [micromolars]', outFile=path+'Substruct_versus_ka_kd'+ext, Title='', col1='bo-', col2='r-', counter=0)
	multi_mplot_str(substrucs_names, np.array(substrucs_data['kd_list'])*1.0E6, label_x= 'Construct', label_y='Ka x  1.0E6 / Kd [micromolars]', outFile=path+'Substruct_versus_ka_kd'+ext, Title='', col1='ro-', col2='r-', legends=['Ka','Kd'], counter=1)
	mplt.cla()
	mplt.clf()

	# Binding energies: seperate plots
	mplot(substrucs_names, substrucs_data['dG0_specific_list'], label_x= 'Construct', label_y='Free energy [kcal/mol]', outFile=path+'Substruct_versus_sp_dG_binding'+ext, Title='', legends=None, col='ro-')
	mplot(substrucs_names, substrucs_data['dG0_non_specific_list'], label_x= 'Construct', label_y='Free energy [kcal/mol]', outFile=path+'Substruct_versus_nonsp_dG_binding'+ext, Title='', legends=None, col='bo-')
	mplot(substrucs_names, substrucs_data['dG0_ka_list'], label_x= 'Construct', label_y='Free energy [kcal/mol]', outFile=path+'Substruct_versus_total_dG'+ext, Title='', legends=None, col='mo-')
	
	# Binding energies: all on one plot
	multi_mplot_str(substrucs_names, substrucs_data['dG0_specific_list'], label_x= 'Construct', label_y='Free energy [kcal/mol]', outFile=path+'Substruct_versus_dG'+ext, Title='', col1='ro-', col2='r-', counter=0)
	multi_mplot_str(substrucs_names, substrucs_data['dG0_non_specific_list'], label_x= 'Construct', label_y='Free energy [kcal/mol]', outFile=path+'Substruct_versus_dG'+ext, Title='', col1='bo-', col2='b-', counter=1)
	multi_mplot_str(substrucs_names, substrucs_data['dG0_ka_list'], label_x= 'Construct', label_y='Free energy [kcal/mol]', outFile=path+'Substruct_versus_dG'+ext, Title='', col1='mo-', col2='m-', counter=2,legends=['dG sp', 'dG non-sp', 'dG Total'])
	mplt.cla()
	mplt.clf()

	# Plot overall P for substructures before and after SAM
	mplot(substrucs_names, substrucs_data['P_with_Substructs_before_list'], label_x= 'Substructure', label_y='Probability', outFile=path+'Substruct_versus_P_BB'+ext, Title='', legends=None, col='ro-')
	mplot(substrucs_names, substrucs_data['P_with_Substructs_after_list'], label_x= 'Substructure', label_y='Probability', outFile=path+'Substruct_versus_P_AB'+ext, Title='', legends=None, col='bo-')

 	# Plot overall P for substructures before and after SAM: all on one plot
	multi_mplot_str(substrucs_names, substrucs_data['P_with_Substructs_before_list'], label_x= 'Substructure', label_y='Probability', outFile=path+'Substruct_versus_P_BB_AB'+ext, Title='', col1='ro-', col2='r-', counter=0)
	multi_mplot_str(substrucs_names, substrucs_data['P_with_Substructs_after_list'], label_x= 'Substructure', label_y='Probability', outFile=path+'Substruct_versus_P_BB_AB'+ext, Title='', col1='bo-', col2='b-', counter=1, legends=['Before binding','After Binding'])
	mplt.cla()
	mplt.clf()

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
		f= open('plot_titles', 'w')
		globVar = 1
	else:
		f= open('plot_titles', 'a')
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

		nConf_trascLen			= []
		nConf_substruct_trascLen= []
		nConf_bound_trascLen	= []

		E_ave_trascLen_BB		= []
		E_ave_trascLen_AB		= []
		E_sum_trascLen_BB		= []
		E_sum_trascLen_AB		= []
		E_min_trascLen_BB		= []
		E_max_trascLen_BB		= []
		E_min_trascLen_AB		= []
		E_max_trascLen_AB		= []

		E_substruct_trascLen_BB	= []
		E_substruct_trascLen_AB	= []

		P_trascLen_BB			= []
		P_trascLen_AB			= []
		E_trascLen_BB			= []
		E_trascLen_AB			= []

		P_substruct_trascLen_BB	= []
		P_substruct_trascLen_AB	= []

		P_substruct_sum_trascLen_BB	= []
		P_substruct_sum_trascLen_AB	= []

		P_min_trascLen_BB		= []
		P_min_trascLen_AB		= []
		P_max_trascLen_BB		= []
		P_max_trascLen_AB		= []
		P_ave_trascLen_BB		= []
		P_ave_trascLen_AB		= []
		P_sum_trascLen_BB		= []
		P_sum_trascLen_AB		= []


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

			# Plot free energies vs prbablities: before & after
			multi_mplot(FreeEnergies_BB, Probabilities_BB, label_x='Free energy [kcal/mol]', label_y='Probability', outFile= path+'FreeEnergy_versus_P_BB_AB1'+ext)
			multi_mplot(FreeEnergies_AB, Probabilities_AB, label_x='Free energy [kcal/mol]', label_y='Probability', outFile= path+'FreeEnergy_versus_P_BB_AB2'+ext)
			mplt.cla()
			mplt.clf()

			multi_mplot(Probabilities_BB, FreeEnergies_BB, label_y='Free energy [kcal/mol]', label_x='Probability', outFile= path+'P_versus_FreeEnergy_BB_AB1'+ext)
			multi_mplot(Probabilities_AB, FreeEnergies_AB, label_y='Free energy [kcal/mol]', label_x='Probability', outFile= path+'P_versus_FreeEnergy_BB_AB2'+ext)
			mplt.cla()
			mplt.clf()

			# Prob. difference: Probabilities_AB - Probabilities_BB  # TEX : r'$\alpha > \beta$'
			diff_P = Probabilities_AB - Probabilities_BB
			mplot(np.array(range(len(FreeEnergies_BB))), diff_P, label_x='Conformer', label_y= r'$dP$', outFile= path+'FreeEnergy_versus_P_diff'+ext)

			# Histograms
			plot_hist(Probabilities_BB, outfile=path+'P_hist_before'+ext, label_x= 'Probability', label_y='Frequency', Title= wTitle('Probability distribution', 'P_hist_before'+ext), legends=None, col='r-', histtype='step')
			plot_hist(Probabilities_AB, outfile=path+'P_hist_after'+ext, label_x= 'Probability', label_y='Frequency', Title=wTitle('Probability distribution','P_hist_after'+ext), legends=None, col='r-', histtype='step')

			plot_hist(FreeEnergies_BB, outfile=path+'E_hist_before'+ext, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','E_hist_before'+ext), legends=None, col='b-', histtype='step')
			plot_hist(FreeEnergies_AB, outfile=path+'E_hist_after'+ext, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','E_hist_after'+ext), legends=None, col='b-', histtype='step')

			# frequencies
			substructFreq = readSubstrucsData(path, substrucs_data['energies_header_list'], inFile='indicies')

			# Create overall plot in the form of subplots


			# stats:
			# overall frequences per lengths
			trascLen.append(int(dir_name[:3]))
			trascLen_name.append(dir_name[:3])

			nConf_trascLen.append(len(FreeEnergies_BB))
			print 'len=',len(FreeEnergies_BB)
			nConf_bound_trascLen.append(np.sum(substructFreq))
			print 'sum=',np.sum(substructFreq)
			nConf_substruct_trascLen.append(substructFreq)
			print 'list', substructFreq

			E_ave_trascLen_BB.append(np.sum(FreeEnergies_BB)/len(FreeEnergies_BB))
			E_ave_trascLen_AB.append(np.sum(FreeEnergies_AB)/len(FreeEnergies_AB))
			E_sum_trascLen_BB.append(np.sum(FreeEnergies_BB))
			E_sum_trascLen_AB.append(np.sum(FreeEnergies_AB))
			E_min_trascLen_BB.append(np.min(FreeEnergies_BB))
			E_max_trascLen_BB.append(np.max(FreeEnergies_BB))
			E_min_trascLen_AB.append(np.min(FreeEnergies_AB))
			E_max_trascLen_AB.append(np.max(FreeEnergies_AB))
	
			E_substruct_trascLen_BB.append(-np.log(np.array(substrucs_data['P_with_Substructs_before_list']) * np.array(substrucs_data['Q_before_list']))*R*T)
			E_substruct_trascLen_AB.append(-np.log (np.array(substrucs_data['P_with_Substructs_after_list']) * np.array(substrucs_data['Q_before_list']))*R*T)
			
			print E_substruct_trascLen_BB
			print E_substruct_trascLen_AB

			P_trascLen_BB.append(Probabilities_BB)
			P_trascLen_AB.append(Probabilities_AB)
			E_trascLen_BB.append(FreeEnergies_BB)
			E_trascLen_AB.append(FreeEnergies_AB)

			P_substruct_trascLen_BB.append(substrucs_data['P_with_Substructs_before_list'])
			P_substruct_trascLen_AB.append(substrucs_data['P_with_Substructs_after_list'])

			# 	Sum of probabilities of bound (substructures)
			P_substruct_sum_trascLen_BB.append(np.sum(substrucs_data['P_with_Substructs_before_list']))
			P_substruct_sum_trascLen_AB.append(np.sum(substrucs_data['P_with_Substructs_after_list']))

			# Sum of bound and unbound ==1
			P_sum_trascLen_BB.append(np.sum(Probabilities_BB))
			P_sum_trascLen_AB.append(np.sum(Probabilities_AB))
			P_ave_trascLen_BB.append(np.sum(Probabilities_BB)/len(Probabilities_BB))
			P_ave_trascLen_AB.append(np.sum(Probabilities_AB)/len(Probabilities_AB))
			# min and max subtruct must include in text
			P_min_trascLen_BB.append(np.min(Probabilities_BB))
			P_min_trascLen_AB.append(np.min(Probabilities_AB))
			P_max_trascLen_BB.append(np.max(Probabilities_BB))
			P_max_trascLen_AB.append(np.max(Probabilities_AB))

			# nfold change per substruct
			# Hist TYPE ... done
			# Frequencies  ... done
			# Subplots
			# Binding data ... done
			# Overall substruct P ... done
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
		nConf_trascLen			= np.array(nConf_trascLen)

		#################################################################################################### Conformers
		# Total No of conformers
		print 'nConf_trascLen=', nConf_trascLen
		mplot(trascLen, nConf_trascLen*1.0E-6, label_x= 'Transcription length', label_y='No. of Conformers (in millions)', outFile=path+'transLen_versus_nConf'+ext, Title=wTitle('No. of conformers at increasing transcription lengths', 'transLen_versus_nConf'+ext), legends=None, col='ro-')


		########################## series plot # No. of conformers per intermediate
		nConf_substruct_trascLen= np.array(nConf_substruct_trascLen)
		# regular
		for i in range(len(nConf_substruct_trascLen)):
			multi_mplot_str(substrucs_names, nConf_substruct_trascLen[i]*1.0E-3, label_x= 'Substructure', label_y='No. of Conformers (in thousands)', outFile=path+'Substructure_transLen_versus_nConf_inter'+ext, Title=wTitle('No. of conformers per intermediate at increasing transcription lengths', 'Substructure_transLen_versus_nConf_inter'+ext), col1=color[i], legends=trascLen_name, counter=i)
		mplt.cla()
		mplt.clf()

		# subplot
		outFile = path+'transLen_versus_nConf_inter_subplot'+ext
		for i in range(len(nConf_substruct_trascLen)):
			multi_subplot(substrucs_names, nConf_substruct_trascLen[i]*1.0E-3, label_x='', label_y='No. of Conformers (in thousands)', outFile=outFile, Title=wTitle('No. of conformers per intermediate at increasing transcription lengths', 'transLen_versus_nConf_inter_subplot'+ext), col1=color[i], legends=trascLen_name,  rows=len(nConf_substruct_trascLen), col=1, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		##########################

		# No of bound conformers
		nConf_bound_trascLen	= np.array(nConf_bound_trascLen)
		mplot(trascLen, nConf_bound_trascLen*1.0E-6, label_x= 'Transcription length', label_y='No. of bound Conformers (in millions)', outFile=path+'transLen_versus_nConf_bound'+ext, Title=wTitle('No. of conformers bound at increasing transcription lengths','transLen_versus_nConf_bound'+ext), legends=None, col='ro-')

		# Fraction of bound conformers
		nConf_trascLen			= np.array(nConf_trascLen, dtype='float64')
		nConf_bound_trascLen	= np.array(nConf_bound_trascLen, dtype='float64')
		fraction_bound			= nConf_bound_trascLen/nConf_trascLen
		mplot(trascLen, fraction_bound, label_x= 'Transcription length', label_y='Fraction of bound Conformers', outFile=path+'transLen_versus_fraction_bound'+ext, Title=wTitle('Fraction of conformers bound at increasing transcription lengths','transLen_versus_fraction_bound'+ext), legends=None, col='ro-')

		print '######################################'

		## printing tables

		# overall stats
		f_stat1=open('stats1', 'w')

		print >> f_stat1, 'Transcription length','\t', 'Total no. of Conformers','\t', 'No. of bound Conformers','\t', 'Fraction of bound', '\t', 'Unbound fraction'
		for i in range(len(trascLen)):
			print >> f_stat1, trascLen[i], '\t', nConf_trascLen[i], '\t', nConf_bound_trascLen[i],'\t', fraction_bound[i], '\t', 1.0-fraction_bound[i]
		f_stat1.close()

		# substruc stats
		f_stat2=open('stats2', 'w')
		for j in range(len(trascLen)):
			print >> f_stat2, 'Transcription length=', trascLen[j]
			print >> f_stat2, 'Substructure','\t','No. of conformers per substructure','\t', 'Contributing fraction from total','\t', 'Contributing fraction from bound'
			for i in range(len(substrucs_names)):
				print >> f_stat2, substrucs_names[i], '\t', nConf_substruct_trascLen[j][i], '\t',  nConf_substruct_trascLen[j][i]/nConf_trascLen[j], '\t', nConf_substruct_trascLen[j][i]/nConf_bound_trascLen[j]
			print  >> f_stat2, '\n'
		f_stat2.close()

		print '######################################'


		########################## series plot # Fraction (from total) of conformers per intermediate
		nConf_substruct_trascLen= np.array(nConf_substruct_trascLen, dtype='float64')
		# regular
		for i in range(len(nConf_substruct_trascLen)):
			multi_mplot_str(substrucs_names, nConf_substruct_trascLen[i]/nConf_trascLen[i], label_x= 'Substructure', label_y='Fraction of Conformers', outFile=path+'Substructure_transLen_versus_Fraction_inter'+ext, Title=wTitle('Fraction of conformers per intermediate at increasing transcription lengths', 'Substructure_transLen_versus_Fraction_inter'+ext), col1=color[i], legends=trascLen_name, counter=i)
		mplt.cla()
		mplt.clf()

		# subplot
		outFile = path+'transLen_versus_nConf_inter_subplot_fraction'+ext
		for i in range(len(nConf_substruct_trascLen)):
			multi_subplot(substrucs_names, nConf_substruct_trascLen[i]/nConf_trascLen[i], label_x='', label_y='Fraction of Conformers', outFile=outFile, Title=wTitle('Fraction of conformers per intermediate at increasing transcription lengths', 'transLen_versus_Fraction_inter_subplot'+ext), col1=color[i], legends=trascLen_name,  rows=len(nConf_substruct_trascLen), col=1, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		########################## series plot # Fraction (from bound) of conformers per intermediate
		# regular
		for i in range(len(nConf_substruct_trascLen)):
			multi_mplot_str(substrucs_names, nConf_substruct_trascLen[i]/nConf_bound_trascLen[i], label_x= 'Substructure', label_y='Fraction (from bound) of Conformers', outFile=path+'Substructure_transLen_versus_Fraction_from_bound_inter'+ext, Title=wTitle('Fraction (from bound) of conformers per intermediate at increasing transcription lengths', 'Substructure_transLen_versus_Fraction_from_bound_inter'+ext), col1=color[i], legends=trascLen_name, counter=i)
		mplt.cla()
		mplt.clf()

		# subplot
		outFile = path+'transLen_versus_nConf_inter_subplot_fraction_from_bound'+ext
		for i in range(len(nConf_substruct_trascLen)):
			multi_subplot(substrucs_names, nConf_substruct_trascLen[i]/nConf_bound_trascLen[i], label_x='', label_y='Fraction (from bound) of Conformers', outFile=outFile, Title=wTitle('Fraction (from bound) of conformers per intermediate at increasing transcription lengths', 'transLen_versus_Fraction_from_bound_inter_subplot'+ext), col1=color[i], legends=trascLen_name,  rows=len(nConf_substruct_trascLen), col=1, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		#################################################################################################### Energies

		# Energy stats
		E_ave_trascLen_BB		= np.array(E_ave_trascLen_BB)
		mplot(trascLen, E_ave_trascLen_BB, label_x= 'Transcription length', label_y='Average free energy [kcal/mol] (before binding)', outFile=path+'transLen_versus_ave_E_BB'+ext, Title=wTitle('Average folding free energy before binding at increasing transcription lengths', 'transLen_versus_ave_E_BB'+ext), legends=None, col='go-')

		E_ave_trascLen_AB		= np.array(E_ave_trascLen_AB)
		mplot(trascLen, E_ave_trascLen_AB, label_x= 'Transcription length', label_y='Average free energy [kcal/mol] (after binding)', outFile=path+'transLen_versus_ave_E_AB'+ext, Title=wTitle('Average folding free energy after binding at increasing transcription lengths','transLen_versus_ave_E_AB'+ext), legends=None, col='ro-')

		E_sum_trascLen_BB		= np.array(E_sum_trascLen_BB)
		mplot(trascLen, E_sum_trascLen_BB, label_x= 'Transcription length', label_y='Total free energies [kcal/mol] (before binding)', outFile=path+'transLen_versus_sum_E_BB'+ext, Title=wTitle('Total folding free energy before binding at increasing transcription lengths', 'transLen_versus_sum_E_BB'+ext), legends=None, col='go-')

		E_sum_trascLen_AB		= np.array(E_sum_trascLen_AB)
		mplot(trascLen, E_sum_trascLen_AB, label_x= 'Transcription length', label_y='Total free energy [kcal/mol] (after binding)', outFile=path+'transLen_versus_sum_E_AB'+ext, Title=wTitle('Total folding free energy after binding at increasing transcription lengths', 'transLen_versus_sum_E_AB'+ext), legends=None, col='ro-')

		E_min_trascLen_BB		= np.array(E_min_trascLen_BB)
		mplot(trascLen, E_min_trascLen_BB, label_x= 'Transcription length', label_y='Minimum free energy [kcal/mol] (before binding)', outFile=path+'transLen_versus_min_E_AB'+ext, Title=wTitle('Minimum folding free energy before binding at increasing transcription lengths', 'transLen_versus_min_E_AB'+ext), legends=None, col='go-')

		E_max_trascLen_BB		= np.array(E_max_trascLen_BB)
		mplot(trascLen, E_max_trascLen_BB, label_x= 'Transcription length', label_y='Maximum free energy [kcal/mol] (before binding)', outFile=path+'transLen_versus_max_E_AB'+ext, Title=wTitle('Maximum folding free energy before binding at increasing transcription lengths', 'transLen_versus_max_E_AB'+ext), legends=None, col='go-')

		E_min_trascLen_AB		= np.array(E_min_trascLen_AB)
		mplot(trascLen, E_min_trascLen_AB, label_x= 'Transcription length', label_y='Minimum free energy [kcal/mol] (after binding)', outFile=path+'transLen_versus_min_E_AB'+ext, Title=wTitle('Minimum folding free energy after binding at increasing transcription lengths', 'transLen_versus_min_E_AB'+ext), legends=None, col='ro-')

		E_max_trascLen_AB		= np.array(E_max_trascLen_AB)
		mplot(trascLen, E_max_trascLen_AB, label_x= 'Transcription length', label_y='Maximum free energy [kcal/mol] (after binding)', outFile=path+'transLen_versus_max_E_AB'+ext, Title=wTitle('Maximum folding free energy after binding at increasing transcription lengths', 'transLen_versus_max_E_AB'+ext), legends=None, col='ro-')

		########################## series plot
		E_substruct_trascLen_BB	= np.array(E_substruct_trascLen_BB)
		E_substruct_trascLen_BB	= np.where(E_substruct_trascLen_BB==np.inf, np.zeros_like(E_substruct_trascLen_BB), E_substruct_trascLen_BB)

		# regular
		for i in range(len(E_substruct_trascLen_BB)):
			multi_mplot_str(substrucs_names, E_substruct_trascLen_BB[i], label_x= 'Substructure', label_y='Free energy [kcal/mol] (before binding)', outFile=path+'transLen_versus_E_sub_BB'+ext, Title=wTitle('E of intermediates before binding at increasing transcription lengths', 'transLen_versus_E_sub_BB'+ext), col1=color[i], legends=trascLen_name, counter=i)
		mplt.cla()
		mplt.clf()

		# subplot
		outFile = path+'transLen_versus_E_sub_BB_subplot'+ext
		for i in range(len(E_substruct_trascLen_BB)):
			multi_subplot(substrucs_names, E_substruct_trascLen_BB[i], label_x= '', label_y='Free energy [kcal/mol] (before binding)', outFile=outFile, Title=wTitle('E of intermediates before binding at increasing transcription lengths', 'transLen_versus_E_sub_BB_subplot'+ext), col1=color[i], legends=trascLen_name, rows=len(E_substruct_trascLen_BB), col=1, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		########################## series plot
		E_substruct_trascLen_AB	= np.array(E_substruct_trascLen_AB)
		E_substruct_trascLen_AB	= np.where(E_substruct_trascLen_AB==np.inf, np.zeros_like(E_substruct_trascLen_AB), E_substruct_trascLen_AB)

		# regular
		for i in range(len(E_substruct_trascLen_AB)):
			multi_mplot_str(substrucs_names, E_substruct_trascLen_AB[i], label_x= 'Substructure', label_y='Free energy [kcal/mol] (after binding)', outFile=path+'transLen_versus_E_sub_AB'+ext, Title=wTitle('E of intermediates after binding at increasing transcription lengths', 'transLen_versus_E_sub_AB'+ext), col1=color[i], legends=trascLen_name, counter=i)
		mplt.cla()
		mplt.clf()

		# subplot
		outFile = path+'transLen_versus_E_sub_AB_subplot'+ext
		for i in range(len(E_substruct_trascLen_AB)):
			multi_subplot(substrucs_names, E_substruct_trascLen_AB[i], label_x= '', label_y='Free energy [kcal/mol] (after binding)', outFile=outFile, Title=wTitle('E of intermediates after binding at increasing transcription lengths', 'transLen_versus_E_sub_AB_subplot'+ext), col1=color[i], legends=trascLen_name, rows=len(E_substruct_trascLen_AB), col=1, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		#################################################################################################### Probabilities

		########################## series plot
		P_trascLen_BB			= np.array(P_trascLen_BB)
		# regular
		for i in range(len(P_trascLen_BB)):
			multi_mplot(E_trascLen_BB[i], P_trascLen_BB[i], label_x='Free energy [kcal/mol]', label_y='P (before binding)', outFile=path+'transLen_versus_P_BB'+ext, Title=wTitle('P before binding at increasing transcription lengths', 'transLen_versus_P_BB'+ext), col1=color[i], legends=trascLen_name)
		mplt.cla()
		mplt.clf()

		# subplot
		outFile = path+'transLen_versus_P_BB_subplot'+'.pdf'#ext
		for i in range(len(P_trascLen_BB)):
			multi_subplot(E_trascLen_BB[i], P_trascLen_BB[i], label_x='Free energy [kcal/mol]', label_y='P (before binding)', outFile=outFile, Title=wTitle('P before binding at increasing transcription lengths', 'transLen_versus_P_BB_subplot'+ext), col1=color_dot[i], legends=trascLen_name, rows=len(P_trascLen_BB), col=1, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		########################## series plot
		# subplot both BB and AB
		outFile = path+'transLen_versus_P_BB_AB_subplot'+'.pdf'#ext
		for i in range(len(P_trascLen_BB)):
			multi_subplot(E_trascLen_BB[i], P_trascLen_BB[i], label_x='Free energy [kcal/mol]', label_y='P', outFile=outFile, Title=wTitle('P binding at increasing transcription lengths', 'transLen_versus_P_BB__AB_subplot'+ext), col1=color_dot[1], legends=trascLen_name, rows=len(P_trascLen_BB), col=1, n=i, ni=0)
			multi_subplot(E_trascLen_AB[i], P_trascLen_AB[i], label_x='Free energy [kcal/mol]', label_y='P', outFile=outFile, Title=wTitle('P binding at increasing transcription lengths', 'transLen_versus_P_BB__AB_subplot'+ext), col1=color_dot[0], legends=trascLen_name, rows=len(P_trascLen_AB), col=1, n=i, ni=1)
			mplt.legend(['Before binding','After binding'])
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()


		'''
		########################## series plot 
		# subplot both BB and AB
?		outFile = path+'transLen_versus_P_BB_AB_subplot'+ext
?		conformers= np.(range(E_trascLen_BB[i]))+1
?
?		for i in range(len(E_trascLen_BB)):
?			multi_subplot(conformers, E_trascLen_BB[i], label_y='Free energy [kcal/mol]', label_x='Conformer', outFile=outFile, Title=wTitle('Conformer Vs E at increasing transcription lengths', 'Conf_versus_E_BB_subplot'+ext), col1=color_dot[1], legends=trascLen_name, rows=len(P_trascLen_BB), col=1, n=i, ni=0)
?			mplt.legend(substrucs_names)
?		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()


?		for i in range(len(E_trascLen_AB)):
?			multi_subplot(conformers, E_trascLen_AB[i], label_y='Free energy [kcal/mol]', label_x='Conformer', outFile=outFile, Title=wTitle('Conformer Vs E increasing transcription lengths', 'Conf_versus_E_AB_subplot'+ext), col1=color_dot[0], legends=trascLen_name, rows=len(P_trascLen_AB), col=1, n=i, ni=1)
?			mplt.legend(substrucs_names)
?		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()
		#####################################
		'''



		########################## series plot
		P_trascLen_AB			= np.array(P_trascLen_AB)
		# regular
		for i in range(len(P_trascLen_AB)):
			multi_mplot(E_trascLen_AB[i], P_trascLen_AB[i], label_x='Free energy [kcal/mol]', label_y='P (after binding)', outFile=path+'transLen_versus_P_AB'+ext, Title=wTitle('P after binding at increasing transcription lengths','transLen_versus_P_AB'+ext), col1=color_dot[i], legends=trascLen_name)
		mplt.cla()
		mplt.clf()

		# subplot
		outFile = path+'transLen_versus_P_AB_subplot'+'.pdf'#ext
		for i in range(len(P_trascLen_AB)):
			multi_subplot(E_trascLen_AB[i], P_trascLen_AB[i], label_x='Free energy [kcal/mol]', label_y='P (after binding)', outFile=outFile, Title=wTitle('P after binding at increasing transcription lengths','transLen_versus_P_AB_subplot'+ext), col1=color_dot[i], legends=trascLen_name, rows=len(P_trascLen_AB), col=1, n=i)
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

		outFile = path+'transLen_versus_E_BB_hist_subplot_500'+ext
		for i in range(len(E_trascLen_BB)):
#			subplot_hist(E_trascLen_BB[i], outfile=path+'E_hist_before'+ext, label_x= r'$dG [kcal/mol]$', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_BB_hist_subplot'+ext), histtype='stepfilled', bins=500, facecolor=colors[i], legends=trascLen_name, n=i)
			subplot_hist(E_trascLen_BB[i], outfile=path+'E_hist_before'+ext, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_BB_hist_subplot'+ext), histtype='stepfilled', bins=500, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_E_AB_hist_subplot_500'+ext
		for i in range(len(E_trascLen_AB)):
			subplot_hist(E_trascLen_AB[i], outfile=path+'E_hist_after'+ext, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_AB_hist_subplot'+ext), histtype='stepfilled', bins=500, facecolor=colors[i], legends=trascLen_name, n=i)
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

		outFile = path+'transLen_versus_E_BB_hist_subplot_150'+ext
		for i in range(len(E_trascLen_BB)):
			subplot_hist(E_trascLen_BB[i], outfile=outFile, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_BB_hist_subplot_150'+ext), histtype='stepfilled', bins=150, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_E_AB_hist_subplot_150'+ext
		for i in range(len(E_trascLen_AB)):
			subplot_hist(E_trascLen_AB[i], outfile=outFile, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_AB_hist_subplot_150'+ext), histtype='stepfilled', bins=150, facecolor=colors[i], legends=trascLen_name, n=i)
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

		outFile = path+'transLen_versus_E_BB_hist_subplot_50'+ext
		for i in range(len(E_trascLen_BB)):
			subplot_hist(E_trascLen_BB[i], outfile=outFile, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_BB_hist_subplot_50'+ext), histtype='stepfilled', bins=50, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_E_AB_hist_subplot_50'+ext
		for i in range(len(E_trascLen_AB)):
			subplot_hist(E_trascLen_AB[i], outfile=outFile, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_AB_hist_subplot_50'+ext), histtype='stepfilled', bins=50, facecolor=colors[i], legends=trascLen_name, n=i)
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

		outFile = path+'transLen_versus_E_BB_hist_subplot_20'+ext
		for i in range(len(E_trascLen_BB)):
			subplot_hist(E_trascLen_BB[i], outfile=outFile, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_BB_hist_subplot_50'+ext), histtype='stepfilled', bins=20, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_E_AB_hist_subplot_20'+ext
		for i in range(len(E_trascLen_AB)):
			subplot_hist(E_trascLen_AB[i], outfile=outFile, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_AB_hist_subplot_50'+ext), histtype='stepfilled', bins=20, facecolor=colors[i], legends=trascLen_name, n=i)
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

		outFile = path+'transLen_versus_E_BB_hist_subplot_10'+ext
		for i in range(len(E_trascLen_BB)):
			subplot_hist(E_trascLen_BB[i], outfile=outFile, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_BB_hist_subplot_50'+ext), histtype='stepfilled', bins=10, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		outFile = path+'transLen_versus_E_AB_hist_subplot_10'+ext
		for i in range(len(E_trascLen_AB)):
			subplot_hist(E_trascLen_AB[i], outfile=outFile, label_x= 'Free energy [kcal/mol]', label_y='Frequency', Title=wTitle('Free energy distribution','transLen_versus_E_AB_hist_subplot_50'+ext), histtype='stepfilled', bins=10, facecolor=colors[i], legends=trascLen_name, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()


		########################## series plot
		P_substruct_trascLen_BB	= np.array(P_substruct_trascLen_BB)
		# regular
		for i in range(len(P_substruct_trascLen_BB)):
			multi_mplot_str(substrucs_names, P_substruct_trascLen_BB[i], label_x= 'Substructure', label_y='P (before binding)', outFile=path+'transLen_versus_P_sub_BB'+ext, Title=wTitle('P of intermediates before binding at increasing transcription lengths', 'transLen_versus_P_sub_BB'+ext), col1=color[i], legends=trascLen_name, counter=i)
		mplt.cla()
		mplt.clf()

		# subplot
		outFile = path+'transLen_versus_P_sub_BB_subplot'+ext
		for i in range(len(P_substruct_trascLen_BB)):
			multi_subplot(substrucs_names, P_substruct_trascLen_BB[i], label_x= '', label_y='P (before binding)', outFile=outFile, Title=wTitle('P of intermediates before binding at increasing transcription lengths', 'transLen_versus_P_sub_BB_subplot'+ext), col1=color[i], legends=trascLen_name, rows=len(P_substruct_trascLen_BB), col=1, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		########################## series plot
		P_substruct_trascLen_AB	= np.array(P_substruct_trascLen_AB)
		#regular
		for i in range(len(P_substruct_trascLen_AB)):
			multi_mplot_str(substrucs_names, P_substruct_trascLen_AB[i], label_x= 'Substructure', label_y='P (after binding)', outFile=path+'transLen_versus_P_sub_AB'+ext, Title=wTitle('P of intermediates after binding at increasing transcription lengths','transLen_versus_P_sub_AB'+ext), col1=color[i], legends=trascLen_name, counter=i)
		mplt.cla()
		mplt.clf()

		#subplot
		outFile = path+'transLen_versus_P_sub_AB_subplot'+ext
		for i in range(len(P_substruct_trascLen_AB)):
			multi_subplot(substrucs_names, P_substruct_trascLen_AB[i], label_x= '', label_y='P (after binding)', outFile=outFile, Title=wTitle('P of intermediates after binding at increasing transcription lengths','transLen_versus_P_sub_AB_subplot'+ext), col1=color[i], legends=trascLen_name,  rows=len(P_substruct_trascLen_AB), col=1, n=i)
		f=mplt.gcf()
		f.set_size_inches(size_a,size_b)	
		mplt.savefig(outFile, dpi=600)
		mplt.cla()
		mplt.clf()

		####################################################



		P_ave_trascLen_BB		= np.array(P_ave_trascLen_BB)
		mplot(trascLen, P_ave_trascLen_BB, label_x= 'Transcription length', label_y='Average Probability (before binding)', outFile=path+'transLen_versus_ave_P_BB'+ext, Title=wTitle('Average P before binding at increasing transcription lengths','transLen_versus_ave_P_BB'+ext), legends=None, col='bo-')

		P_min_trascLen_AB		= np.array(P_min_trascLen_AB)
		mplot(trascLen, P_min_trascLen_AB, label_x= 'Transcription length', label_y='Average Probability (after binding)', outFile=path+'transLen_versus_ave_P_AB'+ext, Title=wTitle('Average P after binding at increasing transcription lengths', 'transLen_versus_ave_P_AB'+ext), legends=None, col='ro-')



		P_sum_trascLen_BB		= np.array(P_sum_trascLen_BB)
		mplot(trascLen, P_sum_trascLen_BB, label_x= 'Transcription length', label_y='Sum of Probability (before binding)', outFile=path+'transLen_versus_sum_P_BB'+ext, Title=wTitle('Sum of P before binding at increasing transcription lengths','transLen_versus_sum_P_BB'+ext), legends=None, col='bo-')

		P_sum_trascLen_AB		= np.array(P_sum_trascLen_AB)
		mplot(trascLen, P_sum_trascLen_AB, label_x= 'Transcription length', label_y='Sum of Probability (after binding)', outFile=path+'transLen_versus_sum_P_AB'+ext, Title=wTitle('Sum of P after binding at increasing transcription lengths', 'transLen_versus_sum_P_AB'+ext), legends=None, col='ro-')




		P_min_trascLen_BB		= np.array(P_min_trascLen_BB)
		mplot(trascLen, P_min_trascLen_BB, label_x= 'Transcription length', label_y='Minimum Probability (before binding)', outFile=path+'transLen_versus_min_P_BB'+ext, Title=wTitle('Minimum P before binding at increasing transcription lengths','transLen_versus_min_P_BB'+ext), legends=None, col='bo-')

		P_min_trascLen_AB		= np.array(P_min_trascLen_AB)
		mplot(trascLen, P_min_trascLen_AB, label_x= 'Transcription length', label_y='Minimum Probability (after binding)', outFile=path+'transLen_versus_min_P_AB'+ext, Title=wTitle('Minimum P after binding at increasing transcription lengths', 'transLen_versus_min_P_AB'+ext), legends=None, col='ro-')

		P_max_trascLen_BB		= np.array(P_max_trascLen_BB)
		mplot(trascLen, P_max_trascLen_BB, label_x= 'Transcription length', label_y='Maximum Probability (before binding)', outFile=path+'transLen_versus_max_P_BB'+ext, Title=wTitle('Maximum P before binding at increasing transcription lengths', 'transLen_versus_max_P_BB'+ext), legends=None, col='bo-')

		P_max_trascLen_AB		= np.array(P_max_trascLen_AB)
		mplot(trascLen, P_max_trascLen_AB, label_x= 'Transcription length', label_y='Maximum Probability (after binding)', outFile=path+'transLen_versus_max_P_AB'+ext, Title=wTitle('Maximum P after binding at increasing transcription lengths','transLen_versus_max_P_AB'+ext), legends=None, col='ro-')


		## needs number formatting/approx, should be all ones
		P_substruct_sum_trascLen_BB		= np.array(P_substruct_sum_trascLen_BB)
		mplot(trascLen, P_substruct_sum_trascLen_BB, label_x= 'Transcription length', label_y='Total Probability (before binding)', outFile=path+'transLen_versus_Total_P_BB'+ext, Title=wTitle('Total P of conformers containing substructs before binding at increasing transcription lengths','transLen_versus_Total_P_BB'+ext), legends=None, col='bo-')

		P_substruct_sum_trascLen_AB		= np.array(P_substruct_sum_trascLen_AB)
		mplot(trascLen, P_substruct_sum_trascLen_AB, label_x= 'Transcription length', label_y='Total Probability (after binding)', outFile=path+'transLen_versus_Total_P_AB'+ext, Title=wTitle('Total P of conformers containing substructs after binding at increasing transcription lengths','transLen_versus_Total_P_AB'+ext), legends=None, col='ro-')


		# % change per length
		# Plot overall P for substructures at various lengths .. done
		# % change over various length per substructure
		# Average free energy and minimal free energy at various lengths  .. done

		x, y = np.meshgrid(np.array(range(len(substrucs_names))), trascLen)

		print x
		print y

		for i in range(len(x)):
			print x[i], y[i], P_substruct_trascLen_AB[i]

		# surface plots
		plot_surf(x,y,P_substruct_trascLen_BB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='P_BB_surf'+ext, N=111)
		plot_surf(x,y,P_substruct_trascLen_AB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='P_AB_surf'+ext, N=111)

		plot_surf(x,y,np.log10(E_substruct_trascLen_BB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='logE_BB_surf'+ext, N=111)
		plot_surf(x,y,np.log10(E_substruct_trascLen_AB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='logE_AB_surf'+ext, N=111)
	
		plot_surf(x,y,E_substruct_trascLen_BB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='E_BB_surf'+ext, N=111)
		plot_surf(x,y,E_substruct_trascLen_AB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='E_AB_surf'+ext, N=111)

		plot_surf(x,y,np.log10(P_substruct_trascLen_BB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='logP_BB_surf'+ext, N=111)
		plot_surf(x,y,np.log10(P_substruct_trascLen_AB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='logP_AB_surf'+ext, N=111)

		# wire plots
		plot_wire(x,y,P_substruct_trascLen_BB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='P_BB_wire'+ext, N=111)
		plot_wire(x,y,P_substruct_trascLen_AB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='P_AB_wire'+ext, N=111)

		plot_wire(x,y,np.log10(E_substruct_trascLen_BB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='logE_BB_wire'+ext, N=111)
		plot_wire(x,y,np.log10(E_substruct_trascLen_AB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='logE_AB_wire'+ext, N=111)
	
		plot_wire(x,y,E_substruct_trascLen_BB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='E_BB_wire'+ext, N=111)
		plot_wire(x,y,E_substruct_trascLen_AB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='E_AB_wire'+ext, N=111)

		plot_wire(x,y,np.log10(P_substruct_trascLen_BB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='logP_BB_wire'+ext, N=111)
		plot_wire(x,y,np.log10(P_substruct_trascLen_AB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='logP_AB_wire'+ext, N=111)

		#####################
		'''
		# Full Data surface plots
		plot_surf(x,y,P_trascLen_BB, xlabel='Conformer', ylabel='Transcription Length', zlabel='Probability', outFile='all_P_BB_surf'+ext, N=111)
		plot_surf(x,y,P_trascLen_AB, xlabel='Conformer', ylabel='Transcription Length', zlabel='Probability', outFile='all_P_AB_surf'+ext, N=111)

		plot_surf(x,y,np.log10(E_trascLen_BB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_logE_BB_surf'+ext, N=111)
		plot_surf(x,y,np.log10(E_trascLen_AB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_logE_AB_surf'+ext, N=111)
	
		plot_surf(x,y,E_trascLen_BB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='E_BB_surf'+ext, N=111)
		plot_surf(x,y,E_trascLen_AB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='E_AB_surf'+ext, N=111)

		plot_surf(x,y,np.log10(P_trascLen_BB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_logP_BB_surf'+ext, N=111)
		plot_surf(x,y,np.log10(P_trascLen_AB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_logP_AB_surf'+ext, N=111)
	
		# Full data wire plots
		plot_wire(x,y,P_trascLen_BB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_P_BB_wire'+ext, N=111)
		plot_wire(x,y,P_trascLen_AB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_P_AB_wire'+ext, N=111)

		plot_wire(x,y,np.log10(E_trascLen_BB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_logE_BB_wire'+ext, N=111)
		plot_wire(x,y,np.log10(E_trascLen_AB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_logE_AB_wire'+ext, N=111)
	
		plot_wire(x,y,E_trascLen_BB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_E_BB_wire'+ext, N=111)
		plot_wire(x,y,E_trascLen_AB, xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_E_AB_wire'+ext, N=111)

		plot_wire(x,y,np.log10(P_trascLen_BB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_logP_BB_wire'+ext, N=111)
		plot_wire(x,y,np.log10(P_trascLen_AB), xlabel='Subtructure', ylabel='Transcription Length', zlabel='Probability', outFile='all_logP_AB_wire'+ext, N=111)

		'''

	
	##### Energy range comparisons

	# comparison of energy ranges
	# % change in P as a function of energy range

	##### Temperature Effects

	# at different temperatures ???



	#mplt.xticks(range(len(substrucs_names)),substrucs_names,rotation=45)
	#mplt.legend(legends, shadow=True, fancybox=True, loc=legend_loc, borderaxespad=0)
