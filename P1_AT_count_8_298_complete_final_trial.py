import matplotlib.pyplot as mplt
import cdecimal
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from RiboSwitch_v15 import Vienna
import sys
import os

# Setting numbers percisions upto 100 decima places
myContext = cdecimal.Context(prec=100, rounding='ROUND_HALF_DOWN')
cdecimal.setcontext(myContext)
D= cdecimal.Decimal
ext='.tiff'
dpi=300

# 26,24 or 7.5,7
# 19,12 or 7,4.4

#mpl.rcParams['font.weight']='bold'
#mpl.rcParams['axes.linewidth']='1.5'
#mpl.rcParams['axes.labelweight'] = 'bold'
#mpl.rcParams['grid.linestyle']= '-'

#mpl.rcParams['xtick.labelsize']='small'
#mpl.rcParams['ytick.labelsize']='small'

#mpl.rcParams['axes.labelsize']='small'
#mpl.rcParams['lines.markeredgewidth']='0'


#mpl.rcParams['axes.titlesize']='medium'


#mpl.rcParams['ztick.labelsize']='small'


size_a,size_b= 7.5,7#26,24
size_c,size_d= 3.3,2.08#19,12
font_size_ratio=6/7.0 # font/A4 page width
subplotHPad=0.5
subplotVPad=0.5
figPad=0.5


font_size= font_size_ratio * size_b
size_a,size_b= 3.3,2.08#19,12

#markeredgewidth_ratio=1.0/12.0

mpl.rcParams['font.weight']='bold'
mpl.rcParams['axes.linewidth']='1.5'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['grid.linestyle']= '-'
mpl.rcParams['font.size']=font_size
mpl.rcParams['savefig.bbox']='tight'
mpl.rcParams['savefig.pad_inches']=figPad
mpl.rcParams['lines.markeredgewidth']=0
mpl.rcParams['legend.borderaxespad']=1.5
mpl.rcParams['mathtext.fontset']='custom'
mpl.rcParams['mathtext.bf']='serif:bold'




def plot_surf(x,y,z, xlabel='', ylabel='', zlabel='', title='', outFile='3Dsurface'+ext, N=1, row=1, col=1, c='r', m='.'):

	if N==0:
		fig 	= mplt.figure()
	else:
		fig 	= mplt.gcf()
	ax 		= fig.add_subplot(row, col, N+1, projection='3d')

	if type(x[0])==str:
		x= np.array(range(len(x)))
		mplt.xticks(x)

	if type(y[0])==str:
		x= np.array(range(len(y)))
		mplt.yticks(y)

	ax.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=2, antialiased=False)#flag)#coolwarm)cmap=mplt.cm.hot, 
	#ax.set_zlim(0, 1)
	ax.set_xlabel(xlabel)#r'$\phi_\mathrm{real}$')
	ax.set_ylabel(ylabel)#r'$\phi_\mathrm{im}$')
	ax.set_zlabel(zlabel)#r'$V(\phi)$')
	mplt.title(title+'\n')
	return


def plot_3Dscatter(x,y,z, xlabel='', ylabel='', zlabel='', title='', outFile='3Dscattter'+ext, N=1, row=1, col=1, c='r', m='.', dtype=''):
	if N==0:
		fig 	= mplt.figure()
	else:
		fig 	= mplt.gcf()
	ax 		= fig.add_subplot(row, col, N+1, projection='3d')

	if type(x[0])==str:
		x= np.array(range(len(x)))
		mplt.xticks(x)

	if type(y[0])==str:
		x= np.array(range(len(y)))
		mplt.yticks(y)

	ax.set_xlim(0, 8)
	ax.set_ylim(0, np.max(y))

	if dtype=='P':
		ax.set_zlim(0, np.max(z)*100.0)
		ax.zaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))
		z=z * 100.0
		ax.set_zlabel(zlabel+r' x 10$^{-2}$ ')
	
		#ax.set_zscale('log')	
	else:
		ax.set_zlim(np.min(z), np.max(z))
		ax.set_zlabel(zlabel)

	#c, m = ('r', 'o')
	#c, m = ('b', '^')
	ax.scatter(x, y, z, c=c, marker=m, alpha=0.5)#, edgecolor=(1.0,0.0,0.0))

	ax.set_xlabel(xlabel)#r'$\phi_\mathrm{real}$')
	ax.set_ylabel(ylabel)#r'$\phi_\mathrm{im}$')
#r'$V(\phi)$')
	mplt.title(title+'\n')
	return



def checkLoop(struc, start, end):
	if struc[start:end].count('(') == struc[start:end].count(')'):
		isLoop = True
	else:
		isLoop = False
	return isLoop

######################################################

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
	return FreeEnergies, Probabilities


#################################################################################

if __name__=='__main__':
	commands			=	sys.argv[1:]
	main_dir			=	os.getcwd()
	mainPath			=	os.getcwd()+'/'+'e_range_8_complete_298'+'/'
	Vienna_path			=	''
	RNAFold_path		=	''
	RNASubOpt_path		=	''
	RNAplot_path		=	''
	seqFilePath			=	''

######################################################

	seqLens=['145_complete_AT', '146_complete_AT', '147_complete_AT', '148_complete_AT', '149_complete_AT', '150_complete_AT', '151_complete_AT', '152_complete_AT']

	all_P1_counts=[]
	all_AT_counts=[]
	all_freeEnergies=[]

	all_FreeEnergies_BB=[]
	all_FreeEnergies_AB=[]

	all_Probabilities_BB=[]
	all_Probabilities_AB=[]

	for seqLen in seqLens:
		# Read subopt output
		len_path			= seqLen+'/'
		outPath				= mainPath+len_path+'RNASubOpt.out'
		p					= Vienna(mainPath, RNAFold_path, RNASubOpt_path, RNAplot_path, seqFilePath, commands)
		strucs, freeEnergies= p.readOutFile(outPath)
#		#####################################################


		FreeEnergies_BB, Probabilities_BB			= read_probabilities(inFile=mainPath+len_path+'probabilities')
		FreeEnergies_AB, Probabilities_AB			= read_all_substruct_data(inFile=mainPath+len_path+'all_substruct_data')


		######################################################
		P1_len		= 8
		P1_start	= 0	# this is the index
		P1_end		= 118 # This is the index (Aptmer is 118 nucleotides)
	
		AT_len		= 11
		AT_start	= 112 #113:123	139:150	# this is the index
		AT_end		= 150 # This is the index (Aptmer is 118 nucleotides)
	
		j=0
	
		f= open(mainPath+len_path+'P1_AT_counts','w')
		P1_counts=[]
		AT_counts=[]
		for struc in strucs:
			for i in range(P1_len):
				start	= P1_start + i 
				end		= P1_end - i  
				isLoop = checkLoop(struc, start, end)	
				if isLoop == True:
					count1	= struc[start:8].count('(')
					count2	= struc[110:end].count(')')
					if (count1 <= count2):
						P1_count0 = count1
					else:
						P1_count0 = count2
		
					P1_count = 0
					for i in range(len(struc[start:8])):
						if struc[start+i]=='(' and struc[end-1]==')' :
							P1_count = P1_count+1
					break
				else:
					count1, count2, P1_count0, P1_count = 0, 0, 0, 0
		
			for i in range(AT_len):
				start	= AT_start + i 
				end		= AT_end - i  
				isLoop = checkLoop(struc, start, end)	
				if isLoop == True:
					count1	= struc[start:end].count('(')
					count2	= struc[start:end].count('(')
					#count1	= struc[start:start+AT_len].count('(')
					#count2	= struc[AT_end-AT_len:end].count(')')
					if (count1 <= count2):
						AT_count = count1
					else:
						AT_count = count2
	
#					AT_count = 0
#					for i in range(len(struc[start:8])):
#						if struc[start+i]=='(' and struc[end-1]==')' :
#							AT_count = AT_count+1
					break
				else:
					count1, count2, AT_count0, AT_count = 0, 0, 0, 0

			print >> f, j, start, end, P1_count, AT_count
			P1_counts.append(P1_count)
			AT_counts.append(AT_count)
			j=j+1

		f.close()
		P1_counts=np.array(P1_counts, dtype='int')
		AT_counts=np.array(AT_counts, dtype='int')
		freeEnergies= np.array(freeEnergies, dtype='float64')



#

		# Per lengths
	#	plot_3Dscatter(P1_counts,AT_counts,freeEnergies, xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nEnergy  [kcal/mol]', title='Transcription length ='+seqLen[:3], N=0, col=1, row=1)
	#	outFile=mainPath+len_path+'3Dscattter_freeEnergies_'+seqLen[:3]+ext
	#	f=mplt.gcf()
	#	#f.set_size_inches(size_c,size_d)		
	#	mplt.savefig(outFile, dpi=dpi)
	#	mplt.cla()
	#	mplt.clf()
		'''
		plot_3Dscatter(P1_counts,AT_counts,FreeEnergies_BB, xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nEnergy  [kcal/mol]', title='Transcription length ='+seqLen[:3], N=0, col=1, row=1,c='b')
		outFile=mainPath+len_path+'3Dscattter_FreeEnergies_BB_'+seqLen[:3]+ext
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)		
		mplt.savefig(outFile, dpi=dpi)
		mplt.cla()
		mplt.clf()

		plot_3Dscatter(P1_counts,AT_counts,FreeEnergies_AB, xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nEnergy  [kcal/mol]', title='Transcription length ='+seqLen[:3], N=0, col=1, row=1, c='r')
		outFile=mainPath+len_path+'3Dscattter_FreeEnergies_AB_'+seqLen[:3]+ext
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)		
		mplt.savefig(outFile, dpi=dpi)
		mplt.cla()
		mplt.clf()
		'''
		plot_3Dscatter(P1_counts,AT_counts,Probabilities_BB, xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length ='+seqLen[:3], N=0, col=1, row=1, c='b', dtype='P')
		outFile=mainPath+len_path+'3Dscattter_Probabilities_BB_'+seqLen[:3]+ext
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)	
		mplt.savefig(outFile, dpi=dpi)
		mplt.cla()
		mplt.clf()

		plot_3Dscatter(P1_counts,AT_counts,Probabilities_AB, xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length ='+seqLen[:3], N=0, col=1, row=1, c='r', dtype='P')
		outFile=mainPath+len_path+'3Dscattter_Probabilities_AB_'+seqLen[:3]+ext
		f=mplt.gcf()
		f.set_size_inches(size_c,size_d)		
		mplt.savefig(outFile, dpi=dpi)
		mplt.cla()
		mplt.clf()


	#	plot_surf(P1_counts[:10],AT_counts[:10],Probabilities_BB[:10], xlabel='\n\nP1', ylabel='\n\nAT', zlabel='\n\nProbability', outFile=mainPath+len_path+'3Dsurface_'+seqLen[:3]+ext, N=111)

		all_P1_counts.append(P1_counts)
		all_AT_counts.append(AT_counts)
		all_freeEnergies.append(freeEnergies)

		all_FreeEnergies_BB.append(FreeEnergies_BB)
		all_FreeEnergies_AB.append(FreeEnergies_AB)

		all_Probabilities_BB.append(Probabilities_BB)
		all_Probabilities_AB.append(Probabilities_AB)


	#i=0
	#for seqLen in seqLens:
	#	#plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_freeEnergies[i], xlabel='P1', ylabel='AT', zlabel='Energy', title='Transcription length ='+seqLen[:3], N=i, row=rows[i], col=cols[i])
	#	plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_freeEnergies[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nEnergy  [kcal/mol]', title='Transcription length ='+seqLen[:3], N=i, row=3, col=3)
	#	i=i+1
	#outFile=outFile=mainPath+'3Dscattter_all_freeEnergies'+ext
	#f=mplt.gcf()
	#f.set_size_inches(size_a,size_b)
	#mplt.savefig(outFile, dpi=dpi)
	#mplt.cla()
	#mplt.clf()

	'''
	i=0
	for seqLen in seqLens:
		#plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_freeEnergies[i], xlabel='P1', ylabel='AT', zlabel='Energy', title='Transcription length ='+seqLen[:3], N=i, row=rows[i], col=cols[i])
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_FreeEnergies_BB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nEnergy  [kcal/mol]', title='Transcription length ='+seqLen[:3], N=i, row=3, col=3, c='b')
		i=i+1
	outFile=outFile=mainPath+'3Dscattter_all_FreeEnergies_BB'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()


	i=0
	for seqLen in seqLens:
		#plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_freeEnergies[i], xlabel='P1', ylabel='AT', zlabel='Energy', title='Transcription length ='+seqLen[:3], N=i, row=rows[i], col=cols[i])
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_FreeEnergies_AB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nEnergy  [kcal/mol]', title='Transcription length ='+seqLen[:3], N=i, row=3, col=3, c='r')
		i=i+1
	outFile=outFile=mainPath+'3Dscattter_all_FreeEnergies_AB'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()
	'''
	'''
	i=0
	for seqLen in seqLens:
		#plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_freeEnergies[i], xlabel='P1', ylabel='AT', zlabel='Energy', title='Transcription length ='+seqLen[:3], N=i, row=rows[i], col=cols[i])
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_Probabilities_BB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length ='+seqLen[:3], N=i, row=3, col=3, c='b', dtype='P')
		i=i+1
	outFile=outFile=mainPath+'3Dscattter_all_Probabilities_BB'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()


	i=0
	for seqLen in seqLens:
		#plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_freeEnergies[i], xlabel='P1', ylabel='AT', zlabel='Energy', title='Transcription length ='+seqLen[:3], N=i, row=rows[i], col=cols[i])
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_Probabilities_AB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length ='+seqLen[:3], N=i, row=3, col=3, c='r', dtype='P')
		i=i+1
	outFile=outFile=mainPath+'3Dscattter_all_Probabilities_AB'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()
	'''


	i=0
	j=0
	for seqLen in seqLens[:3]:
		#plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_freeEnergies[i], xlabel='P1', ylabel='AT', zlabel='Energy', title='Transcription length ='+seqLen[:3], N=i, row=rows[i], col=cols[i])
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_Probabilities_BB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length BB ='+seqLen[:3], N=j, row=3, col=2, c='b', dtype='P')
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_Probabilities_AB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length AB ='+seqLen[:3], N=j+1, row=3, col=2, c='r', dtype='P')
		print 'i=',i , 'j=',j
		i=i+1
		j=j+2
	outFile=outFile=mainPath+'3Dscattter_all_Probabilities_BB_AB_145_146_147'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()

#	i=0
	j=0
	for seqLen in seqLens[3:6]:
		#plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_freeEnergies[i], xlabel='P1', ylabel='AT', zlabel='Energy', title='Transcription length ='+seqLen[:3], N=i, row=rows[i], col=cols[i])
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_Probabilities_BB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length BB ='+seqLen[:3], N=j, row=3, col=2, c='b', dtype='P')
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_Probabilities_AB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length AB ='+seqLen[:3], N=j+1, row=3, col=2, c='r', dtype='P')
		print 'i=',i , 'j=',j
		i=i+1
		j=j+2
	outFile=outFile=mainPath+'3Dscattter_all_Probabilities_BB_AB_148_149_150'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()

#	i=0
	j=0
	for seqLen in seqLens[6:]:
		#plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_freeEnergies[i], xlabel='P1', ylabel='AT', zlabel='Energy', title='Transcription length ='+seqLen[:3], N=i, row=rows[i], col=cols[i])
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_Probabilities_BB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length BB ='+seqLen[:3], N=j, row=3, col=2, c='b', dtype='P')
		plot_3Dscatter(all_P1_counts[i],all_AT_counts[i],all_Probabilities_AB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length AB ='+seqLen[:3], N=j+1, row=3, col=2, c='r', dtype='P')
		print 'i=',i , 'j=',j
		i=i+1
		j=j+2
	outFile=outFile=mainPath+'3Dscattter_all_Probabilities_BB_AB_151_152'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()





	'''

	for seqLen in seqLens[6:]:
		plot_surf(all_P1_counts[i],all_AT_counts[i],all_Probabilities_BB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length BB ='+seqLen[:3], N=j, row=3, col=2, c='b', dtype='P')
		plot_surf(all_P1_counts[i],all_AT_counts[i],all_Probabilities_AB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nProbability', title='Transcription length AB ='+seqLen[:3], N=j+1, row=3, col=2, c='r', dtype='P')
		print 'i=',i , 'j=',j
		i=i+1
		j=j+2
	outFile=outFile=mainPath+'3Dsurface_all_Probabilities_BB_AB_151_152'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()
	'''

	'''
	i=0
	for seqLen in seqLens:
		x,y=np.meshgrid(all_P1_counts[i][:100], all_AT_counts[i][:100])
		plot_surf(x,y,all_FreeEnergies_BB[i][:100], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nEnergy  [kcal/mol]', title='Transcription length ='+seqLen[:3], N=i, row=3, col=3, c='b')
		i=i+1
	outFile=outFile=mainPath+'3Dsurface_all_Enregies_BB'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()


	i=0
	for seqLen in seqLens:
		x,y=np.meshgrid(all_P1_counts[i][:100], all_AT_counts[i][:100])
		plot_surf(x,y,all_FreeEnergies_AB[i], xlabel='\n\nP1  [bp]', ylabel='\n\nAT  [bp]', zlabel='\n\nEnergy  [kcal/mol]', title='Transcription length ='+seqLen[:3], N=i, row=3, col=3, c='r')
		i=i+1
	outFile=outFile=mainPath+'3Dsurface_all_Energies_AB'+ext
	f=mplt.gcf()
	f.set_size_inches(size_a,size_b)
	mplt.savefig(outFile, dpi=dpi)
	mplt.cla()
	mplt.clf()
	'''

# Per Energy
#plot_3Dscatter(P1_counts,AT_counts,z, xlabel='', ylabel='', zlabel='', outFile='3Dscattter'++ext, N=111)

