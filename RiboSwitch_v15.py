import os
import math
import sys
import shutil
import numpy as np
import subprocess
import shlex
import matplotlib.pyplot as mplt
import cdecimal
#from Bio import SeqIO
import time
time_start = time.time()

# Setting numbers percisions upto 100 decima places
myContext = cdecimal.Context(prec=100, rounding='ROUND_HALF_DOWN')
cdecimal.setcontext(myContext)
D= cdecimal.Decimal


class Vienna:
	''' A class for all Vienna package related function '''

	def __init__(self, mainPath, RNAFold_path, RNASubOpt_path, RNAplot_path, seqFilePath, commands, T=298.15, energy_range=8):
		self.mainPath			= mainPath
		self.seqFilePath		= seqFilePath		# Sequence file path
		self.RNAFold_path		= RNAFold_path		# RNAfold program path
		self.RNASubOpt_path		= RNASubOpt_path	# RNASubopt program path
		self.RNAplot_path		= RNAplot_path
		self.T					= D(T)				# Folding Temperature in Kelvin
		self.k					= D(0.0019872041)	#(18) # Boltzmann constant in KCal/mol  
		self.energy_range		= energy_range		# Energy range for which suboptima structures are searched by RNAfold
		self.K_C_conv			= D(273.15)			# conversion factor from Klevin to Celsius 
		return


	def readFastaSeqs0(self, seqFilePath):
		'''Reading a set of sequences in fassta format using BioPython module'''
		records		= SeqIO.parse(FilePath, 'fasta')
		seqs=[]
		seq_ids=[]
		for record in records:
			seq_id	= record.id
			seq		= record.seq
			seq_ids.append(seq_id)
			seqs.append(seq)
		return seq_ids, seqs


	def readFastaSeqs(self, seqFilePath):
		'''Reading a single sequence in Fasta format'''
		f=open(seqFilePath)
		data = f.read()
		seq_id, seq, i= data.split('\n')
		print seq_id[1:], seq
		return seq_id[1:], seq


	def RNAFold(self, seq, outPath):
		''' Folds a given RNA sequence using RNAFold program'''
		''' Command: RNAfold -p -d2 --noLP < test_sequenc.fa > test_sequenc.out '''
		''' Command: RNAfold    -d2 --noLP < test_sequenc.fa > test_sequenc.out '''

		print 'Calculating optimal structure of the sequence ..', seq
		stcommand	= [self.RNAFold_path,
						'-p',
						'-d2',
						'--noLP']
		s			= subprocess.Popen(stcommand, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		stdoutdata, stderrdata= s.communicate(input=seq)
		s.wait()
		f=open(outPath, 'w')
		print outPath

		stdoutdata= str(stdoutdata).split('\n')
		stderrdata= str(stderrdata).split('\n')
		for i in range(len(stdoutdata)):
			f.write(stdoutdata[i])
		for i in range(len(stderrdata)):
			f.write(stderrdata[i])
		f.close()
		return

	def RNASubOpt(self,seq,outPath, Enrgy_sorted='True', dos='False'):
		'''Calculates suboptimal structures via RNAsubopt, using the following options: 
			#-e, --deltaEnergy
			#-s, --sorted
			#-T, --temp=DOUBLE
			#-D, --dos 'Compute density of states instead of secondary structures (default=off) :This option enables the evaluation of the number of secondary 				structures in certain energy bands arround the MFE.'
		'''
		print 'Calculating suboptimal structures of the sequence ..', seq
		stcommand	= [self.RNASubOpt_path,
						'--deltaEnergy',str(self.energy_range),
						'-s',
						'-T', str(self.T - self.K_C_conv),
						'--noLP']
		s			= subprocess.Popen(stcommand, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdoutdata, stderrdata= s.communicate(input=seq)
		s.wait()
		f=open(outPath, 'w')
		print outPath

		stdoutdata= str(stdoutdata).split('\n')
		stderrdata= str(stderrdata).split('\n')
		for i in range(len(stdoutdata)):
			f.write(stdoutdata[i]+'\n') 		# *** changed in linux version
		for i in range(len(stderrdata)):
			f.write(stderrdata[i])
		f.close()
		return


	def readOutFile(self, outFile):
		'''Reads the output file of the RNAsubopt program'''
		print 'Reading Subopt output file ..'
		f		= open(outFile,'r')
		data	= f.read()
		data	= data.split('\n')				# *** changed in linux version
		print data[0]
		seq, MFE, length	= data[0].split()
		structures	= data[1:]
		energies	= []
		strucs		= []
		i=0
		#if structures[-1].split(' ') == []:
		#	structures=structures[:-1]
		#elif structures[-1].split(' ') == '\n':	# *** changed in linux version
		#	structures=structures[:-1]			# *** changed in linux version
		#elif structures[-1].split(' ') == '':		# *** changed in linux version
		#	structures=structures[:-1]			# *** changed in linux version

		print len(structures), 'suboptimal structures are found.'
		print 'last structures=', structures[-2]==''
		structures	= structures[:]
		for struc in structures:
			struc			= struc.strip('\n') # *** changed in linux version
			if structures[i] != '':
#				print len(struc.split(' ')),',',i+1,',',struc
				struc, energy	= struc.split(' ')	# *** changed in linux version
				energy			= float(energy)
				strucs.append(struc)
				energies.append(D(energy))
			i=i+1
		energies	= np.array(energies)
		return strucs, energies


	def mplot(self, x, y, label_x= 'x', label_y='y', outFile='FreeEnergy_versus_P.png'):
		'''A plotting fuction for suing matPlotLib'''
		print 'Creating plots of free energies and probabilities ..'
		mplt.plot(x,y, 'r.')
		#mplt.plot(x,y, 'b.')		
		mplt.xlabel(label_x)
		mplt.ylabel(label_y)
		mplt.savefig(outFile)
		mplt.cla()
		mplt.clf()		
		return


	def RNAplot(self, seq, strucs,format='svg', N=100):
		'''Utlizes RNAplot to draw calculated structure'''
		print 'Plotting RNA structure ..'
		stcommand	= [self.RNAplot_path,
					'-o', format]
		s			= subprocess.Popen(stcommand, stdin=subprocess.PIPE)
		for i in range(N):
			try:
				s.stdin.write(">"+str(i+1)+"\n")
				s.stdin.write(seq+"\n")
				s.stdin.write(strucs[i]+"\n")
			except:
				break
		stdoutdata, stderrdata= s.communicate()
		s.wait()
		return


	def Entropy(self, strucs, A=1.0):
		'''Calculates entropy from the number of unpaired bases'''
		print 'Calculating entropies ..'
		entropies	= []
		for struc in strucs:
			entropy	= D(3.0/2.0)*D(A)*D(math.log(struc.count('.')))
			entropies.append(D(entropy))
		return entropies


	def FreeEnergy(self, entropies, k, a=0.0, b=1.0, c=1.0):
		'''Calculates from entropy
			# entropy is the number of unpaired bases
			# a, b, c are contants
			# k is a loop in the structure'''
		print 'Calculating Free Energies  ..'
		freeEnergies	= D(a)+D(b)( D(k-1.0))+(D(c) * entropies)
		return freeEnergies


	def PartitionFunction(self, freeEnergies):
		'''Calculates the partition function from free energies'''
		print 'Calculating the partition function ..'
		qs= np.exp((-freeEnergies)* (D(1.0)/(self.T * self.k)) )
		#print '%E'%qs
		Q	= np.sum(qs)
		print 'Q =','%G'%Q
		return Q


	def Probability(self, freeEnergies, Q, fileName='probabilities'):
		'''Calculates probablities from free energy(ies) and the partition function'''
		print 'Calculating probablities ..'
		probArray	= (D(1.0)/Q) *  np.exp((-freeEnergies)/D(self.T * self.k))
		f= open(fileName,'w')
		print >> f, 'Structure', 'Free Energy','exp(-E/T)', 'Probability'
		for i in range(len(probArray)):
			print >> f, i, '%G'%freeEnergies[i], '%G'%np.exp(-freeEnergies[i]/self.T), '%G'%probArray[i]
		f.close()
		#average_freeEnergies= np.sum(freeEnergies)/len(freeEnergies)
		#delta_freeEnergies= freeEnergies- average_freeEnergies
		#probArray_diff	= (D(1.0)/Q) * (  np.exp((-freeEnergies/self.T) * self.k) -  np.exp((-average_freeEnergies/self.T) * self.k)  )
		#probArray_delta	= (D(1.0)/Q) *  np.exp((-delta_freeEnergies/self.T) * self.k)
		#print 'Probabilities differences from average =', probArray_diff
		#print 'Probabilities delta free energies = =', probArray_diff
		return probArray

	def FreeEnergy_from_struc(self, seq, struc, baseParingParams):
		# split structure into paired and unpaired
		# zero ernergy for unpaired
		# calc from brackets, based on sequence
		return

	def subStructure_search(self, strucs, substruc):
		indecies=[]
		for struc in strucs:
			index= struc.find(substruc)
			indecies.append(index)
		return np.array(indecies)



################################################################################
### Main Program:

if __name__=='__main__':
	commands			=	sys.argv[1:]
	main_dir			=	os.getcwd()
	mainPath			=	os.getcwd()+'/'#'/' instead for linux					# *** changed in linux version
	Vienna_path			=	'/home/osama/Programs/viennarna_2.2.7-1_amd64/usr/bin/'	# *** changed in linux version
	RNAFold_path		=	Vienna_path+'RNAfold'									# *** changed in linux version
	RNASubOpt_path		=	Vienna_path+'RNAsubopt'									# *** changed in linux version
	RNAplot_path		=	Vienna_path+'RNAplot'									# *** changed in linux version
	#seqFilePath			=	'2gis.fasta'
	#seqFilePath			=	'SAMI_short_trial.fasta'#'SAMI.fasta'
	seqFilePath			=	'SAMI.fasta'
	logFileName			=	'report.log'
	logFile				= open(logFileName,'w')

	# Read sequence and initiate class
	p			= Vienna(mainPath, RNAFold_path, RNASubOpt_path, RNAplot_path, seqFilePath, commands)
	seq_id, seq	= p.readFastaSeqs(seqFilePath)

	# RNAFold: finds the optimal fold
	inPath	=	seqFilePath
	outPath	=	mainPath+'RNAFold.out' 
	p.RNAFold(seq, outPath)

	# Get all suboptimal structures
	inPath	=	seqFilePath
	outPath	=	mainPath+'RNASubOpt.out'
	p.RNASubOpt(seq, outPath)

	# Read subopt output
	outPath	=	mainPath+'RNASubOpt.out'
	strucs, freeEnergies= p.readOutFile(outPath)

	# Calculate Entropies
	entropies	=	p.Entropy(strucs)

	# Calculates free energies using entropies
	#	freeEnergies = FreeEnergy(self, entropies, k, a=0.0, b=1.0, c=1.0):

	# Partition function
	Q			=	p.PartitionFunction(freeEnergies)

	# Porbabilities
	probArray	=	p.Probability(freeEnergies, Q)

	# Writeout data in a table
	print >> logFile, 'Partition function Q = ', '%E'%Q
	print >> logFile, '\n'
	print >> logFile, 'Suboptimal structures data:'
	print >> logFile, 'No.','\t','Free Energy','\t','Probability','\t','Entropy'
	for i in range(len(freeEnergies)):
		print >> logFile, i,'\t','%G'%freeEnergies[i],'\t','%G'%probArray[i],'\t', '%G'%entropies[i]
	print >> logFile, 'Sum of All Probabilities=', '%E'%np.sum(probArray)
	print >> logFile, '\n'
	
	# Plot free energies and prbablities
	p.mplot(freeEnergies, probArray, label_x='Free Energy', label_y='Probability', outFile= 'FreeEnergy_versus_P.png')
	p.mplot(probArray, freeEnergies, label_y='Free Energy', label_x='Probability', outFile= 'P_versus_FreeEnergy.png')
#	p.mplot(freeEnergies, p2, label_x= 'Free Energy', label_y='Probability', outFile= 'FreeEnergy_versus_expET.png')	

	logFile.close()

	os.mkdir('RNAPlot_output')
	os.chdir('RNAPlot_output')
	# Plotting RNA structures
	p.RNAplot(seq, strucs,format='svg', N=len(strucs))


####################################################################


	time_end = time.time()
	print '\n'
	print 'Time of excultion: ', time_end - time_start, 'seconds'

