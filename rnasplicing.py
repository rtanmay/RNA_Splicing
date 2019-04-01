# This program only solves for non overlapping introns and only for certain combinations
# E.g: removal of ana from banana. bna or ban (two combinations)
############################################################################################
# Importing the libraries
from Bio.Seq import Seq
from Bio import SeqIO

############################################################################################
# Function to remove the introns in prerna
def remove_intron(prernastr, intronstr):
	index = 0
	length = len(intronstr)
	while prernastr.find(intronstr) != -1:
		index = prernastr.find(intronstr)
		prernastr = prernastr[0:index] + prernastr[index+length:]
	return prernastr

############################################################################################
"""
Input is in FASTA format
First string is the DNA String
Rest below them are the introns string
"""
counter=0
intron_list=[]
for record in SeqIO.parse("input.txt","fasta"):
	if (counter==0):
		dna_string=record.seq # record.seq returns an Seq object
		counter+=1
	else:
		intron_list.append(str(record.seq))
print("Reading Done!")
countt=0
print("DNA =",dna_string)
for i in intron_list:
	countt=countt+1
	print("Intron %d = %s" % (countt,i))
print("")
prerna=dna_string
prerna=str(prerna)

############################################################################################
"""
Removing the introns
In case of overlapping introns in the prerna, the first occurence of intron is
removed and then the subsequent string is considered for removal
"""
print("Removing the introns")
for i in intron_list:
	prerna=remove_intron(prerna,i) 
print("Removal Done!")
prerna=Seq(prerna)
############################################################################################
# Transcribtion
print("Doing Transcribtion")
prerna=prerna.transcribe()
print("Transcribtion Done!")
maturerna=prerna
print("MaturedRNA = ",maturerna)
print("MaturedRNA length = ",len(maturerna))
print("")
############################################################################################
# Converting the matured rna to protein
print("Translating the matured RNA")
protein=maturerna.translate(to_stop=True)
print("Translation Done!")
############################################################################################
# Output the protein formed
print("Protein=",protein)
############################################################################################