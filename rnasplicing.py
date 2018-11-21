from Bio.Seq import Seq
from Bio import SeqIO

#String to store the maturerna
maturerna=""

#Function to remove the introns in prerna
def remove_intron(mainstring,substring):
	global maturerna
	maturerna=""
	mainstring=str(mainstring)
	substring=str(substring)
	mlength=len(mainstring)
	slength=len(substring)
	i=0
	while i<mlength:
		initialpos=i
		j=0
		while i<mlength and j<slength and substring[j]==mainstring[i]:
			i+=1
			j+=1
		j-=1
		i-=1
		if j==slength-1:
			i=initialpos+j+1
		else:
			maturerna+=mainstring[i]
			i=initialpos+1

"""
The input is in FASTA format
The first string is the DNA String
The rest below them are introns string
"""
counter=0
lis=[]
for record in SeqIO.parse("input.txt","fasta"):
	if counter==0:
		dna_string=record
		counter+=1
	else:
		lis.append(str(record))

"""
Remove this when above starts working 

#Input the DNA string
dna_string=Seq(input())

#Input the number of introns
num_introns=int(input())
i=0
lis=[]
while i<num_introns:
	intron_string=input()
	lis.append(intron_string)
	i+=1
"""


#Get the rna string from dna
prerna=dna_string.transcribe()
prerna=str(prerna)

"""
Removing the introns
In case of overlapping introns in the prerna, the first occurence of intron is
removed and then the subsequent string is considered for removal
"""
for i in lis:
	remove_intron(prerna,i) 
	prerna=maturerna
	maturerna=""
maturerna=Seq(prerna)

#Converting the matured rna to protein
protein=str(maturerna.translate())

#Output the protein formed
print(protein)