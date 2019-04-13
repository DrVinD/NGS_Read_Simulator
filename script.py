##  ---------------------------NGS Read Simulator---------------------------  ##
##  Randomly picks out reads from a reference genomic fasta file and outputs  ##
##  them as a fastq format.                                                   ##
##  Number of Reads = 100,000                                                 ##
##  Read length     = 50bp                                                    ##
##  Error Rate      = 0.01                                                    ##
##  (1% of the time a base is wrongly called by the sequencer)                ##
##  Author  :   Vinayak Dev                                                   ##
##  Date    :   16-03-2019                                                    ##
# coding: utf-8


import time
import os
import random

print("----------------------------------------------------------\n")
print("******************* NGS READ SIMULATOR *******************\n")
print("----------------------------------------------------------\n\n")
fname = input("Enter the Reference Genome(.fa) filename(e.g. hg38.fa): ")
start = time.time()

###--------------------------------------------------------------------------###
## Cleaning & Processing the input(reference) genome in order to pick         ##
## 100,000 reads of 50 bp at random.                                          ##
###--------------------------------------------------------------------------###
s = ""
try:
    with open(fname) as f:
        for lines in f:
            line = str(lines).strip('\n').upper()
            if('>' not in line):  # excluding fasta-format seq identifier line
                s = s+str(line)   # formulated genomic seq in string var s
except IOError:
    print("Could not read file")

s = s.replace('N', '')            # Clean and processed genomic seq in variabl s
print("\n----------Reference genome was cleaned and processed.----------\n")


###--------------------------------------------------------------------------###
## Picking 100,000 50bp-reads at random from processed genomic seq(var 's')   ##
###--------------------------------------------------------------------------###
l = len(s)
print("The size of the input Reference Genome: ",fname, "is",l, "bp.\n")
seq = ""
Reads = []
for i in range(0, 100000):
    x = random.randint(0, l-50)   # 100,000 random position numbers
    tmp = s[x : x+50]             # temporary variable to pick a read of 50bp
    Reads.append(tmp)
    seq+=str(tmp)                 # appending read-seq to a string var 'seq'

print("----------100,000 50bp-Reads picked randomly----------")


###--------------------------------------------------------------------------###
## Putting error rate of 0.01 at base calling                                 ##
###--------------------------------------------------------------------------###
count = 0
flag = []                         #   A/Q 1% is base call error of sequencer #
NT_A = ["T", "G", "C"]            #   error rate = 0.01;  #reads = 100,000;  #
NT_T = ["A", "G", "C"]            #         read length = 50bp , therefore   #
NT_G = ["A", "T", "C"]            #   0.01*50*100,000 = 50,000th bp, probably#
NT_C = ["A", "T", "G"]            #    be sequenced wrongly, to simulate the #
                                  #               sequencer error.           #
while(count < 50000):
    i = random.randint(0, len(seq))# Random number for error position
    if i in flag:
        continue
    else:
        flag.append(i)            # position saved in flag to remove redundancy
        #print(i, seq[i-2:i+2])
        tmp_str = list(seq)
        if tmp_str[i] == 'T':
            tmp_str[i]= NT_T[random.randint(0, 2)]          # Randomly picking
        elif tmp_str[i] == 'A':                             # the nucleotide.
            tmp_str[i]= NT_A[random.randint(0, 2)]
        elif tmp_str[i] == 'G':
            tmp_str[i]= NT_G[random.randint(0, 2)]
        else:
            tmp_str[i]= NT_C[random.randint(0, 2)]
        seq = ''.join(tmp_str)
        count+=1


FINAL_Reads = []
l = 50
i = 0
while i < len(seq):
    tmp = seq[i:i+l]
    i=i+l
    FINAL_Reads.append(tmp)       # Final Reads list

print("\n----------Error integrated succesfully----------\n")


###--------------------------------------------------------------------------###
## Converting reads to fastq format with dummy quality values generated       ##
## according to Illumina Sequencing technology.                               ##
###--------------------------------------------------------------------------###


QS = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"

f = open('Simulated_Reads.fq', 'w')
j=0.000001                        # Read identifier
for r in FINAL_Reads:
    Q=""
    for i in range(0, 50):
        x = random.randint(0, 40)
        Q+=QS[x]                  # Providing Dummy quality values randomly
    f.write("@"+"Illumina Read:"+str(j)+'\n'+str(r)+'\n'+"+"+'\n'+str(Q)+'\n\n')
    j+=0.000001
f.close()
end = time.time()

print("----------Simulated_Reads.fq Fastq file generated succesfully----------")
print("\nTotal Simulation Time in seconds:", end - start)

## Last: hg38 whole genome simulation time on 3.1GHz:i5:macOS:8GBRam was
## 3159.192456960678s = 52.65 mins.
