import numpy as np
#import pandas as pd
import re
from string import printable
import time


def read_dna():
    gen="" #put all dna seriers as one string
    file1=open('DNA.fq','r')
    while True:
        line=file1.readline() #read line by line
       # print(line)
        if len(line)==0:
           break;
        #if len(line)==0: #if reach to the end then break out of loop
          #  break;
       # if line.startswith('A') or line.startswith('T') or line.startswith('G') or  line.startswith('C'):
        regex = re.compile('[+@_!#$%^&*()<>?/\|}{~:]')
        if regex.search(line) is not None:
                {}
        else:
             gen+=line


    return gen


# d is the number of characters in the input alphabet
d = 256
q=101


def compute_hash(dna,pat):
    dna_hash=0
    pat_hash=0
    for i in range(len(pat)):

        pat_hash=(d*pat_hash+ord(pat[i]))%q
        dna_hash=(d*dna_hash+ord(dna[i]))%q
    #print("hashhhhh are",pat_hash,dna_hash)
    return dna_hash,pat_hash

def Rabin_karp(dna,pat):
    #q=3 #prime number using in computation
    dna_len=len(dna)#dna length
    pat_len=len(pat)#pattern length
    match_index=[]
    h=1
    for i in range(pat_len - 1):
        h = (h * d) % q
        #compute the initial hash for p and dan

    dna_hash,pat_hash=compute_hash(dna,pat) # compute the first value of hash for pattern and first window of dna
    start_time = time.time()
    for s in range(dna_len-pat_len+1):   #matching
       # print(dna_hash, pat_hash)
        if dna_hash==pat_hash: # if they are equal then compare char by char
            #print("equalllllll")

            for j in range(pat_len):

                if pat[j]!=dna[j+s]:
                 # print(pat[j],dna[j+s])
                  break
            #print("m is", m)
            if(j==pat_len-1): # this will be equal in case of full match
                 match_index.append(s)

        if s<dna_len-pat_len: # we want to continue our comparsion by taking the next window of DNA
           # print (dna[s])
            dna_hash = (d*(dna_hash-ord(dna[s])*h) + ord(dna[s+pat_len]))%q#new hash for next window by subtract MSB and add LSB
           # print("hash computed")
    end_time = time.time()
    print("time by rabin  =",end_time-start_time)
    return match_index


def kmp(dna, pat):


    #  prefix suffix array
    match_index=[]
    prefix_suffix = [0] * len(pat)
    pat_index = 0  # index for pat[]

    # Preprocess the pattern (calculate lps[] array)
    start_time = time.time()
    compute_prefix_table(pat, prefix_suffix )
   # print(lps)

    i = 0  # index for txt[]

    while i < len(dna):
        if pat[pat_index] == dna[i]:
            i += 1
            pat_index += 1

        if pat_index == len(pat):
            match_index.append(i - len(pat))
            pat_index = prefix_suffix [pat_index - 1]

            # mismatch after j matches
        else:
            if pat[pat_index] != dna[i]:

             if pat_index == 0:
                i += 1
             else:
                 pat_index = prefix_suffix[pat_index - 1]#get how many steps we will skip
    end_time = time.time()
    print("Time by kmp =",end_time-start_time)
    return match_index
def compute_prefix_table(pat, prefix_suffix):
    prefix_suffix[0]=0 # initilize the first index to zero
    length=0
    index=1  # start matching between 0,1
   # M=len("amera")
    while index <(len(pat)):
        if pat[length]==pat[index]: # suffix and prefix are match

             length+=1 #increase the length of matching by 1
             prefix_suffix[index]=length # the length of matching in index
             index += 1  # move to next charecter in pattern
        if pat[length]!=pat[index]: # in case of mismatch then we have more than one option
            if(length==0):
                prefix_suffix[index]=length
                index+=1 #move on to next chatecter in pattern
            else :
                length=prefix_suffix[length-1] # in case of mismatching check the last matching charecter length to be able to skip as possible as we can





#try to use the memorization idea to compute lps array on demand

def main():
    gen=read_dna()
   # print(gen)

    # Your statements here
    pattern="GGCATATGAAAATTTATTACTACAGTGTTTT"
    matches=Rabin_karp(gen, pattern)
    print("Rabin run and pattern is found at positions ")
    for match in matches:
        print(match," "),



    matches= kmp(gen,pattern)
    print("KMP run and pattern is found at positions ")
    for match in matches:
        print(match," "),



if __name__ == '__main__':
    main()