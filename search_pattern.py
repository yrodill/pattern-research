#!/usr/bin/python
#-*- coding: utf-8 -*-


import sys
import argparse

"""
Bothorel Benoît
January 2019
Launch the script with :
python3 search_pattern.py -algo [str] -pattern [str] -file [.txt or .fasta]         
"""
from copy import copy,deepcopy


"""
Parse a fasta/txt file format
Returns a string containing the sequence
"""
def parse_file(file):
    string =""
    f=open(file,'r')
    seqs = f.readlines()
    
    for i in range(len(seqs)):
        if seqs[i][0] != ">":
            temp=""
            for j in range(len(seqs[i])):
                if seqs[i][j] != "\n":
                    temp += seqs[i][j]
            string += temp

    return string

"""
Function used to make the search
according to the arguments passed in entry
"""
def search(algo,pattern,seq):
    if algo.lower() == "naive":
        naive(pattern,seq)
    elif algo.lower() == "mp":
        p=MP_precalcul(pattern)
        MP_search(pattern,seq,p)
    elif algo.lower() == "bm":
        BM_search(pattern,seq)
    else:
        print ("\nWrong algo !")

"""
Naive search : first implementation without a while loop
Finds the motive match(s) in the sequence
"""
def naive_search(pattern,seq):
    count=0
    pos=[]
    for i in range(0,len(seq)-len(pattern)+1):
        if pattern == seq[i:i+len(pattern)]:
            count+=1
            pos.append(i)
    for index in pos:
        print ("Match at : {}".format(index))
    print ("Total: {}".format(count))

"""
Implementation of the naive search according to the pseudocode given
Finds the patterns in the sequence and display the index
"""
def naive(pattern,seq):
    count=0
    for i in range(0,len(seq)-len(pattern)+1):
        j=0
        while j<=len(pattern)-1 and seq[i+j] == pattern[j]:
            j+=1
        if j == len(pattern):
            print ("Match at : {}".format(i))
            count+=1
    print ("Total: {}".format(count))

"""
Determine the prefixes of a pattern
Returns a list of integers
"""
def MP_precalcul(pattern):
    prefixes=[]
    prefixes.append(-1)

    j=0
    for i in range(0,len(pattern)-1):
        while j > -1 and pattern[i] != pattern[j+1]:
            j=prefixes[j]
        j+=1
        prefixes.append(j)
    return prefixes

"""
Morris-Pratt algorithm to find pattern match(s)
"""
def MP_search(pattern,sequence,prefixes):
    count=0
    j=0
    for i in range(0,len(sequence)):
        while j>-1 and seq[i] != pattern[j]:
            j=prefixes[j]
        j+=1
        if j == len(pattern):
            print ("Match at : {}".format(i-len(pattern)+1))
            j=prefixes[j-1]
            count+=1
    print ("Total: {}".format(count))

"""
Creation of the alphabet (useful if you want to use
proteinic sequences)
Returns a list of chars
"""
def getAlphabet(sequence):
    alphabet=[]
    for letter in sequence:
        if letter not in alphabet:
            alphabet.append(letter)
    return alphabet

"""
Returns the table of the wrong characters
"""
def badCharTable(pattern,sequence):
    table={}
    alphabet=getAlphabet(sequence)
    for l in alphabet:
        if l not in pattern:
            table[l]=len(pattern)
        else:
            for i in range(len(pattern)):
                if pattern[i]==l and l not in table.keys():
                    for j in range(len(pattern)):
                        if j>i and pattern[j] == l:
                            table[l]=j-i-1
                            break
                        else:
                            table[l]=len(pattern)-i-1
                 
    return table

"""
Returns the suffixes' table of the pattern
"""
def getSuffixesTable(pattern):
    m=len(pattern)
    suffixe = [0] * m
    suffixe[m - 1] = m
    g = m - 1
    previous_i = 0 #sauvegarde de l'indice i du tour de boucle précédent
    for i in range(m - 2, -1, -1):     
        if (i > g and suffixe[i + m - 1 - previous_i] < i - g): #la valeur est deja dans la table
            suffixe[i] = suffixe[i + m - 1 - previous_i]
        else:
            if i < g:
                g = i
            previous_i = i
            while (g >= 0 and pattern[g] == pattern[g + m - 1 - previous_i]):
                g -= 1 #g est decremente tant que les lettres sont identiques 2 a 2
            suffixe[i] = previous_i - g 
    return suffixe

"""
Returns the good suffixes table
"""
def getGoodSuffixes(pattern):
    m=len(pattern)
    suffixe = getSuffixesTable(pattern)
    bmGs = [0] * m
    for i in range(m):
        bmGs[i] = m
    for i in range(m-1, -1, -1):
        if suffixe[i] == i+1: #si il existe un suffixe dans le pattern
            for j in range(m - 1 - i):
                if (bmGs[j] == m):
                    bmGs[j] = m - 1 - i #la valeur de saut est changee pour correspondre au suffixe
    for i in range(m - 1):
        bmGs[m - 1 - suffixe[i]] = m - 1 - i #utilisation de la table de suffixe
    return bmGs

"""
Boyer and Moore algorithm
"""
def BM_search(pattern,sequence):
    count=0
    j=0
    i=0
    bmBc=badCharTable(pattern,sequence)
    bmGs=getGoodSuffixes(pattern)
    while i <= len(sequence) - len(pattern):
        j=len(pattern)-1
        while j>=0 and sequence[i+j] == pattern[j]:
            j=j-1
        if j<0 :
            print ("Match at : {}".format(i))
            count+=1
            i = i+bmGs[0]
        else:
            i=i+max(bmGs[j],(bmBc[sequence[i+j]]-len(pattern)+1+j))

    print ("Total : {}".format(count))

"""
Main function to define the arguments needed
and to run the program
""" 
if __name__=="__main__":

    parser=argparse.ArgumentParser()
    parser.add_argument('-algo',help="Choose the algorithm to process your sequence from file. Possible value :'naive' or 'mp'")
    parser.add_argument('-pattern',help="Enter a string representing the pattern that you will check. Example : 'ATG'")
    parser.add_argument('-file',help="Choose your fasta or txt file to be processed")
    args=parser.parse_args()

    # init
    if len(sys.argv)<3:
        print ("Wrong input files.\nUsage : -algo [str] -pattern [str] -file [.txt or .fasta]")               
    else:
        algo=args.algo
        pattern=args.pattern.upper()
        seq=parse_file(args.file)
        search(algo,pattern,seq)