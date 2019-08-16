#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 09:28:56 2018

@author: Hyeyeon Hwang (The LaSalle Lab at University of California, Davis)

Version 3
Updated August 16, 2019

"""
# This script converts Bismark format files to Permeth and DSS format files.

# Input: Bismark format:
# inData[0..5] = [chromosome] [start position] [end position] [% methylated] 
#               [count methylated] [count unmethylated]

# Output: Permeth format:
# outData[0..8] = [chromosome] [start position] [end position] [% methylated-sumReads] 
#               [placeholder00] [strand] [placeholder01] [placeholder02] [color] 

# Output: DSS format:
# outData[0..3] = [chr = chromosome] [pos = start position] 
#               [N = count methylated + count unmethylated] [X = count methylated]

import csv
import errno
import gzip
import io
import multiprocessing as mp
import os
import sys
import time

# write header for each Permeth file
def writeHeader_permeth(ph_chrNum):
    with open(outPermDir + "PerMeth_" + trackName+ "_chr" + ph_chrNum + ".bed","w") as permOutHeader:
        writer = csv.writer( permOutHeader, delimiter = '\t', lineterminator = '\n')
        
        headerData = [None] * 9
        headerData[0] = "track name=" + trackName + "chr" + ph_chrNum+ \
            " description=" + trackName + "chr" + ph_chrNum + \
            " useScore=0 itemRgb=On db=" + genome
        
        writer.writerow(headerData)

# write header for each DSS file
def writeHeader_dss(dh_chrNum):
    with open(outDssDir + trackName+ "_chr" + dh_chrNum + ".DSS.txt","w") as dssOutHeader:
        writer = csv.writer( dssOutHeader, delimiter = '\t', lineterminator = '\n')
        
        headerData = [None] * 4
        headerData[0] = 'chr'
        headerData[1] = 'pos'
        headerData[2] = 'N'
        headerData[3] = 'X'
        
        writer.writerow(headerData)

# write data rows in Permeth format               
def writeData_permeth(pd_chrNum):
    with open(outPermDir + "PerMeth_" + trackName + "_chr" + pd_chrNum + ".bed","a") as permOutFile:
        writer = csv.writer( permOutFile, delimiter = '\t', lineterminator = '\n')
    
        for r in chrms[pd_chrNum]:
            inData = r.decode().split("\t")
            
            inData[3] = int(round(float(inData[3])))
            percMeth = float(inData[3]) / 100 
            inData[2] = int(inData[2])
            sumReads = int(inData[4]) + int(inData[5])
            inData[3] = str(inData[3]) + "-" + str(sumReads)
              
            # RGB color depends on percent methylated 
            rgb = "0,0,0"         # black
            if percMeth > 0 and percMeth <= 0.6:
                rgb = "27,74,210" # blue
            elif percMeth > 0.6 and percMeth <= 0.8:
                rgb = "27,210,57" # green
            elif percMeth > 0.8:
                rgb = "210,27,27" # red
            
            outData = [None] * 9
            outData[0] = inData[0]  # chromosome
            outData[1] = inData[1]  # start position
            outData[2] = inData[2]  # end position
            outData[3] = inData[3]  # % methylated
            outData[4] = 0          # placeholder00
            outData[5] = "+"        # strand 
            outData[6] = 0          # placeholder 01
            outData[7] = 0          # placeholder 02
            outData[8] = rgb        # color
            
            writer.writerow(outData)

# write data rows in DSS format
def writeData_dss(dd_chrNum):

    with open(outDssDir + trackName + "_chr" + dd_chrNum + ".DSS.txt","a") as dssOutFile:
        writer = csv.writer( dssOutFile, delimiter = '\t', lineterminator = '\n')
    
        for r in chrms[dd_chrNum]:
            inData = r.decode().split("\t")
            sumReads = int(inData[4]) + int(inData[5])      

            outData = [None] * 4
            outData[0] = inData[0] # chr = chromosome
            outData[1] = inData[1] # pos = start position
            outData[2] = sumReads  # N = count methylated + count unmethylated
            outData[3] = inData[4] # X = count methylated
            
            writer.writerow(outData)

#*********************** main **************************
startTime = time.time()

inFile = sys.argv[1]
genome = sys.argv[2]
cores = int(sys.argv[3])
trackName = inFile.split("_")[0]
wdir = os.getcwd() 

# create Permeth directory if it doesn't exist
outPermDir = wdir + "/Permeth/"  
try:
    os.makedirs(outPermDir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

# create DSS directory if it doesn't exist 
outDssDir = wdir + "/DSS/"   
try:
    os.makedirs(outDssDir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

# if sample files are in different directories
#   if genome == 'rheMac8':
#       genomeDir = 'Genome/rheMac8/'
#   elif genome == 'hg38':
#       genomeDir = 'Genome/hg38/
#   inDir = wdir + genomeDir   
#   with gzip.open(inDir + inFile, "r") as bisInFile:

with gzip.open(inFile, mode="r") as bisInFile:
    reader = io.BufferedReader(bisInFile)

    # validChrm = chromosomes of interest = [ 1..22 or 1..20, X, Y, M ] 
    validChrm = []
    
    if genome == "hg38" or genome == "hg19":
        print("genome in hg38")
        for i in range(1, 23): # 1..22
            validChrm.append(str(i))
    if genome == "rheMac8" or genome == "rn6":
        for i in range(1, 21): # 1..20
            validChrm.append(str(i))
    if genome == "mm9" or genome == "mm10":
        for i in range(1, 20): # 1..19
            validChrm.append(str(i))
            
    validChrm.append("X") 
    validChrm.append("Y")
    validChrm.append("M")

    # chrms is a dictionary with keys = each valid chromosome number
    # chrms[key = chrmNum] contains all the rows with chrmNum in Bismark file
    chrms = {}
    # chrmNumList is list of unique valid chromosomes found in Bismark file
    chrmNumList = []
    
    for row in reader:
        # rowSplit = row.decode().split("\t")
        chrmNum = row.decode().split("\t")[0].split("chr")[1]
        if chrmNum in validChrm:
            validChrmNum = chrmNum
            # add each row of a valid chromosome in the Bismark file to chrms[key = chrmNum]
            chrms.setdefault(validChrmNum,[]).append(row)
            if validChrmNum not in chrmNumList:
                chrmNumList.append(validChrmNum)
    
    # multiprocessing using pool.map()  
    pool = mp.Pool(processes=cores)
    # Permeth 
    pool.map(writeHeader_permeth, chrmNumList)
    pool.map(writeData_permeth, chrmNumList)
    # DSS
    pool.map(writeHeader_dss, chrmNumList)
    pool.map(writeData_dss, chrmNumList)

endTime = time.time()
print("Time = " + str(endTime-startTime))

