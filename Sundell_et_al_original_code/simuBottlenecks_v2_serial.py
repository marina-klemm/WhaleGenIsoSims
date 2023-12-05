
#!usr/bin/python
#
# Filename: simuBottlenecks_v2.py
# Author : Argeopop (Martin Heger, Juhana Kammonen, Paivi Onkamo, Tarja Sundell)
# Purpose : simuPOP simulation of Finnish population histories ("Pullonkaulasimulaatiot 2.0"):
#
# February 2011 -> This script will use a configuration text file so that there is no need to modify
#                  an entire script for every different model used.
# 

import sys

import random
import math
from simuOpt import setOptions
setOptions(optimized=False, alleleType='short')
from simuPOP import *
from simuPOP.sampling import drawRandomSamples

# GENERATION OPTIONS (1 gen = 10 years)

staticPhaseEnd = 200 # 9000 BP - begin steady growth at latest
firstBalance = 400 # balancing growth before population peak
firstDecline	= 525 # 5750 BP - decline begins at population peak
firstMinimum = 640 # the decline evens out before bottleneck
secondBalance	= 690 # 4100 BP - first bottleneck begins - 
firstRecovery	= 850 # 3800 BP - first bottleneck ends
thirdBalance	= 900 # population growth evens out again 
secondDecline	= 950 # decline begins at population peak

firstCheck      = 200 # generation of first real "checkpoint" - population split @ 7000 BP
secondCheck     = 400 #
thirdCheck      = 401
fourthCheck    = 525 # population peak @ 5750 BP
fifthCheck     = 635
sixthCheck     = 691 # population minimum @ 4090 BP
seventhCheck     = 900 # fourth checkpoint - second population minimum
eighthCheck     = 960 #
ninthCheck      = 1099 # fifth checkpoint - present day, final population

#NOTE : There are actually six checkpoints, the very first being already in initOps[]

seq_file = "crs_based_sequences.txt" # the file containing the DNA-sequences with which the individuals are initialized
saami_seq_file = "saamiBgpop_sequences.txt"
# necessary filenames are saved as string objects and the respective files are opened:

# Arguments are read from configuration file
argv = open("conf_file.txt","r")
argList = []

while True:
 line = argv.readline()
 if line == "":
    break
 line = line.rstrip("\n")
 argList.append(line)

print argList
    
sampf = "bn_samples_%s_%s.txt" % (argList[1], sys.argv[1])
typef = "bn_mtHaplotype_frequencies_%s_%s.txt" % (argList[1], sys.argv[1])
#groupf = "bn_mtHaplogroup_frequencies_%s_%s.txt" % (argList[2], argList[1])
Y_typef = "bn_YSTR_frequencies_%s_%s.txt" % (argList[1], sys.argv[1])

YSTR_arlequinf = "bn_YSTR_arlequin_chunk_%s_%s.txt" % (argList[1], sys.argv[1])
mt_arlequinf = "bn_mt_arlequin_chunk_%s_%s.txt" % (argList[1], sys.argv[1])
seq_arlequinf = "bn_seq_arlequin_chunk_%s_%s.txt" % (argList[1], sys.argv[1])

run_statsf = "bn_run_stats_07%s.txt" % (argList[1])

YSTR_outputf = "bn_907_frequencies_%s.txt" % (sys.argv[1])
mt_outputf = "bn_832_frequencies_%s.txt" % (sys.argv[1])
mySampler_statsf = "bn_907_832_stats_07%s.txt"  % (argList[1])

# load the background populations
archaicEuropeanBgPop = loadPopulation("archaicEuropeanBackgroundPop_0.txt")
archaicScandinavianBgPop = loadPopulation("archaicScandinavianBackgroundPop_0.txt")
saamiBgPop = loadPopulation("saamiBackgroundPop_0.txt")

# the variables limiting ages of individuals:

maleMinMating = 2 # 20 years
maleMaxMating = 6 # 60 years

femaleMinMating = 2 # 20 years
femaleMaxMating = 4 # 40 years

maxAge = 6 # equals 60 yrs - maximum age of an individual

numGens = 1100 # number of generations - need to fix it, just do it here =)

# the global DNA-sequence list which associates a list index with a specific DNA-sequence
# NOTE: integer numbers are associated with bases as follows: [0 = A | 1 = T | 2 = C | 3 = G]
seqList = list()

# THIS IS THE GLOBAL LIST FOR CHECKPOINT GENERATION INDICES IN THE SIMULATION - NEED TO CHANGE THEM? JUST DO IT HERE:
genList = [0, 200, 400, 401, 525, 691, 900, 960, 1099]

# INTERNAL MIGRATION PROBABILITIES:
femnear = float(argList[2])
femfar = float(argList[3])
malnear = float(argList[4])
malfar = float(argList[5])

# BEGIN USER DEFINED FUNCTIONS (in CASE-INSENSITIVE alphabetical order):
def appendAverages(source, target): # appends averages from source list to target list

    global genList
    
    sp_counter = 0 # same purpose as generation counter
    gen_counter = 0
    while sp_counter < 3:
        while gen_counter < len(genList): # for all sampled generations
            if sp_counter == 2 and gen_counter == 0:
                gen_counter = gen_counter + 1
            else:
                sumList = list() # list for allele counts in generations indicated by gen_counter
                for i in range(0, len(source)):
                    #print ("%d | %d | %d | %d") % (sp_counter, gen_counter, source[i][0], source[i][1])
                    if source[i][0] == sp_counter and source[i][1] == gen_counter:

                        sumList.append(source[i][2])
                        #print "Append was made."
                    
                #NOTE: from now on the generation is indicated only by the average's index in the target list:
                if len(sumList) > 0:
                    average = (sum(sumList)/len(sumList))
                    target.append(average)

                gen_counter = gen_counter + 1

        sp_counter = sp_counter + 1
        gen_counter = 0 # generation counter is reset
        
    return True

def appendAveragesTotal(source, target): # special appending function for bottleneck simulations

    global genList
    sumList = list()
    marker = 0
    for i in range(0, len(genList)):
        for j in range(marker, marker+10): # for each generation:
            sumList.append(source[j])
        average = (sum(sumList)/len(sumList))
        target.append(average)
        marker = marker + 10 # increment 10 steps and start over, repeat len(genList) times

    return True
    
def appendCount(counters, target_allele, ploidy, chrom, ind):
	
    x = list()
    loci = range(len(target_allele)) # a table of loci
                    
    for locus in loci:
        x.append(ind.allele(locus, ploidy, chrom))

    if ploidy == 0 and chrom == 0: # checks whether it is sequence data we are handling

        found = False
        
        for i in range(0, len(seqList)):
            if found == False:
                if seqList[i] == x:
                    x = i # assigns sequence callsign to x
                    found = True
                
    check = exists(counters, x)
    if check >= 0:
            fix = counters[check]
            oldCount = getY(fix)
            newCount = oldCount + 1
            fix = setY(fix, newCount)
            counters[check] = fix
    else:
            counters.append((x, 1))

    return True

def cleanUp(pop): # cleans up the individuals:
                  # individuals in mating ages are marked for VSP split,
                  # individuals whose age exceeds maxAge are removed from population.
                  # individuals who have not experienced migration are saved to a list to be used later

    global maleMinMating, maleMaxMating, femaleMinMating, femaleMaxMating, maxAge
                  
    for i in range(0, pop.numSubPop()):

        for idx, individual in enumerate(pop.individuals(i)):

            myAge = individual.age
            
            if myAge > maxAge:
                individual.group = 3 # this individual will end up in a group that neither mates nor copies in this generation (=dies)
                
            else:
                ## BO: MALE/FEMALE is more readable ...
                if individual.sex() == MALE:
                    if myAge < maleMinMating:
                        individual.group = 0
                    if maleMinMating <= myAge <= maleMaxMating: # if in mating age:
                        individual.group = 1 # this guy is in mating age
                    if maleMaxMating < myAge <= maxAge:
                        individual.group = 2 # this guy has past matingAge

                if individual.sex() == FEMALE:
                    if myAge < femaleMinMating:
                        individual.group = 0
                    if femaleMinMating <= myAge <= femaleMaxMating: # if in mating age:
                        individual.group = 1 # this girl is in mating age
                    if femaleMaxMating < myAge <= maxAge:
                        individual.group = 2 # this girl has past matingAge
        
    return True

def computeMigrators(rate, popSize):

    if rate == 0.0: # this means rate has been defined zero in the conf file -> model with no external migration
        return 0

    elif int(rate * popSize) == 0: # rate makes number of migrators less than 1

        if random.random() < rate: # at a very small probability...
            return 1 # ...add 1 migrator
        else:
            return 0

    else:
        return int(rate * popSize)
        
def demo(gen, pop):  # demographics function to control population growth (actually birth rate now).
                     # number of bottlenecks and the severities are as agreed on 01.09.2009. Check the generation "names" from
                     # the beginning of the script

    if gen <= staticPhaseEnd:
        return [eval(argList[9]), eval(argList[10])]
    
    elif gen <= firstBalance:
        return [eval(argList[11]), eval(argList[12])] # (SAA, MUU) growth until max
    
    elif gen < firstDecline:
        return [eval(argList[13]), eval(argList[14]), eval(argList[15])] #(SAA, SW, NE)
        # popsize evens out, drops slightly
        # bottleneck events around this time

    elif gen < firstMinimum:
        return [eval(argList[16]), eval(argList[17]), eval(argList[18])]
    
    elif gen < secondBalance:
        return [eval(argList[19]), eval(argList[20]), eval(argList[21])]

    elif gen < firstRecovery:
        return [eval(argList[22]), eval(argList[23]), eval(argList[24])]

    elif gen < thirdBalance:
        return [eval(argList[25]), eval(argList[26]), eval(argList[27])]
    
    elif gen < secondDecline:
        return [eval(argList[28]), eval(argList[29]), eval(argList[30])]
    
    else: #final rise - growth from bottleneck until present
        return [eval(argList[31]), eval(argList[32]), eval(argList[33])]

def dumpSequences(): # dumps sequences from seqList in the end of sample file
                     # has to be implemented as a separate function so that the dump is done only when sampling for the run has finished

    global seqList

    output_file = open(sampf, "a")

    linelist = list()

    for i in range(len(seqList)):
        linelist.append("%d\t" % i)
        for j in range(len(seqList[i])):
            if seqList[i][j] == 0:
                linelist.append("A")
            if seqList[i][j] == 1:
                linelist.append("C")
            if seqList[i][j] == 2:
                linelist.append("G")
            if seqList[i][j] == 3:
                linelist.append("T")
            if seqList[i][j] == 4:
                linelist.append("-")
        linelist.append("\n")
            
    output_file.writelines(linelist)
    output_file.close
    
    return True

def exists(counters, a): # checks whether object a exists in list 'counters' and
                         # returns the index of first encountered object a in the list
    length = len(counters)
    if length == 0: # if length of the list is zero a obviously doesn't exist...
        return -1   # ...and so we return -1
    else:
        for j in range(0, length): # search through the whole list
            if getX(counters.__getitem__(j)) == a: # if object a is encountered...
                return j # return the index of this FIRST object a in the list
        return -1 # if object a is not encountered during the loop, return -1

def frequencies(samples, gen, targetString, target_filename, ploidy, chrom, switch): # loci is a range of target loci,
                                                              # switch indicates whether its sequence-data we're looking at

    loci = eval(targetString) # evaluates parameter string

    output_file = open(target_filename, "a")

    linelist = list()
    linelist.append("SAMPLES @ GEN %d\n" % (gen))
    output_file.writelines(linelist)

    for i in range(0,len(samples)): # for all samples
        smpPop = samples[i] # take i-th pop from the list
        smpPop_total = samples[i].clone() # clone the i-th pop from the list for later use
        subs = smpPop.subPopSizes() # subs will become a list of subpopulation sizes in smpPop
        linelist = list() # list of strings to be appended to the text file specified above

        linelist.append('STATS FOR SAMPLE ' + str(i+1) + ':\n')

        total_counters = list() # creates a list for total allele counters
        total_size = smpPop_total.popSize()

        neiList_all = list() # this is a list for reserving the computed Nei's gene diversities for subpops
        
        for j in range(0,smpPop.numSubPop()): # for every subPop
            
            counters = list() # creates a list for allele counters
            for idx,ind in enumerate(smpPop.individuals(j)): # for every individual in this subPop

                appendCount(counters, loci, ploidy, chrom, ind)    

            linelist.append('Subpopulation\tAllele\tIndividuals\tFreq. in sample\n')
            
            #Let's print the results then:
            allTot = len(counters) # number of different alleles in sample, equals to the length of 'counters'
            singTot = singletons(counters)
            freqSum = 0 # sum of total frequencies. Should be 1 in the end.
            neiList = list() # list for computing Nei's gene diversity

            for idx,fromTuple in enumerate(counters): # for every tuple in counters
                carriers = float(getY(fromTuple)) # get the number of allele carriers and float it
                sampleSize = float(subs[j]) # get size of j-th (=present) sample subpopulation and float it
                freq = carriers/sampleSize # the frequency for allele found counters

                neiList.append(pow(freq, 2)) # appends the alleleFreq freq squared into neiList
                
                freqSum = freqSum + freq # update the sum of frequencies
                linelist.append('\t' + str(j) + '\t' + str(getX(fromTuple)) + '\t' + str(int(carriers)) + '\t\t' + str(freq) + '\n')
                
            # Compute Nei's gene diversity for this subpop:
            nei = neisGD(sampleSize, neiList) 
            neiList_all.append(nei) # append computed GD into a specific list for later use

            linelist.append('Different haplotypes in subpop: ' + str(allTot) + '\n')
            linelist.append('Singleton haplotypes in subpop: ' + str(singTot) + '\n')
            linelist.append('Nei\'s gene diversity for haplotypes in subpop: ' + str(nei) + '\n')
            linelist.append('Sum of frequencies in subpop: ' + str(freqSum) + '\n')
            linelist.append('|----------------------------------------------|' + '\n')

        output_file.writelines(linelist)

        linelist = list() # this author saw that it was better to create a new list for these lines 
        smpPop_total.mergeSubPops() # merges the cloned pop's subpops into single pop with 1500 individuals
        for idx,ind in enumerate(smpPop_total.individuals(0)): # there's only one subpop (0) so this suffices

            appendCount(total_counters, loci, ploidy, chrom, ind)
    
        # Let's print these results too:
        allTot_total = len(total_counters)
        singTot_total = singletons(total_counters)
        freqSum = 0 # sum of total frequencies. Should be 1 in the end.
        neiList = list() # this resets the neiList in the above for loop
        linelist.append('TOTAL STATS FOR SAMPLE %d:\nFreq. in sample\t\tIndividuals\tHaplotype\n' % (i+1))
    
        for idx,fromTuple in enumerate(total_counters): # for every tuple in total_counters
            carriers = float(getY(fromTuple))
            freq = carriers/(float(total_size))

            neiList.append(pow(freq, 2)) #append the freq for found allele squared in neiList
        
            freqSum = freqSum + freq
            linelist.append(str(freq) + '\t\t' + str(int(carriers)) + '\t' + str(getX(fromTuple)) + '\n')

        # Compute Nei's gene diversity and FST for total sample:
        nei_total = neisGD(total_size, neiList)
        #nei_average = (sum(neiList_all)/len(neiList_all))

        #fst_total = (nei_total - nei_average)/(nei_total)    
    
        linelist.append('Different haplotypes in sample: ' + str(allTot_total) + '\n')
        linelist.append('Singleton haplotypes in sample: ' + str(singTot_total) + '\n')
        linelist.append('Nei\'s gene diversity for haplotypes in sample: ' + str(nei_total) + '\n')
        #linelist.append('Total Fst value in sample: ' + str(fst_total) + '\n')
        linelist.append('Sum of frequencies in sample: ' + str(freqSum) + '\n')
        linelist.append('|----------------------------------------------|' + '\n')

        output_file.writelines(linelist)

    output_file.close()    
    return True

# function that dumps mitochondrial data as one chunk in arlequin format
# NOTE: data produced by this function does not fit Arlequin until fileChopper_bn.py has been executed on arl_file
def freqArlequinDump(samples, arl_filename, gen):

    output_file = open(arl_filename, "a")
    
    # let us compute the mitochondrial haplotype frequencies from samples
    sampFreqs = list() # list of lists
                   
    for i in range(0, len(samples)):

        newSampFreqList = list()
        sampFreqs.append(newSampFreqList) # to be populated with tuples of format (allele, (sample, freqInSample))
                                          # NOTE: freqInSample is absolute, not relative
                              
        for idx, ind in enumerate(samples[i].individuals(0)):
    
            haplotype = ind.allele(0)
            check = exists(sampFreqs[i], haplotype) # check if haplotype has already been seen in samples.
            if check >= 0: # allele exists in list, increment the old counter:
                fix = sampFreqs[i].__getitem__(check)
                oldCount = getY(getY(fix))
                newCount = oldCount + 1
                fix = (haplotype, (i, newCount)) # replaces list element fix with updated data
                sampFreqs[i].__setitem__(check, fix)
            else: # allele does not exist in the list, add new counter for the allele to the list:
                sampFreqs[i].append((haplotype, (i, 1)))
    
    linelist = list()

    for i in range(0, len(samples)):
        
        linelist.append('[Profile]\n')
        linelist.append('\tTitle=\"Frequency data (gen %d) - Sample %d\"\n' % (gen,i))
        linelist.append('\tNbSamples=%d\n' % 1)
        linelist.append('\tGenotypicData=%d\n' % 0)
        linelist.append('\tLocusSeparator=TAB\n')
        linelist.append('\tDataType=FREQUENCY\n')
        linelist.append('\tFrequency=REL\n')

        linelist.append('[Data]\n')
        linelist.append(' [[Samples]]\n')

        sampleSize = samples[i].subPopSize()
            
        linelist.append('\tSampleName=\"Subpopulation %d\"\n' % 0)
        linelist.append('\tSampleSize=%d\n' % sampleSize)
        linelist.append('\tSampleData={\n')

        # data is now retrieved from sampFreqs:

        for element in sampFreqs[i]:
            haplotype = getX(element)
            sampleNumber = getX(getY(element)) # gets the first object of element  
            if sampleNumber == i: # if first object from (sample, freqInSample) matches current sample
                absFreq = getY(getY(element))
                linelist.append('%d\t' % haplotype + str(float(absFreq)/sampleSize) + "\n")

        linelist.append('}\n')
                                
    output_file.writelines(linelist)
    output_file.close()
    
    return True

def getX((x,y)):
    return x # returns first object of pair (x,y)

def getY((x, y)):
    return y # returns second object of pair (x,y)
    

def getYSTRs(string): # special "makeString" -function to get all 16 YSTR-loci
    
    #IMPORTANT NOTE: this script will not recognize numbers which contain "e" in their notation
    
    chars = list() # list for saving numeral characters of a string
    YSTRs = list() # list of YSTR-alleles contained by the string in order by locus
    tabCounter = 0 # counts the number of tabs encountered
    YSTRcounter = 0
    begins = False # flag that tells when the numbers of interest begin

    for char in string:

        if begins == False:
            if char == '\t':
                tabCounter = tabCounter + 1
                
            if tabCounter == 3:
                begins = True
                
        else: # heads up! the number of interest begins
            
            if char.isdigit(): #if character read is a digit or a separating them
                chars.append(char)

            elif char == ' ' : # this would mean that present YSTR-number has ended
                YSTR = int(makeString(chars)) # YSTR will become an int value of YSTRs in this locus
                YSTRs.append(YSTR) # ...and this value is appended into the list
                chars = list() # finally 'chars'-list is reset
                YSTRcounter = YSTRcounter + 1
                if YSTRcounter == 16: # which would indicate all loci for this individual are read
                    return YSTRs # returns a list like (14,13,25,16,17,....)

    return True

def initGroups(pop): # initializes infoField 'group' for all individuals in pop

    for i in range(0, pop.numSubPop()):
        for idx, ind in enumerate(pop.individuals(i)):
            if (ind.sex() == MALE) and (maleMinMating <= ind.age <= maleMaxMating):
                ind.group = 1

            elif (ind.sex() == MALE) and (maleMaxMating <= ind.age <= maxAge): # this individual has passed mating age
                ind.group = 2
                
            elif (ind.sex() == FEMALE) and (femaleMinMating <= ind.age <= femaleMaxMating):
                ind.group = 1

            elif (ind.sex() == FEMALE) and (femaleMaxMating <= ind.age <= maxAge): # this individual has passed mating age
                ind.group = 2

            else: # this individual is too young to mate:
                ind.group = 0 # group 0 will contain only individuals too young to mate - thus they cannot be anyones parents

    return True

def initPop(pop):

    global seq_file # imports the global-intended filename here
    
    initSex(pop)

    # initialize info field 'age' and divide pop into virtual subPops by gender and age:
    pop.setIndInfo([random.randint(0, maxAge) for x in range(pop.popSize())], 'age')


    # get Saami subpop sequences from seq_file and initialize mitochondria with them:
    sequences = initSeqValues(pop, saami_seq_file)
    initGenotype(pop,
               haplotypes=sequences,
               # proportions relative to the frequency of sequences chosen for initialization (33 pcs, 10 different sequences):  
               prop=[(1.0/7)]*7,               
               loci=range(631), ploidy=[0], subPops=[0])


    # get Muu Suomi subpop sequences from seq_file and initialize mitochondria with them:
    sequences = initSeqValues(pop, seq_file)
    
    initGenotype(pop,
               haplotypes=sequences,
               ## proportions relative to the frequency of sequences chosen for initialization (33 pcs, 10 different sequences):  
               prop=[(1.0/33),(5.0/33),(7.0/33),(1.0/33),(3.0/33),(3.0/33),(1.0/33),(2.0/33),(7.0/33),(3.0/33)],               
               loci=range(631), ploidy=[0], subPops=[1])

    # initialize Y-chromosome:
    initGenotype(pop,
                haplotypes=([14,14,24,16,17,14,15,14,11,10,21,14,12,14,10,19],  # UU30   x 18
                        [14,12,23,16,16,14,31,13,10,10,22,11,11,16,10,20], # UU110  x 8
                        [14,12,23,16,15,14,31,13,10,10,22,11,11,16,10,20], # UU27   x 6
                        [14,12,23,16,16,14,31,13,10,10,23,11,11,16,10,20], # VA2    x 5
                        [14,13,23,16,17,14,15,14,11,10,21,14,12,14,10,19], # UU31   x 5
                        [14,14,24,16,17,14,15,14,10,10,21,14,12,14,10,19], # UU8    x 5
                        [14,12,23,16,15,14,31,13,10,10,21,11,11,16,10,20], # UU61   x 4
                        [14,13,24,16,17,14,15,15,11,10,21,14,12,14,10,19], # UU51   x 4
                        [14,14,24,16,18,14,15,14,11,10,21,14,12,14,10,19], # UU134  x 4
                       ),                
                prop=[(18.0/59), (8.0/59), (6.0/59), (5.0/59), (5.0/59), (5.0/59), (4.0/59), (4.0/59), (4.0/59)],
                
                loci=range(631, 647), ploidy=[1])
           
    return True

def initSeqValues(pop, seq_file): # the goal here is to open a sequence text file and read sequences

    sequences = list() # list of character lists which this function populates and eventually returns
    
    source_file = open(seq_file, "r")
    readList = source_file.readlines()
    source_file.close() # done reading, close source file

    for string in readList:

        baseList = list() # create new list
        for char in string: # for every character in this string
            if char != "\n": # if not end of line
                base = eval(char) # evaluated char should be a base number (4 if deleted allele)
                baseList.append(base) # append base number into 'charList'
            
        sequences.append(baseList) # finally append the list per se into 'sequences'        

    # 'sequences' is now a ready list of character lists, let's return it:
    return sequences

def insertMigration(pop, param):

    gen = pop.dvars().gen

    # CONTINUOUS SMALL GENE FLOW ('drizzle'):
    # ---------------------------------------
    
    if gen <= 400: # 2-subpop-phase

        if param[0] == 0: # In this case function call concerns Saami subpop:

            migrPop = pop.extractSubPops(subPops=[param])            
            migrators = computeMigrators(float(argList[36]), pop.popSize())
            print 'Migrators to Saami: %d' % migrators
            backgroundPop = saamiBgPop       
            insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])

            migrators = computeMigrators(float(argList[38]), pop.popSize())       
            backgroundPop = archaicEuropeanBgPop       
            insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])


        else: # param[0] == 1: # Muu Suomi (SW from 400 onward)

            migrPop = pop.extractSubPops(subPops=[param])
            migrators = computeMigrators(float(argList[38]), pop.popSize())
            backgroundPop = archaicEuropeanBgPop       
            insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])

    else: # gen > 400 (3-subpop-phase)

        if param[0] == 0: # In this case function call concerns Saami subpop:

            migrPop = pop.extractSubPops(subPops=[param])           
            migrators = computeMigrators(float(argList[36]), pop.popSize())
            print 'Migrators to Saami: %d' % migrators
            backgroundPop = saamiBgPop       
            insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])

            migrators = computeMigrators(float(argList[38]), pop.popSize())       
            backgroundPop = archaicEuropeanBgPop       
            insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])

        elif param[0] == 1: # SW from 400 onward

            migrPop = pop.extractSubPops(subPops=[param])            
            migrators = computeMigrators(float(argList[38]), pop.popSize())     
            backgroundPop = archaicEuropeanBgPop       
            insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])

            if gen > 750: # add Archaic Scandinavian migration

                migrPop = pop.extractSubPops(subPops=[param])            
                migrators = computeMigrators(float(argList[37]), pop.popSize())     
                backgroundPop = archaicScandinavianBgPop       
                insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])
                

        else: # param[0] == 2: # NE from 400 onward

            migrPop = pop.extractSubPops(subPops=[param])            
            migrators = computeMigrators(float(argList[38]), pop.popSize())     
            backgroundPop = archaicEuropeanBgPop       
            insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])

            migrators = computeMigrators(float(argList[36]), pop.popSize())       
            backgroundPop = saamiBgPop       
            insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])

    # MIGRATION WAVES:
    # ----------------

    ## Typical Comb Ware migration wave (CA 6000 BP) from Archaic European background population:

    if gen in [499, 500, 501] and param[0] == 2:

        migrPop = pop.extractSubPops(subPops=[param[0]]) # NE-subpop
        migrators = computeMigrators(float(argList[34]), pop.popSize())
        backgroundPop = archaicEuropeanBgPop
        insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])

    ## Corded Ware migration wave (CA 5200 BP) from Archaic European background population:

    if gen in [579, 580, 581] and param[0] == 1:

        migrPop = pop.extractSubPops(subPops=[param[0]]) # SW-subpop
        migrators = computeMigrators(float(argList[35]), pop.popSize())
        backgroundPop = archaicEuropeanBgPop
        insertMigrators(pop, migrPop, migrators, backgroundPop, param[0])
        
    return True

def insertMigrators(pop, migrPop, migrators, backgroundPop, subPop):

    gen = pop.dvars().gen

    if migrators == 0:

        return True # no migration happens

    else: # migrators was greater than zero

        print "Destination subpop: %d" % subPop
        print "Migrators: %d" % migrators
        initialSplitter = ProductSplitter(splitters=[SexSplitter(), InfoSplitter('group', values=[0, 1, 2, 3])])
        interSplitter = CombinedSplitter(splitters=[initialSplitter], vspMap=[(0, 4), (1, 5), (2, 6), (3, 7)])

        finalSplitter = CombinedSplitter(splitters=[interSplitter, SexSplitter()])

        backgroundPop.setVirtualSplitter(finalSplitter)

        samples = drawRandomSamples(backgroundPop, sizes=migrators, numOfSamples=1) # about the same as in earlier simulations
        migrPop.addIndFrom(samples[0])
        migrPop.mergeSubPops()

        pop.removeSubPops(subPops=[subPop])

        if gen <= 400: # 2-subpop phase

            if subPop == 0:
                sparePop = pop.extractSubPops(subPops=[0])
                pop.addIndFrom(migrPop)
                pop.removeSubPops(subPops=[0])
                pop.addIndFrom(sparePop)

            else: # subPop == 1:
                # subPops in proper order, nothing to do but add migrPop
                pop.addIndFrom(migrPop)

        else: # gen > 400 # 3-subpop phase

            if subPop == 0:
                firstSparePop = pop.extractSubPops(subPops=[0]) # spare SW
                pop.removeSubPops(subPops=[0]) # remove SW
                pop.addIndFrom(migrPop) # add Saami
                secondSparePop = pop.extractSubPops(subPops=[0]) # spare NE
                pop.removeSubPops(subPops=[0]) # remove NE
                pop.addIndFrom(firstSparePop) # add SW
                pop.addIndFrom(secondSparePop) # add NE

            elif subPop == 1:
                sparePop = pop.extractSubPops(subPops=[1]) # spare NE
                pop.removeSubPops(subPops=[1]) # remove NE
                pop.addIndFrom(migrPop) # add SW
                pop.addIndFrom(sparePop) # add NE

            else: # subPop == 2
                # subPops in proper order, nothing to do but add migrPop
                pop.addIndFrom(migrPop)

        return True

def lethalEvent(pop): # the idea here is to affect a portion of the population with
                      # a some kind of an lethal effect that "kills" individuals by simply removing them.
                      # This is a step closer to the "more realistic" simulation model.

    gen = pop.dvars().gen

    #Could be something like this:
    indices = range(pop.popSize()) # a list of individual indices
    numAffected = int(.15 * pop.popSize()) # this would make the affection ratio 15% of current population

    affected = random.sample(indices, numAffected) # picks numAffected indices

    #finally remove selected individuals
    pop.removeIndividuals(affected)

    return True

def loadBgPops(pop):

    gen = pop.dvars().gen

    if gen == 100:
        
        archaicEuropeanBgPop = loadPopulation("archaicEuropeanBackgroundPop_200.txt")
        archaicScandinavianBgPop = loadPopulation("archaicScandinavianBackgroundPop_200.txt")
        saamiBgPop = loadPopulation("saamiBackgroundPop_200.txt")
        
    if gen == 300:
        
        archaicEuropeanBgPop = loadPopulation("archaicEuropeanBackgroundPop_400.txt")
        archaicScandinavianBgPop = loadPopulation("archaicScandinavianBackgroundPop_400.txt")
        saamiBgPop = loadPopulation("saamiBackgroundPop_400.txt")
        

    if gen == 500:
        
        archaicEuropeanBgPop = loadPopulation("archaicEuropeanBackgroundPop_600.txt")
        archaicScandinavianBgPop = loadPopulation("archaicScandinavianBackgroundPop_600.txt")
        saamiBgPop = loadPopulation("saamiBackgroundPop_600.txt")
        

    if gen == 700:
        
        archaicEuropeanBgPop = loadPopulation("archaicEuropeanBackgroundPop_800.txt")
        archaicScandinavianBgPop = loadPopulation("archaicScandinavianBackgroundPop_800.txt")
        saamiBgPop = loadPopulation("saamiBackgroundPop_800.txt")

    if gen == 900:

        archaicEuropeanBgPop = loadPopulation("archaicEuropeanBackgroundPop_1000.txt")
        archaicScandinavianBgPop = loadPopulation("archaicScandinavianBackgroundPop_1000.txt")
        saamiBgPop = loadPopulation("saamiBackgroundPop_1000.txt")

    return True

def makeString(chars): #concatenates characters of list 'char' into one string
    
    numberString = ''
    for char in chars:
        numberString = numberString + char
        
    return numberString

def myDump(samples, target_filename, gen): # our own dumper =)
                                           # parameters: list of samples, file handle to target file, current generation number

    output_file = open(target_filename, "a")

    linelist = list() # list for "writelines"

    linelist.append("SAMPLES @ GEN %d\n" % (gen))
    
    for num,sample in enumerate(samples): # for every sample
        linelist.append("RAW DATA FROM SAMPLE " + str(num) + "\n")
        for sp in range(0, sample.numSubPop()): # for all subpops in sample
            linelist.append("Subpopulation " + str(sp) + ":\n" + "ind\tsex\thaptype\thapgroup\tY_haptype\n")

            for idx,ind in enumerate(sample.individuals(sp)): # for every individual in subpop
                
                gender = ind.sex() # if you want gender as character string please use sexChar() instead.
                
                ## access the DNA-sequence. map the sequence as a new haplotype:

                DNAsequence = list(ind.genotype(ploidy=0, chroms=0))

                haplotype = sequenceMap(DNAsequence) # calls mapper function, establishing this sequences haplotype number.
                                                     # sequence is saved in global list 'seqList'

                # insert a line containing data into linelist.
                linelist.append(str(idx) + "\t" + str(gender) + "\t" + str(haplotype) + "\t")

                for locus in range(0,16): # for loci with indices 0-15 in Y_chrom
                    linelist.append(str(ind.allele(locus, ploidy=1, chrom=2)) + " ")

                linelist.append("\n") # line feed, carriage return

    output_file.writelines(linelist)

    output_file.close()
    
    return True

def mySampler(source_file, target_file, sampleSize, clue): # takes a specific sample from a myDump-generated population dump 
                                              # clue is to indicate whether we're counting mt or YSTR data                                              
    # BEGIN PART ONE:
    #================
    
    input_file = open(source_file, "r")
    readList = input_file.readlines() # readlist for samples
    input_file.close()

    sp_number = -1 #indicates current subpopulation number
    sample_counter = -1 # counts the samples to indicate
    interval = True # boolean value that indicates that we're currently in a header interval of the file
    irr_sample = True # boolean value indicating an irrelevant sample.

    sampList = list() # list of lists of individuals (index indicating sample number)
    for i in range(0, 10):
        indList = list() # new list for individuals in sample i
        sampList.append(indList)

    # Lines of sample file are walked through.
    for string in readList:

        if "SAMPLES @ GEN 1099" in string:
            irr_sample = False # and will remain False for the rest of this phase.           

        if irr_sample == False:

            if interval: 

                if "RAW DATA FROM SAMPLE" in string:
                    interval = True
                    chars = readNumbers(string, 0)
                    numberString = makeString(chars)
                    sample_number = int(numberString)

                if "Subpopulation" in string:
                    chars = readNumbers(string, 0)
                    numberString = makeString(chars)
                    sp_number = int(numberString)

                if "ind\tsex" in string:
                    interval = False # heads up! Interval section has ended and following string contains information of interest

            else:

                    if "RAW DATA FROM SAMPLE" in string: # new interval section begins.
                        interval = True
                        chars = readNumbers(string, 0)
                        numberString = makeString(chars)
                        sample_number = int(numberString)

                    elif "Subpopulation" in string:
                        chars = readNumbers(string, 0)
                        numberString = makeString(chars)
                        sp_number = int(numberString)

                    elif "ind\tsex" in string:
                        interval = False # heads up! Interval section has ended and following string contains information of interest

                    else: # If control reaches this point the target string is certain to contain individual numeric data:

                        if clue == "mt":
                            
                            chars = readNumbers(string, 2) # read haplotype from the list
                            
                            if (chars != None) and (len(string) < 500): # while end of data is not reached
                                
                                numberString = makeString(chars)
                                haplotype = int(numberString)
                                target_list = sampList[sample_number] # sets the individual list corresponding current sample as target
                                target_list.append((haplotype, sp_number))

                        else: # clue == "YSTR"

                            YSTR_list = getYSTRs(string) # read YSTRs from this string

                            if (len(string) < 500): # while end of data is not reached

                                target_list = sampList[sample_number]
                                target_list.append((YSTR_list, sp_number))
                        
    # sample the lists ('sampleSize' samples):
    randomSamples = list()

    for i in range(0,len(sampList)):

        newRandomList = random.sample(sampList[i], sampleSize) # pick sampleSize individuals randomly from sampList[i]
        randomSamples.append(newRandomList) # append the random samples into a list of random samples.

    output_file = open(target_file, "a")
        
    for i in range(0,len(randomSamples)): # for every sample

        linelist = list() # new list for output lines

        total_size = len(randomSamples[i])
        total_counters = list()
        neiList = list() # helping list for computing Nei's gene diversity

        # following loop walks through all individuals of the sample and counts the number of every haplotype    
        for pair in randomSamples[i]:
                x = pair[0] # takes the first element of the pair (haplotype or YSTRs dependig of switch)
                check = exists(total_counters, x) # check if element x has already been seen in this sample.
                if check >= 0: # element exists in list, increment the old counter:
                    fix = total_counters[check]
                    oldCount = fix[1]
                    newCount = oldCount + 1
                    fix = (x, newCount)
                    total_counters[check] = fix
                else: # element does not exist in the list, add new counter for the element to the list:
                    total_counters.append((x, 1))
        
        # Let's print the results:
        allTot_total = len(total_counters)
        singTot_total = singletons(total_counters)
        freqSum = 0 # sum of total frequencies. Should be 1 in the end.

        linelist.append('STATS FOR mySampler SAMPLE %d:\nAllele\tIndividuals\tFreq. in sample\n' % i)
        
        for idx,fromPair in enumerate(total_counters): # for every pair in total_counters
            carriers = float(fromPair[1])
            freq = carriers/(float(total_size))

            neiList.append(pow(freq, 2)) #append the freq for found allele squared in neiList
            
            freqSum = freqSum + freq
            linelist.append(str(fromPair[0]) + '\t' + str(int(carriers)) + '\t\t' + str(freq) + '\n')

        # Compute Nei's gene diversity and FST for total sample:
        nei_total = neisGD(total_size, neiList)
        #nei_average = (sum(neiList_all)/len(neiList_all))
            
        linelist.append('Different haplotypes in sample: ' + str(allTot_total) + '\n')
        linelist.append('Singleton haplotypes in sample: ' + str(singTot_total) + '\n')
        linelist.append('Nei\'s gene diversity for haplotypes in sample: ' + str(nei_total) + '\n')
        linelist.append('Sum of haplotype frequencies in sample: ' + str(freqSum) + '\n')
        linelist.append('|----------------------------------------------|' + '\n')

        output_file.writelines(linelist)

    output_file.close()

def myStatsCount(mt_file, Y_file, mt_sampleSize, Y_sampleSize, stats_file): # handles data reading and file output for mySampler samples.
                                                                            # Goal is to get data output similar to that of JU Palo et al. 2009.

    # mySampler has created a 832 ind sample data file 
    # We now go ahead and read this file:
    
    input_file = open(mt_file, "r")

    #read all lines of this file into one list:
    readlist = input_file.readlines()

    mt_haplotypes_l = list() # list for saving number of alleles from every subpop. format (sp,count) where sp=subPop & count=count 
    mt_singletons_l = list()
    mt_diversities_l = list()

    sample_counter = -1 # the sample counter (there's 10 samples)

    # following for loop walks through the lines of the input file. Crucial data is stored into lists to be used later
    for string in readlist:

            if "Different haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                mt_haplotypes_l.append((count))
                
            if "Singleton haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                mt_singletons_l.append((count))

            if "Nei\'s gene diversity for haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                mt_diversities_l.append((count))

    input_file.close()
    
    #average for number of mt haplotypes in this run (average of 10 samples as specified by Martin on 2.10.2009)
    mt_hapTot_average = (sum(mt_haplotypes_l)/len(mt_haplotypes_l))

    #average for mt gene diversities in this run
    mt_divTot_average = (sum(mt_diversities_l)/len(mt_diversities_l))

    #average for mt singletons in this run (NOTE: singletons not present in Palo data)
    mt_singTot_average = (sum(mt_singletons_l)/len(mt_singletons_l))

    # mySampler also creates a 907 ind sample data file for YSTR-data
    # we open this now:

    input_file = open(Y_file, "r")

    #read all lines of this file into one list:
    readlist = input_file.readlines()

    Y_haplotypes_l = list() # list for saving number of alleles from every subpop. format (sp,count) where sp=subPop & count=count 
    Y_singletons_l = list()
    Y_diversities_l = list()

    sample_counter = -1 # reset back to "NULL" value

    # following for loop walks through the lines of the input file. Crucial data is stored into lists to be used later
    for string in readlist:

            if "Different haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                Y_haplotypes_l.append((count))
                
            if "Singleton haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                Y_singletons_l.append((count))

            if "Nei\'s gene diversity for haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                Y_diversities_l.append((count))

    input_file.close()

    # this block computes the numbers we need and puts them into lists conveniently
    # so later all we need to do is just access the lists and print them

    #average for number of haplotypes in this run (average of 10 samples as specified by Martin on 2.10.2009)
    Y_hapTot_average = (sum(Y_haplotypes_l)/len(Y_haplotypes_l))

    #average for diversities in this run
    Y_divTot_average = (sum(Y_diversities_l)/len(Y_diversities_l))

    #average for singletons in this run (NOTE: singletons not present in Palo data)
    Y_singTot_average = (sum(Y_singletons_l)/len(Y_singletons_l))
        
    #these two stupid-looking lines create output file for reading:
    output_file = open(stats_file, "a")
    output_file.close()

    #these lines open the file again, this time for reading:
    output_file = open(stats_file, "r")

    checklist = output_file.readlines()
    output_file.close()

    if len(checklist) == 0: # when true the output file was empty and we need to print header row into first line of file:
        writelist = list() #new list for writable data
        output_file = open(stats_file, "a")
        writelist.append("Run")

        writelist.append("\tY_N") # sample size is introduced        
        writelist.append("\tY_A")
        writelist.append("\tY_GD")
        writelist.append("\tY_S\t") # singletons are count, not thrown to waste despite they are missing from Palo et al. data

        # same is applied for mitochondria...

        writelist.append("\tmt_N") # sample size is introduced        
        writelist.append("\tmt_A")
        writelist.append("\tmt_GD")
        writelist.append("\tmt_S\n") # singletons are count, not thrown to waste despite they are missing from Palo et al. data
        
        output_file.writelines(writelist) # writes header into output file
        output_file.close()

    writelist = list() # writelist is reset
    output_file = open(stats_file, "a")

    writelist.append("%s\t" % (sys.argv[1])) # this appends the run number into the beginning of a data row

    writelist.append(str(Y_sampleSize) + "\t")

    writelist.append(str(Y_hapTot_average) + "\t")

    writelist.append("%f" %(Y_divTot_average) + "\t")

    writelist.append(str(Y_singTot_average) + "\t\t")

    writelist.append(str(mt_sampleSize) + "\t")

    writelist.append(str(mt_hapTot_average) + "\t")

    writelist.append("%f" % (mt_divTot_average) + "\t")

    writelist.append(str(mt_singTot_average) + "\n")

    output_file.writelines(writelist)

    output_file.close()
    return True

def neisGD(n, neiList): # computes Nei's gene diversity

    nei = n*(1 - sum(neiList))/(n-1)
    return nei

def overspillLethalEvent(pop, numRemovables):

    if numRemovables <= 0:

        return True

    else:
        numRemovables = numRemovables + 50
        print 'Removables: %d' % numRemovables
        indices = range(pop.popSize()) # a list of individual indices
        print 'Indices: %d' % len(indices)
        affected = random.sample(indices, int(numRemovables))
        pop.removeIndividuals(affected)
        
    return True
    

def printSize(pop): # prints the current year (BP = before present) and current popSize into the console

    gen = pop.dvars().gen

    if sys.argv[1] == "1":
        popSizeFile = open("popSizes.txt", "a")
    writeList = list()

    line = 'Year ca. : %d AD (gen %d) | Population size is: ' % (-9000 + (pop.dvars().gen+1)*10, pop.dvars().gen) + str(pop.popSize())
    print line
    if sys.argv[1] == "1":
        writeList.append(line + "\n")

    if gen <= 400: # 2-subpop phase:
        
        line = 'Saami: %d | Muu Suomi: %d' % (pop.subPopSize(subPop=[0]), pop.subPopSize(subPop=[1]))
        print line
        if sys.argv[1] == "1":
            writeList.append(line + "\n")

    else: # gen > 400 (3-subpop phase):

        line = 'Saami: %d | SW: %d | NE: %d' % (pop.subPopSize(subPop=[0]), pop.subPopSize(subPop=[1]), pop.subPopSize(subPop=[2]))
        print line
        if sys.argv[1] == "1":
            writeList.append(line + "\n")

    if sys.argv[1] == "1":    
        popSizeFile.writelines(writeList)
        popSizeFile.close()
    
    return True # I now see why this has to be done - if only 'print' was called, program execution would halt altogether

def printNumSubPops(pop):

    print 'Num. of subPops is: % d' % pop.numSubPop()

    return True

def printVsp(pop): # prints the distribution of virtual subpops

    global finalSplitter
    
    for i in range(6):
        print finalSplitter.name(i)

    return True

def readNumbers(string, tabs): #parameters are a string of characters and the number of tabs to reach the number to
                               #be read from the string (2 for haplotype).

    #IMPORTANT NOTE: this script will not recognize numbers which contain "e" in their notation
    
    chars = list() # list for saving numeral characters of a string
    tabCounter = 0 # counts the number of tabs encountered
    begins = False # flag that tells when the numbers of interest begin
    
    if tabs == 0: # just check the number in the string
        chars = list() # list for saving numeral characters of a string
        for idx,char in enumerate(string):
            if char.isdigit() or char == '.': #if character read is a digit or a '.' separating them
                chars.append(char)

        return chars

    else:
        for char in string:

            if begins == False:
                if char == '\t':
                    tabCounter = tabCounter + 1
                
                if tabCounter == tabs:
                    begins = True
                
            else: # heads up! the number of interest begins
            
                if char.isdigit() or char == '.': #if character read is a digit or a '.' separating them
                    chars.append(char)

                else: # this would mean that the numbers of interest have ended as char is not a digit or "."
                    return chars # list 'chars' is ready so we return it

def removeOverspill(pop):

    gen = pop.dvars().gen
    
    if gen <= staticPhaseEnd:
        gen += 1 # increment gen
        nextSubPopSizes = [eval(argList[9]), eval(argList[10])] # (SAA, MUU) growth until max
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed 
        overspillLethalEvent(pop, numRemovables)

    elif gen <= firstBalance:
        gen += 1 # increment gen
        nextSubPopSizes = [eval(argList[11]), eval(argList[12])] # (SAA, MUU) growth until max
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed 
        overspillLethalEvent(pop, numRemovables)

    
    elif gen < firstDecline:
        gen += 1
        nextSubPopSizes = [eval(argList[13]), eval(argList[14]), eval(argList[15])] #(SAA, SW, NE)
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < firstMinimum:
        gen += 1
        nextSubPopSizes = [eval(argList[16]), eval(argList[17]), eval(argList[18])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < secondBalance:
        gen += 1
        nextSubPopSizes = [eval(argList[19]), eval(argList[20]), eval(argList[21])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < firstRecovery:
        gen += 1
        nextSubPopSizes = [eval(argList[22]), eval(argList[23]), eval(argList[24])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < thirdBalance:
        gen += 1
        nextSubPopSizes = [eval(argList[25]), eval(argList[26]), eval(argList[27])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < secondDecline:
        gen += 1
        nextSubPopSizes = [eval(argList[28]), eval(argList[29]), eval(argList[30])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    else: #final rise - growth from bottleneck until present
        gen += 1
        nextSubPopSizes = [eval(argList[31]), eval(argList[32]), eval(argList[33])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    return True

def sample(pop): # samples the current population.
                 # we can add any functionality we want concerning the population (random sampling, freq counting, ...)
                 # then we can simply call this function with pyOperator in the middle of the simulation:

    # at this point we sample only males, so females are removed before sampling:
    smpPop = pop.clone()
    removables = list() # list of indices of individuals to be removed
    
    for idx, ind in enumerate(smpPop.individuals()):
        if ind.sex() == FEMALE:
            removables.append(idx) # this individual's index is added into removables
    smpPop.removeIndividuals(removables)

    # we would like sampleSize to be 1000 but sometimes that much isn't left in the sample population.
    # this fixes it:
    sampleSizes = list()
    for i in range(0, pop.numSubPop()):
        if smpPop.subPopSize(i) < 1000:
            sampleSizes.append(smpPop.subPopSize(i))
        else:
            sampleSizes.append(1000) # if 1000+ individuals in subpop, samplesize is 1000

    samples = drawRandomSamples(smpPop, sizes=sampleSizes, numOfSamples=10) # about the same as in earlier simulations
    myDump(samples, sampf, pop.dvars().gen)
    frequencies(samples, pop.dvars().gen, str(range(0,631)), typef, 0, 0, "seq") # calls frequencies for DNA-sequences
    #frequencies(samples, pop.dvars().gen, range(blah, blah), group_file) # calls frequencies for locus 1 (haplogroups)
    frequencies(samples, pop.dvars().gen, str(range(0,16)), Y_typef, 1, 2, "YSTR") # calls frequencies for YSTR-data
    YSTRArlequinDump(samples, YSTR_arlequinf, pop.dvars().gen) # calls Arlequin dumper for microsatellites
    freqArlequinDump(samples, mt_arlequinf, pop.dvars().gen) # calls Arlequin dumper for mt_dna
    seqArlequinDump(samples, seq_arlequinf, pop.dvars().gen)

    return True

def sequenceMap(sequence): # The goal here is to check whether the DNA-sequence given as first parameter
                                # already exists in global seqList - and if not, this new sequence is appended to the list.
                                # Thus, the 600-base-long sequences become associated with their absolute indices in the list:

    global seqList
    
    if len(seqList) == 0: # if the list is empty
            seqList.append(sequence)
            return 0
        
    if sequence in seqList:
        return seqList.index(sequence)
    else:
        seqList.append(sequence)
        return (len(seqList)-1) # gives back index (=haplotype callsign) of this new sequence

def seqArlequinDump(samples, target_filename, gen):

    global seqList

    arl_file = open(target_filename, "a")

    linelist = list()

    for i in range(0, len(samples)):
        
        linelist.append('[Profile]\n\n')
        linelist.append('\tTitle="mtDNA data (gen %d)- sample %d"\n\n' % (gen,i))
        linelist.append('\t\tNbSamples=%d\n' % 1)
        linelist.append('\t\tGenotypicData=%d\n' % 0)
        linelist.append('\t\tDataType=DNA\n')
        linelist.append('\t\tLocusSeparator=NONE\n')

        linelist.append('[Data]\n\n')
        linelist.append('\t[[Samples]]\n\n')

        for j in range(0, samples[i].numSubPop()):

            sampleSize = samples[i].subPopSize(j)
            
            linelist.append('\tSampleName=\"Subpopulation %d\"\n' % j)
            linelist.append('\tSampleSize=%d\n' % sampleSize)
            linelist.append('\tSampleData={\n')

            sampFreqs = list()
                              
            for idx, ind in enumerate(samples[i].individuals(j)): 
                appendCount(sampFreqs, range(631), 0, 0, ind)
            # then retrieve data from sampFreqs:

            for element in sampFreqs:
                haplotype = getX(element)
                count = getY(element) # gets the first object of element  
                linelist.append('%d\t%d\t' % (haplotype, count))

                for base in seqList[haplotype]:
                    if base == 0:
                        linelist.append("A")
                    if base == 1:
                        linelist.append("C")
                    if base == 2:
                        linelist.append("G")
                    if base == 3:
                        linelist.append("T")
                    if base == 4:
                        linelist.append("-")

                linelist.append('\n')

            linelist.append('\t\t}\n')

        # ACTIVATE IN SUBPOP PHASE:
        # Divide subpopulations into groups:
                                    
        #linelist.append(' [[Structure]]\n')
        #linelist.append('\tStructureName=\"ESA/LSA grouping\"\n')
        #linelist.append('\tNbGroups=%d\n' % 2)
        #linelist.append('\tGroup={\n')

        #linelist.append('\"Subpopulation %d\"\n' % 0)
        #linelist.append('\"Subpopulation %d\"\n' % 1)
        #linelist.append('\"Subpopulation %d\"\n' % 2)

        #linelist.append('}\n\n')
        #linelist.append('\tGroup={\n')

        #for j in range(3, samples[i].numSubPop()):
            #linelist.append('\"Subpopulation %d\"\n' % j)

        #linelist.append('}\n')
                                
    arl_file.writelines(linelist)
    arl_file.close()

    return True

def setX((x, y), new):
    return (new, y) # sets the first object of (x,y) to value 'new' and returns the pair

def setY((x, y), new):
    return (x, new) # sets the second object of (x,y) to value 'new' and returns the pair
# end tuple object control functions


def singletons(counters): # checks the number of singletons in list counters
                          # returns the number of singletons as an int value
    singTot = 0                      
    for k in counters: # for all counters in list
        if getY(k) == 1: # check whether the counter is showing 1...
            singTot = singTot + 1 # ...and if so, increment variable singTot
            
    return singTot # return the final value

# Note : This is 100% simuBottlenecks***.py -specific function with
#        no parameter or error checking. Please consult the authors in case
#        you have questions.
def statsCount(mt_file, ystr_file, stats_filename):

    global genList # function uses this
    
    # BEGIN PART ONE: Open MITOCHONDRIAL HAPLOTYPE frequencies data file for reading, read data into lists,
    #                 compute data using the lists, print data into output_file
    #=========================================================================================

    # open haplotype frequencies file for reading:
    input_file = open(mt_file, "r")

    #read all lines of this file into one list:
    readlist = input_file.readlines()

    mt_alleles = list() # list for saving number of alleles from every generation. format (generation,count) 
    mt_singletons = list()
    mt_diversities = list()

    mt_alleles_total = list()
    mt_singletons_total = list()
    mt_diversities_total = list()

    genc = -1 # generation counter - "-1" is the first zero value. Later reset as "0" 
    sp_counter = -1
    final = False
    # following for loop walks through the lines of the input file. Data of interest is stored into lists to be used later
    for string in readlist:

            if "SAMPLES @ GEN" in string:
                genc = genc + 1
             
            if "SAMPLES @ GEN 0" in string:
                sp_MAX = 2

            if "SAMPLES @ GEN 200" in string:
                sp_MAX = 2

            if "SAMPLES @ GEN 400" in string:
                sp_MAX = 2
                
            if "SAMPLES @ GEN 401" in string:
                sp_MAX = 3

            if "SAMPLES @ GEN 525" in string:
                sp_MAX = 3

            if "SAMPLES @ GEN 691" in string:
                sp_MAX = 3

            if "SAMPLES @ GEN 900" in string:
                sp_MAX = 3

            if "SAMPLES @ GEN 960" in string:
                sp_MAX = 3
                
            if "SAMPLES @ GEN 1099" in string:
                sp_MAX = 3
                final = True
                
            if "Subpopulation" in string: # We use this kind of if-clauses for distinguishing where we are going in the file
                sp_counter = sp_counter + 1
                if sp_counter == sp_MAX or final == True:
                    sp_counter = 0
                    final = False # we don't need final any more

            if "Different haplotypes in subpop" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                mt_alleles.append((sp_counter, genc, count))
                
            if "Singleton haplotypes in subpop" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                mt_singletons.append((sp_counter, genc, count))

            if "Nei\'s gene diversity for haplotypes in subpop" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                mt_diversities.append((sp_counter, genc, count))

            if "Different haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                mt_alleles_total.append((count))
                
            if "Singleton haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                mt_singletons_total.append((count))

            if "Nei\'s gene diversity for haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                mt_diversities_total.append((count))
 
            
    input_file.close()

    # this block computes the numbers we need and puts them into lists conveniently
    # so later all we need to do is just access the lists and print them

    mt_allele_averages = list()
    mt_singleton_averages = list()
    mt_diversity_averages = list()

    appendAverages(mt_alleles, mt_allele_averages) #walks through list and appends the contents
    appendAverages(mt_singletons, mt_singleton_averages) #etc. etc.
    appendAverages(mt_diversities, mt_diversity_averages)

    mt_allTot_averages = list()
    mt_singTot_averages = list()
    mt_divTot_averages = list()

    appendAveragesTotal(mt_alleles_total, mt_allTot_averages)
    appendAveragesTotal(mt_singletons_total, mt_singTot_averages)
    appendAveragesTotal(mt_diversities_total, mt_divTot_averages)

    # END PART ONE 

    # BEGIN PART TWO: Open Y_CHROMOSOMAL haplotype frequencies data file for reading, read data into lists,
    #                   compute data using the lists, print ALL data into stats file
    #=========================================================================================

    # open haplotype frequencies file for reading:
    input_file = open(ystr_file, "r")

    #read all lines of this file into one list:
    readlist = input_file.readlines()

    Y_alleles = list() # list for saving number of alleles from every subpop. format (sp,count) where sp=subPop & count=count 
    Y_singletons = list()
    Y_diversities = list()

    Y_alleles_total = list()
    Y_singletons_total = list()
    Y_diversities_total = list()

    genc = -1 # generation counter - "-1" is the first zero value. Later reset as "0" 
    sp_counter = -1
    final = False
    # following for loop walks through the lines of the input file. Data of interest is stored into lists to be used later
    for string in readlist:

            if "SAMPLES @ GEN" in string:
                genc = genc + 1
             
            if "SAMPLES @ GEN 0" in string:
                sp_MAX = 2

            if "SAMPLES @ GEN 200" in string:
                sp_MAX = 2

            if "SAMPLES @ GEN 400" in string:
                sp_MAX = 2
                
            if "SAMPLES @ GEN 401" in string:
                sp_MAX = 3

            if "SAMPLES @ GEN 525" in string:
                sp_MAX = 3

            if "SAMPLES @ GEN 691" in string:
                sp_MAX = 3

            if "SAMPLES @ GEN 900" in string:
                sp_MAX = 3

            if "SAMPLES @ GEN 960" in string:
                sp_MAX = 3
                
            if "SAMPLES @ GEN 1099" in string:
                sp_MAX = 3
                final = True
   
            if "Subpopulation" in string: # We use this kind of if-clauses for distinguishing where we are going in the file
                sp_counter = sp_counter + 1
                if sp_counter == sp_MAX or final == True:
                    sp_counter = 0
                    final = False

            if "Different haplotypes in subpop" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                Y_alleles.append((sp_counter, genc, count))
                
            if "Singleton haplotypes in subpop" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                Y_singletons.append((sp_counter, genc, count))

            if "Nei\'s gene diversity for haplotypes in subpop" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                Y_diversities.append((sp_counter, genc, count))

            if "Different haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                Y_alleles_total.append((count))
                
            if "Singleton haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                Y_singletons_total.append((count))

            if "Nei\'s gene diversity for haplotypes in sample" in string:
                chars = readNumbers(string, 0)
                count = float(makeString(chars)) # makes a string from chars and floats it
                Y_diversities_total.append((count))

                
    input_file.close()

    # this block computes the numbers we need and puts them into lists conveniently
    # so later all we need to do is just access the lists and print them

    Y_allele_averages = list()
    Y_singleton_averages = list()
    Y_diversity_averages = list()

    appendAverages(Y_alleles, Y_allele_averages) #walks through list 'alleles' and appends the contents into 'allele_averages'
    appendAverages(Y_singletons, Y_singleton_averages) #etc. etc.
    appendAverages(Y_diversities, Y_diversity_averages)

    Y_allTot_averages = list()
    Y_singTot_averages = list()
    Y_divTot_averages = list()
    
    appendAveragesTotal(Y_alleles_total, Y_allTot_averages)
    appendAveragesTotal(Y_singletons_total, Y_singTot_averages)
    appendAveragesTotal(Y_diversities_total, Y_divTot_averages)


    # So far so good. Then let's use the data and write it into output file (tab-separated values -style)

    #these two stupid-looking lines create output file for reading:
    output_file = open(stats_filename, "a")
    output_file.close()

    #these lines open the file again, this time for reading:
    output_file = open(stats_filename, "r")
    checklist = output_file.readlines()
    output_file.close()

    if len(checklist) == 0: # when true the output file was opened first time and we need to print header row into first line of file:
        writelist = list() #new list for writable data
        output_file = open(stats_filename, "a")
        writelist.append("Run\t")
        
        for i in genList: # mt_alleles for all gens
            writelist.append("mt_nA_SAA_%d" % i + "\t")

        for i in genList: # mt_alleles for all gens
            if i in [0,200,400]: # MUU for the first 3 samples
                writelist.append("mt_nA_MUU_%d" % i + "\t")

            else:
                writelist.append("mt_nA_SW_%d" % i + "\t")

        for i in genList: # mt_alleles for all gens

            if i not in [0,200,400]: # if not 2-subpop phase
                writelist.append("mt_nA_NE_%d" % i + "\t")

        for i in genList: # mt_alleles for all gens
            writelist.append("mt_nA_total_%d" % i + "\t")

        for i in genList: # mt_singletons for all gens
            writelist.append("mt_nS_SAA_%d" % i + "\t")

        for i in genList: # mt_singletons for all gens
            if i in [0,200,400]:
                writelist.append("mt_nS_MUU_%d" % i + "\t")
                
            else:
                writelist.append("mt_nS_SW_%d" % i + "\t")

        for i in genList: # mt_singletons for all gens
            if i not in [0,200,400]:
                writelist.append("mt_nS_NE_%d" % i + "\t")

        for i in genList: # mt_singletons for all gens
            writelist.append("mt_nS_total_%d" % i + "\t")

        for i in genList: # mt_GD for all gens
            writelist.append("mt_GD_SAA_%d" % i + "\t")

        for i in genList: # mt_GD for all gens
            if i in [0,200,400]:
                writelist.append("mt_GD_MUU_%d" % i + "\t")
            else:
                writelist.append("mt_GD_SW_%d" % i + "\t")

        for i in genList: # mt_GD for all gens
            if i not in [0,200,400]:
                writelist.append("mt_GD_NE_%d" % i + "\t")

        for i in genList: # mt_GD for all gens
            writelist.append("mt_GD_total_%d" % i + "\t")

        for i in genList: # Same for YSTR....
            writelist.append("Y_nA_SAA_%d" % i + "\t")

        for i in genList:
            if i in [0,200,400]:
                writelist.append("Y_nA_MUU_%d" % i + "\t")

            else:
                writelist.append("Y_nA_SW_%d" % i + "\t")

        for i in genList:
            if i not in [0,200,400]:
                writelist.append("Y_nA_NE_%d" % i + "\t")

        for i in genList:
            writelist.append("Y_nA_total_%d" % i + "\t")
            
        for i in genList:
            writelist.append("Y_nS_SAA_%d" % i + "\t")

        for i in genList:
            if i in [0,200,400]:
                writelist.append("Y_nS_MUU_%d" % i + "\t")
            else:
                writelist.append("Y_nS_SW_%d" % i + "\t")

        for i in genList:
            if i not in [0,200,400]:
                writelist.append("Y_nS_NE_%d" % i + "\t")

        for i in genList:
            writelist.append("Y_nS_total_%d" % i + "\t")
            
        for i in genList:
            writelist.append("Y_GD_SAA_%d" % i + "\t")

        for i in genList:
            if i in [0,200,400]:
                writelist.append("Y_GD_MUU_%d" % i + "\t")
            else:
                writelist.append("Y_GD_SW_%d" % i + "\t")

        for i in genList:
            if i not in [0,200,400]:
                writelist.append("Y_GD_NE_%d" % i + "\t")
            
        for i in genList:
            if i == 1099: # if we are looking at the last generation
                writelist.append("Y_GD_total_%d" % i + "\n") # put line feed to prevent extra TAB at the end of the line
            else:
                writelist.append("Y_GD_total_%d" % i + "\t")

        output_file.writelines(writelist) # writes header into output file
        output_file.close()

    writelist = list() # writelist is reset
    output_file = open(stats_filename, "a")
    writelist.append("%s\t" % (sys.argv[1])) # this appends the run number into the beginning of a data row

    # FOLLOWING for-LOOPS APPEND THE CRUCIAL NUMBERS INTO writelist (order is very important) 
    for i in range(0, len(mt_allele_averages)):
        writelist.append(str(mt_allele_averages[i]) + "\t")
        
    for i in range(0,len(mt_allTot_averages)):
        writelist.append(str(mt_allTot_averages[i]) + "\t")

    for i in range(0, len(mt_singleton_averages)):
        writelist.append(str(mt_singleton_averages[i]) + "\t")

    for i in range(0,len(mt_singTot_averages)):
        writelist.append(str(mt_singTot_averages[i]) + "\t")
        
    for i in range(0, len(mt_diversity_averages)): # the diversities are round to 6 decimals
        writelist.append("%f" %(mt_diversity_averages[i]) + "\t")

    for i in range(0,len(mt_divTot_averages)):
        writelist.append("%f" %(mt_divTot_averages[i]) + "\t")

    for i in range(0,len(Y_allele_averages)):
        writelist.append(str(Y_allele_averages[i]) + "\t")

    for i in range(0,len(Y_allTot_averages)):
        writelist.append(str(Y_allTot_averages[i]) + "\t")

    for i in range(0,len(Y_singleton_averages)):
        writelist.append(str(Y_singleton_averages[i]) + "\t")

    for i in range(0,len(Y_singTot_averages)):
        writelist.append(str(Y_singTot_averages[i]) + "\t")

    for i in range(0,len(Y_diversity_averages)):
        writelist.append("%f" % (Y_diversity_averages[i]) + "\t")
        
    for i in range(0,len(Y_divTot_averages)): # the diversities are round to 6 decimals
        if i == (len(Y_divTot_averages) - 1):
            writelist.append("%f" %(Y_divTot_averages[i]) + "\n")
        else:
            writelist.append("%f" %(Y_divTot_averages[i]) + "\t")
            
    output_file.writelines(writelist)

    output_file.close()
    
    return True

# function that dumps microsatellite data as one chunk in arlequin format
# NOTE: data produced by this function does not fit Arlequin until fileChopper_bn.py has been executed on arl_file

def YSTRArlequinDump(samples, arl_filename, gen):

    output_file = open(arl_filename, "a")

    linelist = list()

    for i in range(0, len(samples)):
        
        linelist.append('[Profile]\n\n')
        linelist.append('\tTitle="Microsatellites (gen %d) - sample %d"\n\n' % (gen, i))
        linelist.append('\t\tNbSamples=%d\n' % 1)
        linelist.append('\t\tDataType=MICROSAT\n')
        linelist.append('\t\tGenotypicData=%d\n' % 0)
        linelist.append('\t\tLocusSeparator=TAB\n')

        linelist.append('[Data]\n\n')
        linelist.append('\t[[Samples]]\n\n')

        sampleSize = samples[i].subPopSize()           

        linelist.append('\t\tSampleName="Microsatellites"\n')
        linelist.append('\t\tSampleSize=%d\n' % sampleSize)
        linelist.append('\t\tSampleData={\n')

        for idx,ind in enumerate(samples[i].individuals()):
            linelist.append(str(idx) + "\t1\t") # hard-coding the value 1 will suffice at this point
            for locus in range(0,16):
                linelist.append(str(ind.allele(locus, ploidy=1, chrom=2)) + "\t")
            linelist.append("\n") # line feed

        linelist.append('}\n')

    output_file.writelines(linelist)
    output_file.close
        
#END USER DEFINED FUNCTIONS

#Here's our population:
pop = Population(size=[250, 250],
                 ploidy=2,
                 loci=[631, 0, 16], # 631 loci in mitochondrial DNA, one for each base - 16 loci in Y-chromosome
                 chromTypes=[CUSTOMIZED, CHROMOSOME_X, CHROMOSOME_Y],
                 infoFields=['age', 'group', 'migrate_to'], ancGen=0
                 )

# set virtual subpopulations, male and female (necessary for sex-specific migration)
initialSplitter = ProductSplitter(splitters=[SexSplitter(), InfoSplitter('group', values=[0, 1, 2, 3])])
interSplitter = CombinedSplitter(splitters=[initialSplitter],
                                 vspMap=[(0, 4), (1, 5), (2, 6), (3, 7)])

finalSplitter = CombinedSplitter(splitters=[interSplitter, SexSplitter()])

pop.setVirtualSplitter(finalSplitter)

# The virtual subpops should now have been defined as follows (where x is [non-virtual] subpop index) :

# (x, 0) : Male, group = 0 or Female, group = 0 (indicating males and females that are not in mating age)
# (x, 1) : Male, group = 1 or Female, group = 1 (males and females in mating age)
# (x, 2) :
# (x, 3) :
# (x, 4) : Male (all male individuals of subpop x)
# (x, 5) : Female (all female indivduals of subpop x)

pop.dvars().gen = 0 # intialize the generation number

# Mutation constants for MatrixMutator:
u = float(argList[6])
kappa = float(argList[7])

# Evolve the POPULATION (self-evolving populations are a simuPOP 1.0.0+ feature)
# ----------------------
pop.evolve(

    initOps = [PyOperator(func=initPop),
               PyOperator(func=initGroups),
               PyOperator(func=sample),
               PyOperator(func=printVsp)],
    
    preOps = [
	# split into two subpops at generation 401 (7000 BP)
	SplitSubPops(subPops=[1], proportions=[0.333,1-0.333], at=401),
        PyOperator(func=loadBgPops),
        InfoExec('age += 1', begin=2),
        PyOperator(func=cleanUp, begin=2),

        # Migration from neighbouring populations in the 2-subpop-phase:
        PyOperator(insertMigration, param=[0], begin=0, end=400),
        PyOperator(insertMigration, param=[1], begin=0, end=400),

        # Migration from neighbouring populations in the 3-subpop-phase:
        PyOperator(insertMigration, param=[0], begin=401),
        PyOperator(insertMigration, param=[1], begin=401),
        PyOperator(insertMigration, param=[2], begin=401),

        # Internal migration (2-subpop-phase):
        Migrator(rate=[
            [0, malnear], # males from subpop 0 to subpops 0 and 1
            [0, femnear], # females from subpop 0 to subpops 0 and 1
            [malnear, 0], # males from subpop 1 to subpops 1 and 0
            [femnear, 0]  # females from subpop 1 to subpops 1 and 0
            ],
            mode=BY_PROPORTION,
            subPops=[(0,4), (0,5), (1,4), (1,5)],
            end=400
            ),

        # Internal migration (3-subpop-phase):
        Migrator(rate=[
            [0, malnear, malfar], # males from subpop 0 to subpops [0, 1, 2]
            [0, femnear, femfar], # females from subpop 0 to subpops [0, 1, 2]
            [malnear, 0, malnear], # ...
            [femnear, 0, femnear],  # ...
            [malfar, malnear, 0],
            [femfar, femnear, 0], # females from subpop 2 to subpops [0, 1, 2]
            ],
            mode=BY_PROPORTION,
            subPops=[(0,4), (0,5), (1,4), (1,5), (2,4), (2,5)],
            begin=401
            ),

        PyOperator(func=removeOverspill), # randomly reduce population size so that parent generation fits.
        PyOperator(func=printNumSubPops),
        
        # mutation call for mitochondrial DNA sequences. Takes deletions into account too:
        MatrixMutator(rate=[
            [0, u/4, (u/4)*kappa, u/4, 0],
            [u/4, 0, u/4, (u/4)*kappa, 0],
            [(u/4)*kappa, u/4, 0, u/4, 0],
            [u/4, (u/4)*kappa, u/4, 0, 0],
            [0, 0, 0, 0, 1]
            ],
                      loci = range(631)),
        
        # mutation call for Y_chromosomal microsatellites
        StepwiseMutator(rates=float(argList[8]), loci=range(631,647)), # mutation rates equal
        
        # debug printing:
        #PyEval(r"'Generation: %d\n' %(gen) "),
        PyOperator(func=printSize),

        PyOperator(func=sample, at=[firstCheck, secondCheck, thirdCheck, fourthCheck, fifthCheck, sixthCheck, seventhCheck, eighthCheck, ninthCheck]),
        PyOperator(func=dumpSequences, at=ninthCheck),
        ],

    # mating scheme similar to that of Savonian expansion simulations:
    matingScheme = HeteroMating(matingSchemes=[
                ## BO: I see why you have to use setMatingScheme, which is not necessary when you use
                ## the new ALL_AVAIL feature in 1.0.1
                ##
                ## Juhana: setMatingScheme has been removed from script

                 # age <= maxAge, copy to the next generation (weight=-1).
                 # no need to worry about too old individuals, they're removed before mating.
                 CloneMating(subPops=[(ALL_AVAIL, 0), (ALL_AVAIL, 1), (ALL_AVAIL, 2)],
                             ops=[
                             ## BO: because CloneGenoTransmitter does not
                             ## by default handle customized chromosomes, an
                             ## explicit list of chromosomes have to be used.
                             CloneGenoTransmitter(chroms=[0,2])
                             ],
                ## BO: I do not see weight = -1, this is important because other wise the number of cloned
                ## individuals will be determined by relative number of individuals compared to randomMating guys.
                             weight=-1),

                 # random mating for individuals in mating ages (vsp 1 in every subpop)
                 RandomMating(subPops=[(ALL_AVAIL, 1)],
                              ops=[MitochondrialGenoTransmitter(),
                                   MendelianGenoTransmitter(),
                                   #IdTagger(), PedigreeTagger(), ParentsTagger()],
                                   ],
                              numOffspring=(POISSON_DISTRIBUTION, 2))
                 ]
                , subPopSize=demo),

    postOps=[
	TicToc(at=[0,1099], output=">>>serial_running_times.txt"),
        PyOperator(func=lethalEvent),
        Stat(popSize=True, subPops=[(ALL_AVAIL, 0), (ALL_AVAIL, 1), (ALL_AVAIL, 2), (ALL_AVAIL, 3)]),
        #PyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in subPopSize]))")
        ],
    
    finalOps = [PyOutput('Simulation complete\n'), Dumper(max=10, structure=True)],
    
    gen = numGens
    )

# finally count stats of interest and perform mySampler on the output files:
statsCount(typef, Y_typef, run_statsf)
mySampler(sampf, YSTR_outputf, 907, "YSTR")
mySampler(sampf, mt_outputf, 832, "mt")
myStatsCount(mt_outputf, YSTR_outputf, 832, 907, mySampler_statsf)
