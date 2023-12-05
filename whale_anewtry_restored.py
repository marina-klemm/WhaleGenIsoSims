#This script uses simuPOP to simulate the evolution of Aotearoa New Zealand's
#population of tohorā Southern Right Whales and its recovery from a very
#harsh bottleneck.
#Code originally written by Dr Emma Carroll (EC) and later edited by
#Dr Marina Klemm (MK). This code also uses some chunks shared by other
#simuPOP users, cited where needed.

#Contact: marinaklemm@gmail.com
#Github: marina-klemm


# =============================================================================
# Loading all packages
# =============================================================================
import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from simuPOP.utils import export
import random
import numpy as np
from operator import add
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP.utils import Trajectory, simulateForwardTrajectory, export, Exporter
import simuPOP.demography as demo
from simuPOP.utils import viewVars
import matplotlib
import csv
import logging


# =============================================================================
# MK: Population numbers used as a reference here
# =============================================================================
# Use Jackson et al, 2016 to estimate trajectory in 20 years intervals
# Starting in 1830 and ending in 2030: (divided into 25 years)
    
#1830-1855: ~ 30,000 individuals
#1855-1880: ~ 3,000
#1880-1905: ~ 200
#1905-1930: ~ 150
#1930-1955: ~ 40
#1955-1980: ~ 150
#1980-2005: ~ 1000 
#2005-2030: ~ 2,139 (2009)
#2010-2030: ~ 2500 (maybe, at least?)

# Ways of thinking/calculating generations:
    # Gen 1: 1830-1855 maximum amount of whales
    # Gen 5: 1930-1955 lowest amount (bottleneck)
    # Gen 7: 1980-2005 recovering 
    # Gen 10: 2010-2030 current(ish) generation
    
# Or:
    # generations: 100,000 years/25 = 4,000
    # Gen 4,000: 1830-1855 maximum amount of whales
    # Gen 4,005: 1930-1955 lowest amount (bottleneck)
    # Gen 4,007: 1980-2005 recovering 
    # Gen 4,010: 2010-2030 current(ish) generation
    
# To allow trajectory to be calculated, it has to have a number of generations 
#between 80 and 1000, from my tries. So, I will start at 75; the first 75 generations
#will not expand, it will be a "burn in" time.

#Gen 75: 1830-1855: ~ 30,000 individuals
#Gen 76: 1855-1880: ~ 3,000
#Gen 77: 1880-1905: ~ 200
#Gen 78: 1905-1930: ~ 150
#Gen 79: 1930-1955: ~ 40
#Gen 80: 1955-1980: ~ 150
#Gen 81: 1980-2005: ~ 1000 
#Gen 82: 2005-2030: ~ 2,139 (2009)
#Gen 83: 2010-2030: ~ 2500 (maybe, at least?)



#FIXME
# Problem here:
    #after exporting each generation into a different output file, I could see
    #that there was no change in the number of individuals, a.k.a., every single
    #generation had 30,000 individuals.
    #So, the model needs to be improved, maybe with dumping of dead individuals?
    
#FIXME
# Another thing:
    #some individuals with ind_id 0 show up after a few generations AND
    #individuals with no father_id or mother_id also show up across the evolution,
    #which is problematic since we cannot trust the mtDNA information then.
  

# =============================================================================
# MK Model 8
# All the other deleted models are in Removed_Models_Whale.py
# =============================================================================

#Gen 1000: 1830-1855: ~ 30,000 individuals
#Gen 1001: 1855-1880: ~ 3,000
#Gen 1002: 1880-1905: ~ 200
#Gen 1003: 1905-1930: ~ 150
#Gen 1004: 1930-1955: ~ 40
#Gen 1005: 1955-1980: ~ 150
#Gen 1006: 1980-2005: ~ 1000 
#Gen 1007: 2005-2030: ~ 2,139 (2009)
#Gen 1008: 2010-2030: ~ 2500 (maybe, at least?)

#Gen 75: 1830-1855: ~ 30,000 individuals
#Gen 76: 1855-1880: ~ 3,000
#Gen 77: 1880-1905: ~ 200
#Gen 78: 1905-1930: ~ 150
#Gen 79: 1930-1955: ~ 40
#Gen 80: 1955-1980: ~ 150
#Gen 81: 1980-2005: ~ 1000 
#Gen 82: 2005-2030: ~ 2,139 (2009)
#Gen 83: 2010-2030: ~ 2500 (maybe, at least?)


model8 = demo.MultiStageModel([
        demo.LinearGrowthModel(
            #up until staticPhaseEnd
        T=75,
        N0=[10000, 10000], 
        NT=[15000, 15000] 
        ),
       #Changing the instant change to a exponential change:
    demo.ExponentialGrowthModel(
    T=1, #first decline
    NT=[(1500), (1500)]  
    ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   NT=[(100), (100)]  
   ),
  demo.ExponentialGrowthModel(
  T=1, #1 generation
   NT=[(75), (75)] 
  ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
  # N0=[(75), (75)], 
   NT=[(20), (20)]
   ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   #N0=[(20), (20)], 
   NT=[(75), (75)]
   ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   #N0=[(75), (75)], 
   NT=[(500), (500)]
   ),
   demo.LinearGrowthModel(
   T=2, #2 generations
  # N0=[(500), (500)], 
   NT=[(1250), (1250)]
   )])
    
    # get a visual presentation of the demographic model
model8.plot('log/demoModel.png',
    title='Model 8: A simple model with more exponential changes and subpopulations')

model8.init_size #returns the initial population size
model8.info_fields
model8.num_gens 


# =============================================================================
# Model 9: trying to make sure the popSize changes over pop.evolve
# =============================================================================
#From the vignette:
#An instant population growth model that evolves a population from size N0 to NT 
#for T generations with population size changes at generation G to NT.
#If G is a list, multiple population size changes are allowed. In that case, a list 
#(or a nested list) of population size should be provided to parameter NT. 
#Both N0 and NT supports fixed (an integer), dynamic (keep passed poulation size) 
#and proportional (an float number) population size. Optionally, one or more operators 
#(e.g. a migrator) ops can be applied to population. Required information fields by 
#these operators should be passed to parameter infoFields. If removeEmpty option is set to True, 
#empty subpopulation will be removed. This option can be used to remove subpopulations.


model9 = demo.InstantChangeModel(T=83, #generations
                                 N0=[(15000, 15000)], #initial size
                                 G=[76, 77, 78, 79, 80, 81, 82], #generations where the size changes
                                 NG=[[(1500), (1500)],[(100), (100)], [(75), (75)],
                                     [(20), (20)], [(75), (75)], [(500), (500)],
                                     [(1050), (1050)]],
                                 #NG=[3000, 200, 150, 40, 150, 1000, 2139], 
                                 ops=[], 
                                 infoFields=[], 
                                 removeEmptySubPops=False)

model9.plot('log/demoModel.png',
    title='Model 9: A simple model with instant changes')

model9.init_size #returns the initial population size
model9.info_fields
model9.num_gens 


# =============================================================================
# Initial code (by EC, edited by MK, unneeded things removed)
# =============================================================================
## EC
## this code is to generate simulated whale populations that have mtDNA, SNP genotypes
## and stable isotope profiles. the code generates two subpopulations in sim that 
## correspond to whale wintering grounds 
## virtual subpopulations (VSPs) are used to represent feeding ground and age-sex classes

## set parameters for simulation - breeding ground which are subpopulations in simupop terms

minMatingAge = 6 ## minimum age at first reproduction
maxMatingAge = 50 ## max age of reproduction
gen = 82 ## for trajectory, based on the model6
nb_loci = 100 ## number of loci to simulate
scenario_id = "1"

## setting up the feeding ground variables
## mean (mean_) and variance (variance_) set for both C and N for two feeding grounds
## deviant proportion: proportion of males that will go to non-natal wintering 
# ground for one winter/breeding opportunity
mean_C = [16.7, 20.5] 
variance_C = [3.24, 2.89]
mean_N = [5.5, 8.7]
variance_N = [0.25, 0.49]
deviant_proportion = 0.1

## Sample count is number of samples taken per wintering ground
numberOfFeedingGrounds = 1
sample_count = 60

# Needed a way to ensure that the simulations begin with whales that are related 
# and show correlation between 
# SI and genetic data. Did this by rapidly expanding population, creating many offspring
# For the first 10 generations, we expand the next generation by 7% (this leads 
# to a rough doubling after 10 years).
# After the population_growing_period, each subpopulation size is kept constant
population_growth_rate = 1.07 
population_growing_period = 10 # in years

#sub_population_size = [15000, 15000] ## MK, same as size, changed for trajectory

def postop_processing(pop):
    for i in range(0, pop.numSubPop()):
        for individual in pop.individuals(i):
            # The feeding ground is fixed at birth (inherited from mother)
            # The C and N values are sampled from a distribution based on the 
            # feeding ground each year
            # The 'feeding_ground' info field is a float. We cannot use that as 
            # an array index so convert to an int
            feeding_ground = int(individual.info('feeding_ground'))
            individual.setInfo(np.random.normal(mean_C[feeding_ground], np.sqrt(variance_C[feeding_ground])), 'carbon')
            individual.setInfo(np.random.normal(mean_N[feeding_ground], np.sqrt(variance_N[feeding_ground])), 'nitrogen')

            # print("Individual ", individual.info('ind_id'), " has native breeding ground ", 
            # individual.info('native_breeding_ground'), " and is currently at breeding ground ", i)
            # Migration
            # Initially, set the migrate_to to the current population of the individual
            individual.setInfo(i, 'migrate_to')
            # If the individual is a male, then we can optionally migrate them using the migrate_to info field
            if individual.sex() == sim.MALE and individual.info('age') >= minMatingAge:
                # If the individual has already migrated, always move them back
                if individual.info('native_breeding_ground') != i:
                    #print("Moving individual ", individual.info('ind_id'), " back to their native breeding ground ", individual.info('native_breeding_ground'), " from temporary breeding ground ", i)
                    individual.setInfo(individual.info('native_breeding_ground'), 'migrate_to')
                # Otherwise, migrate them to another population with a probabilistic model
                elif np.random.uniform() < deviant_proportion:
                    # Individual will migrate.
                    new_population = (i + 1) % 2
                    #print("Individual ", individual.info('ind_id'), " will migrate to ", new_population)
                    individual.setInfo(new_population, 'migrate_to')
    return True

def init_native_breeding_grounds(pop):
    # Assign the native breeding ground to each individual. I don't know how to do this except by doing it individually
    # Fortunately, we can just inherit this maternally, so it only has to be run once
    for i in range(0, pop.numSubPop()):
        for individual in pop.individuals(i):
            individual.setInfo(i, 'native_breeding_ground');
    return True

def configure_new_population_size(gen, pop):
    # It is critical to specify the sub population sizes independently of each other. Each sub-population may be a different size
    if (gen < population_growing_period):
        return [pop.subPopSize(0) * population_growth_rate, pop.subPopSize(1) * population_growth_rate]
    else:
        return [pop.subPopSize(0), pop.subPopSize(1)]

def runSimulation(scenario_id, sub_population_size, minMatingAge, maxMatingAge, gen):
    '''
    sub_population_size   A vector giving the population sizes for each sub-population. The subpopulations determine which breeding ground an individual belongs to
    minMatingAge          minimal mating age.
    maxMatingAge          maximal mating age. Individuals older than this are effectively dead
    years                 number of years to simulate
    '''

# scenario_id describes the batch of files to load
# The mitochondrial DNA will be in mtdna_<scenario_id>
# The SNP DNA will be in snp_<scenario_id>

# Read the mitochondrial haplotype frequencies. There's a bit to unpack here
# We read the lines into an array, and for each one, call split() on it to get one element per column.
# However, we do not want this - we want the transpose, where haplotype_frequencies[0] is a vector of
# all the frequencies for population 0, and haplotype_frequencies[1] is the corresponding vector for
# population 2. list(map(list, zip(*t))) will achieve this transformation for us.
# While we are at it, we also convert the strings into floats.
mitochondrial_file = "mtdna_" + scenario_id + ".txt"
with open(mitochondrial_file, "r") as fd:
    haplotype_frequencies = list(map(list, zip(*[list(map(float, line[0:-1].split())) for line in fd])))

#if len(haplotype_frequencies) != len(sub_population_size):
 #   raise ValueError('The number of populations in the population size vector and the number of populations deduced from the haplotype file are different')

if len(haplotype_frequencies) != len(model8.init_size):
    raise ValueError('The number of populations in the population size vector and the number of populations deduced from the haplotype file are different')


# Now read the SNP data. This builds a 2D array indexed as snp[locus][population]
snp_file = "snp_" + scenario_id + ".txt"
with open(snp_file, "r") as fd:
    snp = [list(map(float, line[0:-1].split())) for line in fd]

sub_population_count = len(model8.init_size)
print(sub_population_count, "subpopulations detected")

# Now we can create the population. We want to give each population a population name, 
#starting from A
sub_population_names = list(map(chr, range(65, 65+sub_population_count)))



# =============================================================================
# MK: Trajectory
# =============================================================================

traj = simulateForwardTrajectory(N=model9.init_size, 
                                 beginGen = 0,
                                 endGen = model9.num_gens, 
                                 beginFreq = [0, 1], 
                                 endFreq = [[0, 1], [0, 1]],
                                 fitness = None)

traj.func() #This is to test that traj.func() exists, which is needed for the
# pop.evolve function. I found that less than 80 generations or more than 6,000
# individuals cause the function to crash, with the message:
# Error message: Cell In[119], line 1
    #traj.func()
#AttributeError: 'NoneType' object has no attribute 'func'


#Values that worked:
    #N=[150, 150], beginGen = 75, endGen = 100, beginFreq = [0.5, 0.5], endFreq = [[0.1, 0.2], [0.2, 0.3]]
    #N=[1000,1000], beginGen = 0, endGen = 87, beginFreq = [0.2, 0.3], endFreq = [[0.1, 0.11], [0.2, 0.21]]
    #N=[15000, 15000], beginGen = 0, endGen = 87, beginFreq = [0.1, 0.2], endFreq = [[0.1, 0.2], [0.1, 0.2]]
    #N=[15000, 15000], beginGen = 0, endGen = 4000, beginFreq = [0.1, 0.2], endFreq = [[0.1, 0.2], [0.1, 0.2]]
    #N=[10000, 10000], beginGen = 0, endGen = 83, beginFreq = [0.2, 0.3], endFreq = [[0.1, 0.2], [0.1, 0.2]]
    
#Now that it worked, I can use the original code and add freqFunc=traj.func() into the matingScheme options.

pop = sim.Population(size = model9.init_size, 
                             ploidy=2,
                             loci=[nb_loci, 1],
                             ancGen=2,#Number of the most recent ancestral generations 
                             #to keep during evolution, i.e., ancGen=2 keep parental and
                             #grandparental generations coexisting with the newest one.
                             infoFields=['age', 'ind_id', 'father_id', 'mother_id', 'nitrogen', 'carbon', 'feeding_ground', 'native_breeding_ground', 'migrate_to'],
                             subPopNames = sub_population_names,
                             chromTypes=[sim.AUTOSOME, sim.MITOCHONDRIAL])

sub_population_names = tuple(sub_population_names)



# Create an attribute on each individual called 'age'. Set it to a random number between 0 and maxMatingAge
# Note that size is a vector - the size of each population. We have to sum these to get the total number of individuals
#individual_count = sum(sub_population_size)
individual_count = sum(model8.init_size)

# Assign a random age to each individual
pop.setIndInfo([random.randint(0, maxMatingAge) for x in range(individual_count)], 'age')
# Assign a random feeding ground to each individual
pop.setIndInfo([random.randint(0, numberOfFeedingGrounds-1) for x in range(individual_count)], 'feeding_ground')


# Currently we have these virtual subpopulations:
# age < minMatingAge (juvenile)
# age >= minMatingAge and age < maxMatingAge + 0.1 (age <= maxMatingAge) (mature)
# age >= maxMatingAge (dead)
#
# Ideally we would want something like this:
# 1) Immature
# 2) Receptive female (every 3 years)
# 3) Non-receptive female
# 4) Mature male
# 5) Dead
# 
# Note that we use a cutoff InfoSplitter here, it is also possible to
# provide a list of values, each corresponding to a virtual subpopulation.
pop.setVirtualSplitter(sim.CombinedSplitter([
    sim.ProductSplitter([sim.SexSplitter(),
                             sim.InfoSplitter('age', cutoff=[minMatingAge, maxMatingAge + 0.1], names=['juvenile', 'mature', 'dead'])])],
                             vspMap = [[0], [1], [2], [3], [4], [5], [0, 1, 3, 4], [1,4]],
                              names = ['Juvenile Male', 'Mature Male', 'Dead Male', 'Juvenile Female', 'Mature Female', 'Dead Female', 'Not dead yet', 'Active']))

sim.dump(pop)

# =============================================================================
# MK: Check if VSPs are initiallized properly
# =============================================================================
pop.numVirtualSubPop()
pop.subPopName([0,0])
pop.subPopName([0,1])
pop.subPopName([1,7])
pop.subPopName([1,4])


##############################################################################
# =============================================================================
# Printing 10 alleles:
# =============================================================================


pop.evolve(
    initOps = [sim.InitSex(), sim.IdTagger(), sim.PyOperator(func=init_native_breeding_grounds)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
    # increase age by 1
    preOps = [sim.InfoExec('age += 1'),
 #Export population in each generation
   Exporter(format='csv', infoFields=('age', 'ind_id', 'father_id', 'mother_id', 'nitrogen', 'carbon', 'feeding_ground', 'native_breeding_ground', 'migrate_to'), 
           output="!'dump_gen_%d.csv' % gen", step=1, begin=75)
           ],
    matingScheme = sim.HeteroMating([
        # age <= maxAge, copy to the next generation (weight=-1)
        # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
        # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
        sim.CloneMating(subPopSize=model9, 
                        ops=[sim.CloneGenoTransmitter(chroms=[0,1])],
                            subPops=[(sub_population, 6) for sub_population in range(0, sub_population_count)], #MK:6 here is because
                            #two of the eight virtual subpopulations are deceased.
                            weight=-1), #MK: if weights are negative, they are multiplied to their parental subpopulation;
            #For example: if parental pop = (500, 1000), and weight = -2, next 
            #generation pop= (1000, 2000). 
            #For weight -1, it keeps the number of individuals from the parental generation.
            #ALSO: if there is a mix of negative and positive weights, the negative will be processed first.
        # Then we simulate random mating only in VSP 1 (ie reproductively mature individuals) within subpopulation (breeding/winter grounds)
        sim.RandomMating(subPopSize=model9,
            ops=[sim.MitochondrialGenoTransmitter(),
                                  sim.MendelianGenoTransmitter(),
                                  sim.IdTagger(),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['feeding_ground']),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['native_breeding_ground']),
                                  sim.PedigreeTagger()],
                             subPops=[(sub_population, 7) for sub_population in range(0, sub_population_count)],
                             weight=1),
        sim.ControlledRandomMating(subPopSize=model9,
                                   freqFunc=traj.func(),
                                   weight=1)] #MK: we decided to keep the same weight as the
        #mitochondrial transmitter.
        ),
    postOps = [

    # Determine the isotopic ratios in individuals
    sim.PyOperator(func=postop_processing),
    sim.Migrator(mode=sim.BY_IND_INFO),
        # count the individuals in each virtual subpopulation
        #sim.Stat(popSize=True, subPops=[(0,0), (0,1), (0,2), (1,0), (1, 1), (1, 2)]),
        # print virtual subpopulation sizes (there is no individual with age > maxAge after mating)
        #sim.PyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in subPopSize]))")

        # Alternatively, calculate the Fst
        # FIXME: How does this actually work? Does it work for > 2 populations? I don't really understand it yet
        # ELC: it is a calculation that partitions variance among and between populations, and can be calculated as a 
        # global statistic or on a pairwise basis. We use it as an indication of genetic differentiation.

        sim.Stat(structure=range(1), subPops=sub_population_names, suffix='_AB', step=10),
        sim.Stat(numOfMales=True, 
                 begin = 75, step = 1),     
        sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10), #Print Fst every 10 steps
        sim.Stat(alleleFreq=[1, 2, 3, 4, 5, 6, 7, 8, 9, 100], vars=['alleleFreq_sp'], step=10), #added this now, to
        #calculate the allele frequencies in selected loci
       sim.PyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % ("
           "subPop[0]['alleleFreq'][1][1], subPop[0]['alleleFreq'][2][1], subPop[0]['alleleFreq'][3][1],"
           "subPop[0]['alleleFreq'][4][1], subPop[0]['alleleFreq'][5][1], subPop[0]['alleleFreq'][6][1],"
           "subPop[0]['alleleFreq'][7][1], subPop[0]['alleleFreq'][8][1], subPop[0]['alleleFreq'][9][1],"
           "subPop[0]['alleleFreq'][100][1], subPop[1]['alleleFreq'][1][1], subPop[1]['alleleFreq'][2][1],"
           "subPop[1]['alleleFreq'][3][1], subPop[1]['alleleFreq'][4][1], subPop[1]['alleleFreq'][5][1],"
           "subPop[1]['alleleFreq'][6][1], subPop[1]['alleleFreq'][7][1], subPop[1]['alleleFreq'][8][1],"
           "subPop[1]['alleleFreq'][9][1], subPop[1]['alleleFreq'][100][1])", step=1, begin = 75),
      # sim.PyEval(r"'%.2f' % haploFreq[0]", step=10)

        #to print out the allele frequencies in selected loci
       
    ],
    finalOps= sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
    gen = model9.num_gens
    )

pop.vars()

# print out population size and allele frequency
for idx, name in enumerate(pop.subPopNames()):
    print('%s (%d): %.4f' % (name, pop.subPopSize(name), 
        pop.dvars(idx).alleleFreq[0][0]))
    
    

viewVars(pop.vars(), gui=False)

sim.dump(pop)


ped = sim.Pedigree(pop);
print("This is the pedigree stuff")


# Now sample the individuals
sample = drawRandomSample(pop, sizes=[sample_count]*sub_population_count)

# Print out the allele frequency data
sim.stat(sample, alleleFreq=sim.ALL_AVAIL)
frequencies = sample.dvars().alleleFreq;
with open('freq.txt', 'w') as freqfile:
    index = 0
    for locus in frequencies:
        if (locus == nb_loci):
            continue
        if (len(frequencies[locus]) < 2):
            continue
        print(index, end=' ', file=freqfile)
        index = index + 1
        for allele in frequencies[locus]:
            print(frequencies[locus][allele], end=' ', file=freqfile)
        print(file=freqfile)
        
# =============================================================================
# MK: Print out mtDNA frequency data in a way that is comparable to the original
# file.        
# =============================================================================
# Print out the mtDNA frequency data
samplemtDNA_count = 6000 #to grab all the possible final haplotypes; this is limited
#by te highest number of individuals in a subpopulation, which is 6,053 in this case
#for subpopulation 1. Once the number of inidividuals is corrected based on the bottleneck,
#this value here has to be changed.
samplemtDNA = drawRandomSample(pop, sizes = [samplemtDNA_count]*sub_population_count)
sim.stat(pop, alleleFreq=sim.ALL_AVAIL)

frequencies =samplemtDNA.dvars().alleleFreq;
last_locus = nb_loci  # Index of the last locus

with open('freq_mtDNA.txt', 'w') as freqfile:
    index = 0
    for locus in frequencies:
        if locus != last_locus:
            continue

        if len(frequencies[locus]) < 2:
            continue

        print(index, end=' ', file=freqfile)
        index += 1

        for allele in frequencies[locus]:
            print(frequencies[locus][allele], end=' ', file=freqfile)

        print(file=freqfile)

#The output has 56 values that add up to 1, which is what is expected for a 
#mtDNA haplotype frequency data.

# We want to remove monoallelic loci. This means a position in the genotype for which all individuals have the same value in both alleles
# To implement this we will build up a list of loci that get ignored when we dump out the file. Generally speaking, if we add all the values up
# then either they will sum to 0 (if all individuals have type 0) or to the number of individuals * 2 (if all individuals have type 1)
geno_sum = [0] * (nb_loci + 1) * 2;
for individual in sample.individuals():
    geno_sum = list(map(add, geno_sum, individual.genotype()))
final_sum = list(map(add, geno_sum[:(nb_loci+1)], geno_sum[(nb_loci+1):]))

monoallelic_loci = [];
for i in range(0, nb_loci):
    if final_sum[i] == 0 or final_sum[i] == sample_count*sub_population_count*2:
        monoallelic_loci = [i] + monoallelic_loci
monoallelic_loci = sorted(monoallelic_loci, reverse=True)

nb_ignored_loci = len(monoallelic_loci)
# Generate the two files
with open('mixfile.txt', 'w') as mixfile:
    with open('haploiso.txt', 'w') as haplofile:
        print(sub_population_count, nb_loci - nb_ignored_loci, 2, 1, file=mixfile)
        print("sex, haplotype, carbon, nitrogen, native_ground", file=haplofile);
        for i in range(0, nb_loci - nb_ignored_loci):
            print('Loc', i+1, sep='_', file=mixfile);
        for individual in sample.individuals():
            genotype = individual.genotype();
            print(1 if individual.sex() == 1 else 0,
                  genotype[nb_loci],
                  individual.info('carbon'),
                  individual.info('nitrogen'),
                      int(individual.info('native_breeding_ground')),
                  file=haplofile, sep=' ')
            print(int(individual.info('native_breeding_ground')+1), end=' ', file=mixfile)
            for i in range(0, nb_loci):
                if i not in monoallelic_loci:
                    print(genotype[i]+1, genotype[i+nb_loci+1]+1, ' ', end='', sep='', file=mixfile)
            print(file=mixfile);
#return sample #MK: return outside function











