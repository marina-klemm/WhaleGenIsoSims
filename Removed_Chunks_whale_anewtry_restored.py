# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:28:45 2023

@author: marin
"""
# =============================================================================
# Using mtDNA_1.txt as the initial frequency
# =============================================================================


#Coding for beginFreq
#From utils.py:
# beginFreq
#            The initial allele frequency of involved loci in all subpopulations.
#            It can be a number (same frequency for all loci in all
#            subpopulations), or a list of frequencies for each locus (same
#            frequency in all subpopulations), or a list of frequencies for each
#            locus in each subpopulation in the order of ``loc0_sp0``,
#            ``loc1_sp0``, ..., ``loc0_sp1``, ``loc1_sp1``, ... and so on.

#So, beginFreq needs a pair of frequencies (2 loci for all alleles), or all the allele
#frequencies for all the loci, or even all those pairs of frequencies for every 
#subpopulation.

#Trying to achieve that for the mtDNA_1.txt frequencies:
file_path = "mtDNA_1.txt"
haplotype_freq_data = np.loadtxt(file_path, delimiter='\t')

# Extract frequencies for each subpopulation
subpop1_freq = haplotype_freq_data[:, 0]  # first column is for subpopulation 0
subpop2_freq = haplotype_freq_data[:, 1]  # second column is for subpopulation 1


subpop1_freq = subpop1_freq.tolist()
len(subpop1_freq) #I have 58 points, but 100 loci, so I will randomly select 42
#of them to add to this simulation

# Select 42 additional datapoints randomly
additional_datapoints = np.random.choice(subpop1_freq, size=42, replace=False)

# Combine the original and additional datapoints to have a total of 100
combined_datapoints = subpop1_freq + additional_datapoints.tolist()
len(combined_datapoints) #100
nb_loci

#Since this is a haplotype list, I need to add a frequency = 0 to represent
#the second allele

# Format each datapoint
#formatted_list = [f'"{value}", "0"' for value in combined_datapoints]
formatted_list = [(float(value), 0) for value in combined_datapoints]

# Print the resulting formatted list
print(formatted_list)
len(formatted_list)



# Format each datapoint
formatted_list = [value for datapoint in formatted_list for value in (f"{datapoint}", '0')]

# Print the resulting formatted list
print(formatted_list)

formatted_list = [0.1, 0.9]  # TODO fix this 


traj = simulateForwardTrajectory(N=model8.init_size, fitness=None, 
    beginGen=0, 
    endGen=model8.num_gens, beginFreq=formatted_list,
#TODO
#add Emma's contemporary frequencies here, copy the same one twice since there is only
#one population with constant gene flow    
    endFreq=[[0, 1], [0, 1]])

#FIXME
#I keep getting the error:
#  Initial frequency should be provided for each locus (nLoci) or each locus at 
# each subpopulation (nLoci * len(N)).
#I've asked BO Peng and I'm waiting for his response.
#In the meantime, I could do just median values:

haplotype_array = np.array(haplotype_frequencies) #convert for easier manipulation
# Calculate median for each row
medians = np.median(haplotype_array, axis=1)
#medians = [0.00588362, 0.00716048]
#So, beginFreq will be: 
    
beginFreq = [[0.00588362, 0], [0.00716048, 0]]

# Set up logging configuration
logging.basicConfig(level=logging.INFO)

# Create a logger for traj function
logger = logging.getLogger("__simuPOP_simulateForwardTrajectory__")

#TODO select allele frequencies
traj = simulateForwardTrajectory(N=[1500, 1500], fitness=None, 
    #FIXME: select the generation where the bottleneck starts, so I can add the 
    #initial size to the trajectory (since it cannot handle > 6,000 individuals)
    beginGen=0, 
    endGen= 82, #gen, 
    #beginFreq=[0.2, 0.3],
    beginFreq = [[0.00588362, 0], [0.00716048, 0]],
    endFreq=[[0.1, 0], [0.2, 0]],
    logger=logger)

#ValueError: Invalid frequency range 0.100000 - 0.000000
#I think it will not let me simulate trajectory using an allele that is always
#zero. Wait until Bo replies, but maybe I will need to use SNPs instead.



# =============================================================================
# Down here: the code merged with forwardTrajectory.py that worked
# =============================================================================

pop.evolve(
    initOps = [sim.InitSex(), sim.IdTagger(), sim.PyOperator(func=init_native_breeding_grounds)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] + 
                   #this means that the mtDNA is the very last locus
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
    # increase age by 1
    preOps = [sim.InfoExec('age += 1')],
    matingScheme=sim.ControlledRandomMating(
        ops=[sim.Recombinator(rates=0.01)],
        loci=5, alleles=1, freqFunc=traj.func()),
    postOps = [
        sim.Stat(structure=range(1), subPops=sub_population_names, suffix='_AB', step=10),
        sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10)
    ],
    gen = 80
)



traj = simulateForwardTrajectory(N=[2000, 4000], fitness=None,
    beginGen=0, endGen=100, beginFreq=[0.2, 0.3],
    endFreq=[[0.1, 0.11], [0.2, 0.21]])
# 
#traj.plot('log/forwardTrajectory.png', set_ylim_top=0.5,
#    plot_c_sp=['r', 'b'], set_title_label='Simulated Trajectory (forward-time)')
pop = sim.Population(size=[2000, 4000], loci=10, infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.8, 0.2], subPops=0),
        sim.InitGenotype(freq=[0.7, 0.3], subPops=1),
        sim.PyOutput('Sp0: loc2\tloc5\tSp1: loc2\tloc5\n'),
    ],
    matingScheme=sim.ControlledRandomMating(
        ops=[sim.Recombinator(rates=0.01)],
        loci=5, alleles=1, freqFunc=traj.func()),
    postOps=[
        sim.Stat(alleleFreq=[2, 5], vars=['alleleFreq_sp'], step=20),
        sim.PyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\n' % (subPop[0]['alleleFreq'][2][1],"
            "subPop[0]['alleleFreq'][5][1], subPop[1]['alleleFreq'][2][1],"
            "subPop[1]['alleleFreq'][5][1])", step=20),
        
    ],
    gen = 101
)


# =============================================================================
# Code that replaced randomMating for ControlledrandomMating, now obsolete
# =============================================================================
pop.evolve(
    initOps = [sim.InitSex(), sim.IdTagger(), sim.PyOperator(func=init_native_breeding_grounds)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], 
                                     freq=haplotype_frequencies[i], 
                                     loci=[nb_loci]) for i in range(0, sub_population_count)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], 
                                     freq=[snp[n][i], 1-snp[n][i]], 
                                     loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
    # increase age by 1
    
    preOps = [sim.InfoExec('age += 1')],
    #matingScheme = sim.HeteroMating([
    matingScheme = sim.ControlledRandomMating(
        #MK: removing options until it stops breaking:
            
        # age <= maxAge, copy to the next generation (weight=-1)
        # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
        # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
        # sim.CloneMating(ops=[sim.CloneGenoTransmitter(chroms=[0,1])],
          #                  subPops=[(sub_population, 6) for sub_population in range(0, sub_population_count)],
           #                 weight=-1),
        
        # Then we simulate random mating only in VSP 1 (ie reproductively mature individuals) within subpopulation (breeding/winter grounds)
       # sim.RandomMating(ops=[sim.MitochondrialGenoTransmitter(),
        #                          sim.MendelianGenoTransmitter(),
         #                         sim.IdTagger(),
          #                        sim.InheritTagger(mode=sim.MATERNAL, infoFields=['feeding_ground']),
           #                       sim.InheritTagger(mode=sim.MATERNAL, infoFields=['native_breeding_ground']),
            #                      sim.PedigreeTagger()],
             #                subPops=[(sub_population, 7) for sub_population in range(0, sub_population_count)],
              #               weight=1)],
            #subPopSize=configure_new_population_size,
               #                               ),

                              subPopSize=configure_new_population_size,
                              freqFunc=traj.func()),
    
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
        sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10)
    ],
    gen = 80
)

sim.dump(pop)

#sim.dump(pop, width=3, loci=[], subPops=[(sim.ALL_AVAIL, sim.ALL_AVAIL)], max=1000, structure=False);
#return



# =============================================================================
# Does it work with the heteromating, trajectory, bottleneck, a.k.a., 
# merging them all together?
# =============================================================================
# Emma's initial code, once again. It works!
pop.evolve(
    initOps = [sim.InitSex(), sim.IdTagger(), sim.PyOperator(func=init_native_breeding_grounds)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
    # increase age by 1
    preOps = [sim.InfoExec('age += 1')],
    matingScheme = sim.HeteroMating([
        # age <= maxAge, copy to the next generation (weight=-1)
        # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
        # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
        sim.CloneMating(ops=[sim.CloneGenoTransmitter(chroms=[0,1])],
                            subPops=[(sub_population, 6) for sub_population in range(0, sub_population_count)], ##6 because it's dividing the pop into the VSPs based on age and removing
                            #the dead ones
                            weight=-1), 
        # Then we simulate random mating only in VSP 1 (ie reproductively mature individuals) within subpopulation (breeding/winter grounds)
        sim.RandomMating(subPopSize=model,
            ops=[sim.MitochondrialGenoTransmitter(),
                                  sim.MendelianGenoTransmitter(),
                                  sim.IdTagger(),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['feeding_ground']),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['native_breeding_ground']),
                                  sim.PedigreeTagger()],
                             subPops=[(sub_population, 7) for sub_population in range(0, sub_population_count)],
                             weight=1),
        sim.ControlledRandomMating(freqFunc=traj.func(),
                                   weight=1)] 
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
        sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10), #Print Fst every 10 steps
        sim.Stat(alleleFreq=[2, 5], vars=['alleleFreq_sp'], step=10), #added this now, to
        #calculate the allele frequencies in selected loci
        sim.PyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\n' % (subPop[0]['alleleFreq'][2][1],"
            "subPop[0]['alleleFreq'][5][1], subPop[1]['alleleFreq'][2][1],"
            "subPop[1]['alleleFreq'][5][1])", step=10) #added this now, to
            #print out the allele frequencies in selected loci
       
    ],
    finalOps= sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
    gen = model4.num_gens
    #80
)


# =============================================================================
# Printing only mitochondrial frequencies:
# =============================================================================

pop.evolve(
    initOps = [sim.InitSex(), sim.IdTagger(), sim.PyOperator(func=init_native_breeding_grounds)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
    # increase age by 1
    preOps = [sim.InfoExec('age += 1')],
    matingScheme = sim.HeteroMating([
        # age <= maxAge, copy to the next generation (weight=-1)
        # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
        # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
        sim.CloneMating(ops=[sim.CloneGenoTransmitter(chroms=[0,1])],
                            subPops=[(sub_population, 6) for sub_population in range(0, sub_population_count)],
                            weight=-1), #if weights are negative, they are multiplied to their parental subpopulation;
            #EX: parental pop: (500, 1000), weight -2, next generation: (1000, 2000). For weight -1, it keeps the number of individuals
            #from the parental generation.
            #ALSO: if there is a mix of negative and positive weights, the negative will be processed first.
        # Then we simulate random mating only in VSP 1 (ie reproductively mature individuals) within subpopulation (breeding/winter grounds)
        sim.RandomMating(subPopSize=model4,
            ops=[sim.MitochondrialGenoTransmitter(),
                                  sim.MendelianGenoTransmitter(),
                                  sim.IdTagger(),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['feeding_ground']),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['native_breeding_ground']),
                                  sim.PedigreeTagger()],
                             subPops=[(sub_population, 7) for sub_population in range(0, sub_population_count)],
                             weight=1),
        sim.ControlledRandomMating(freqFunc=traj.func(),
                                   weight=1)] #we decided to keep the same weight as the
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
                 haploFreq=[], 
                 begin = 70, step = 75, end = 85),     
        sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10), #Print Fst every 10 steps
        sim.Stat(alleleFreq=[1, 2, 3, 4, 5, 6, 7, 8, 9, 100], vars=['alleleFreq_sp'], step=10), #added this now, to
        #calculate the allele frequencies in selected loci
       sim.PyEval(r"'%.2f\t%.2f\n' % ("
           "subPop[0]['alleleFreq'][100][1]," 
           "subPop[1]['alleleFreq'][100][1])", step=10),
      # sim.PyEval(r"'%.2f' % haploFreq[0]", step=10)

        #to print out the allele frequencies in selected loci
       
    ],
    finalOps= sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
    gen = model4.num_gens
    )


###############################################################################
# =============================================================================
# WHERE TO GO NEXT:
# =============================================================================

    # =============================================================================
    # 1. FIXED! How do I store and print allele frequencies for trajectory?
    # =============================================================================

#Edit the PyEval to print out the allele frequencies in each generation, so I can
#export and plot the values. From the forwardTrajectory.py script:
    

pop.dvars().traj
    # Forward-time trajectory plot https://simupop.sourceforge.net/manual_release/build/userGuide_ch7_sec2.html


traj.plot('log/forwardTrajectory.png', set_ylim_top=0.5,
          plot_c_sp=['r', 'b'], set_title_label='Simulated Trajectory (forward-time)')    

    #Bo Peng himself does not recommend using the plot option within simuPOP,
    #So I will work on projecting the trajectory with the output, as suggested by him
    #Using R (or even python, by loading my dataset into it.)
    #He mentions it in:  https://github.com/BoPeng/simuPOP/issues/73
                        #https://github.com/BoPeng/simuPOP/issues/56

    # =============================================================================
    # 2. FIXED! How do I add a bottleneck? Below, the original InstantChangeModel script,
    # which I edited for my population simulation above.
    # =============================================================================
    
    #Maybe useful link?
    #https://simupop.sourceforge.net/manual_release/build/userGuide_ch7_sec3.html?highlight=stats

#Does this work for decreasing the size of the population?
import simuPOP as sim
from simuPOP.demography import *
model = MultiStageModel([
    InstantChangeModel(T=200, 
        # start with an ancestral population of size 1000
        N0=(2000, 'Ancestral'),
        # change population size at 50 and 60
        G=[40, 70], 
        # change to population size 200 and back to 1000
        NG=[(200, 'bottleneck'), (1000, 'Post-Bottleneck')]),
    ExponentialGrowthModel(
        T=80, 
        # split the population into two subpopulations
        N0=[(400, 'P1'), (600, 'P2')],
        # expand to size 4000 and 5000 respectively
        NT=[4000, 5000])]
    )

# model.init_size returns the initial population size
# migrate_to is required for migration
pop = sim.Population(size=model.init_size, loci=1,
    infoFields=model.info_fields)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5])
    ],
    matingScheme=sim.RandomMating(subPopSize=model),
    finalOps=
        sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
    gen=model.num_gens
)


# get a visual presentation of the demographic model
model.plot('log/demoModel.png',
    title='A bottleneck + exponential growth demographic model')

    # =============================================================================
    # 3. How do I store haplotype frequencies so I can print them?
    # A.K.A.: Output statistics
    # =============================================================================


    
    #Useful link?
    #https://simupop.sourceforge.net/manual_release/build/userGuide_ch8_sec4.html?highlight=stats

#To work this out, I will modify Emma's original code that prints out allele
#frequencies:
#haploFreq needs a list of haplotype positions, ALL_AVAIL (used for alleleFreq
#does not work)
sim.stat(pop, haploFreq = [[0,1], [2,3]])
#Does it make sense to have alleles for mtDNA?

viewVars(pop.vars()) #this then prints out the haplotype frequencies for the
#alleles at the specified positions (0,1, 2, 3)

#From https://simupop.sourceforge.net/manual_release/build/refManual_ch3_sec11.html?highlight=rep:
    # I can print out the stat for particular generations, which is my goal.
    # I want stats for:
      
        # Gen 70: 1830-1855 maximum amount of whales
        # Gen 75: 1930-1955 lowest amount (bottleneck)
        # Gen 78: 1980-2005 recovering 
        # Gen 85: 2010-2030 current(ish) generation
        
sim.stat(pop,
         numOfMales=True, 
         haploFreq=[], 
         begin = 70, step = 75, end = 85)        

viewVars(pop.vars())
        
       
        
       
        
       
# =============================================================================
# 21/11 meeting to do list:        
# =============================================================================
       
        
       
#TODO
#Count the number of haplotypes in the last loci, that is where the
#mitochondrial frequencies are.

#TODO
#Output the sim.dump(pop) at every generation, to make sure that there are no
#random individuals are created out of the blue.

#TODO
#Create a report with all the methods I used and choices I made for each point of
#the script, plus what's left to do.

        
        
        
        
        
        
        
        
        
        
        
        
        
        


haplofrequencies = pop.dvars().haploFreq
with open('haplofreq.txt', 'w') as haplofreqfile:
    index = 0
    for locus in haplofrequencies:
        if (locus == nb_loci):
            continue
        if (len(haplofrequencies[locus]) < 2):
            continue
        print(index, end=' ', file=haplofreqfile)
        index = index + 1
        for allele in haplofrequencies[locus]:
            print(haplofrequencies[locus][allele], end=' ', file=haplofreqfile)
        print(file=haplofreqfile)
    
































    # =============================================================================
    # 4. What weight to add to the mating scheme?
    # =============================================================================
    #Useful link: https://github.com/BoPeng/simuPOP/issues/75


# =============================================================================
# Export population - add # once done to not overwrite files by accident
# =============================================================================
#sim.utils.export(pop, format="csv", output="popNov16th.csv", infoFields=['ind_id', 'father_id', 
#                                                'mother_id', 'nitrogen', 'carbon', 
#                                                'feeding_ground', 'native_breeding_ground', 
#                                                'migrate_to'], sexFormatter ={1:'M', 2:'F'})






# =============================================================================
# Final comments on Emma's original code
# =============================================================================

# Plan
# * Add simulation of mitochondrial DNA
# * Break populations into VSPs and assign breeding scheme amongst VSPs
# * SNP


# Subpopulation 1: individuals who breed at wintering ground 1
# Subpopulation 2: individuals who breed at wintering ground 2

# VSP 1: individuals who feed at feeding ground 1
# VSP 2: individuals who feed at feeding ground 2



#   Migration ***
# Only males migrate
# Migrants spend 1 year in the opposing breeding grounds then return to their home population
# On average 1% of the males migrate to the other breeding ground in a given season
