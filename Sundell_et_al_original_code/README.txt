# A small documentation of files in this arcive, JI. Kammonen 23/11/2023:
simuBottlenecks_v2_serial.py	-This is the main simuPOP script for a "Standard" run.
				 You will need a computer with ca. 20GB RAM to be able to finish the run.

conf_file.txt 			-A hard coded list of values that the main script reads. The contents of the lines has not been documented.
				 Frankly, I cannot remember what each line does but there's e.g. sample sizes, mutation and migration rates and
				 variable mathematical functions to model the population growth/decline as per Fig. 1 of Sundell et al. 2013 (attached).
				 This specific file corresponds to simulation scenario "E1" in Sundell et al. 2013 (attached). We concluded that
				 E1 is one of the 4 scenarios where the simulated genetic diversity indices all arrived closest to those
				 observed in actual population level data from Finland.

crs_based_sequences.txt		-A list of 10 "founding" mtDNA haplotypes from the HVS of the H.sapiens mtDNA.
				 The sequences are in simuPOP-readable format where the 4 nucleotide bases are represented by numbers 0,1,2,3.
				 There's no file for the Y-STR haplotypes since these were much simpler and could be defined the main script.

Sundell_et_al_2013.pdf		-For reference, our more recent publication where we used the main script in this archive. Attaching it here
				 as it is a bit challenging to find it online. The publication volume is the "CAA2012 : Proceedings of the 40th
				 Conference in Computer Applications and Quantitative Methods in Archaeology". I apologize for the poor printing
				 quality of some of the figures.
