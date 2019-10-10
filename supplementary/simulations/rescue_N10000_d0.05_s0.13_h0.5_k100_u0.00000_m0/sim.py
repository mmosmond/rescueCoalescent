import pyslim, tskit, msprime, os
import matplotlib.pyplot as plt
import numpy as np
import csv
import time

N0 = 1e4 #intial population size and carrying capacity
Ne = round(N0*4/7) #long-term effective population size (for recapitation)
d = 0.05 #decline rate
s = 0.13 #selection coefficient
h = 0.5 #dominance
B = 2 #number of offspring per parent
L = 2e7 #length of chromosome
L0 = round(L/2) #location of selected site
u = 0 #mutation rate at selected site
rbp = 2e-8 #per basepair recombination rate
U = 6e-9 #per basepair mutation rate at neutral sites
k = 100 #initial number of beneficial alleles
m = 0 #migration rate
datadir = "data/" #location to put output files
nreps = 100 #number of replicates
nwindows = 100 #number of windows across genome to compute stats for
nsamples = 100 #number of chromosomes to sample

X = np.zeros(nreps) #empty vector for number of mutations in each sim

for i in range(nreps):

	script = """
	initialize() {
		
		defineConstant("N0", %d); //initial pop size (integer)
		defineConstant("d", %f); // wildtype decline rate [0,1]
		defineConstant("s", %f); //beneficial selection coefficient ([0,1]; s>d for rescue to be possible) 
		defineConstant("h", %f); //beneficial dominance [0,1]
		defineConstant("B", %d); //offspring per parent (positive integer; must be great than 1 for possible persistence) 
		defineConstant("L", %d - 1); //number of sites (positive integer)
		defineConstant("L0", %d - 1); //site number of beneficial locus (positive integer, L0<L)
		defineConstant("u", %.8f); //mutation rate at beneficial locus [0,1]
		defineConstant("m", %.8f); //migration rate [0,1]
		defineConstant("rbp", %.8f); //recombination rate per base pair [0,1]
		defineConstant("k", %d); //initial number of mutants
		defineConstant("outfile_dynamics", "%sdynamics_%d.txt"); //where to save dynamics
		defineConstant("outfile_tree", "%stree_%d.trees"); //where to save tree

		initializeSLiMModelType("nonWF"); //non Wright Fisher model
		initializeTreeSeq(); //record the tree
		initializeMutationType("m1", h, "f", s); //beneficial mutation characteristics
		m1.mutationStackPolicy = "f"; //keep first mutation
		initializeMutationType("m2", 0.5, "f", 0.0); //neutral mutations (heritability has no affect)
		initializeGenomicElementType("g1", m1, 1.0); //define element g1 to have beneficial mutations
		initializeGenomicElementType("g2", m2, 1.0); //define element g2 to have neutral mutations
		initializeGenomicElement(g1, L0, L0); //element g1 is just one site
		initializeGenomicElement(g2, 0, L0 - 1); // element g2 is everything to the left...
		initializeGenomicElement(g2, L0 + 1, L); // ...and everything to the right of LO
		initializeMutationRate(c(0,u,0), c(L0-1, L0, L)); //mutation rate per site 
		initializeRecombinationRate(rbp); //recombination rate between sites
		writeFile(outfile_dynamics, "t n p"); //start writing to the dynamics file     
	}

	reproduction() { //occurs immediately before early events
		for (i in 1:B) //B matings per parent
			subpop.addCrossed(individual, subpop.sampleIndividuals(1)); //random mating, 1 offspring per pair
	}

	//initialize population
	1 early() {
		sim.addSubpop("p1", N0); //initialize population of wildtypes
		target = sample(p1.genomes, k); //choose k chromosomes...
		for (i in target)
			i.addNewDrawnMutation(m1, L0); //... and give beneficial mutation
		defineConstant("simID", getSeed()); //get the random seed to make sure any future runs use a different seed	
		sim.outputFull("/tmp/slim_" + simID + ".txt"); //output this initial state to use for future runs if needed
	}

	//discrete generations, hard carrying capacity, census and update fitness
	early() {
		//enforce discrete generations
		inds = sim.subpopulations.individuals; //get info on all individuals
		inds[inds.age > 0].fitnessScaling = 0.0; //parents all die (discrete generations)

		//hard carrying capacity by random culling
		off = inds[inds.age == 0]; //offspring
		N = length(off); //total number of offspring
		indices = which(inds.age == 0); //indices of offspring
		if (N > N0) { //if more than N0...
			inds[sample(indices, N-N0)].fitnessScaling = 0.0; //...kill a random subset to reduce N to N0
			off = inds[inds.fitnessScaling > 0]; //get surviving offspring
		}	

		// migration
		if (m>0 & N>0) { //if adapting from migration and some offspring made
			if (runif(1)<m) { //with probability m
				target = sample(off.genomes, 1); //choose a chromosome to add a migrant allele to
	 			target.addNewDrawnMutation(m1, L0); //add the migrant allele	
			}
		}
		
		// census offspring
		N = length(off); //population size
		freq = sum(asInteger(off.genomes.countOfMutationsOfType(m1)>0))/(2*N); //frequency of beneficial mutation
		if ((u==0 & m==0 & freq == 0) | (N==0)) { //if fail to adapt
	        catn(sim.generation-1 + ": " + N); //print generation and population size
	        writeFile(outfile_dynamics, "t n p"); //erase and restart the output file
			catn("all hope was lost in generation " + (sim.generation-1) + " - RESTARTING"); //alert the user
			sim.readFromPopulationFile("/tmp/slim_" + simID + ".txt"); //reinitialize simulation
			setSeed(getSeed() + 1); //change random seed
		}
		else {           
			catn(sim.generation-1 + ": " + N + ", " + freq); //print generation and population size and frequency
			writeFile(outfile_dynamics, sim.generation-1 + " " + N + " " + freq, append=T);
	        if (freq == 1.0 & N == N0) { //if mutation fixed and population recovered
				catn("rescue complete in generation " + (sim.generation-1.0));
	            sim.simulationFinished(); //end simulation
	            sim.treeSeqOutput(outfile_tree); //save tree sequence
	        }
		}

		//fitness scaling
		p1.fitnessScaling = (1.0 - d)/B; //scale fitness so that wildtype has multiplicative growth rate B*(1-d)/B = (1-d) and the mutant B*(1+s)(1-d)/B = (1+s)(1-d) ~ 1+(s-d)
	}

	//backup: end simulation if runs too long
	1000 late () {
		catn("times up, make max gens longer");
		sim.simulationFinished();
	}
	""" %(N0,d,s,h,B,L,L0,u,m,rbp,k,datadir,i,datadir,i)

	# start = time.time()

	os.system("echo '" + script + "' | slim") #run script in SLiM

	ts = pyslim.load("%stree_%d.trees" %(datadir, i)) #load tree-sequences
	X[i] = len(np.unique([i.genotypes for i in ts.variants()][0])) #number of ubnique copies of beneficial allele
	ts = msprime.mutate(ts.recapitate(rbp, Ne=Ne).simplify(), rate=U, keep=True) #recapitate and overlay neutral mutations
	offspring_nodes = (np.where(ts.tables.nodes.time==0.)[0]).astype(np.int32) #chromosomes in offspring (to exclude the parental generation from the samples)
	samples = np.random.choice(offspring_nodes, nsamples, replace=False) #random sample of chromosomes in offspring
	ts = ts.simplify(samples) #simplify to sample only 
		
	windows = np.linspace(0, ts.sequence_length, nwindows+1) #delimit windows
	site_div = ts.diversity([ts.samples()], windows=windows) #compute average pairwise diversity in windows
	tajimasD = ts.Tajimas_D([ts.samples()], windows=windows) #compute average Tajima's D in windows

	site_div = [i[0] for i in site_div] #flatten site_div list
	tajimasD = [i[0] for i in tajimasD] #flatten tajimasD list
	midpoints = windows[:-1]+(windows[1:]-windows[:-1])/2 #midpoint of each window
	distances = midpoints-L0 #distance from selected site
	recombination = [(1-(1-2*rbp)**abs(i))/2*np.sign(i) for i in distances] #recombination rate (signed for plotting)

	# end = time.time()
	# print(end-start)

	# plt.plot(recombination, meanH) #plot heterozygosity
	# plt.axvline(0, linestyle=":") #add location of sweep
	# plt.show ()

	# save to file
	csvData = zip(recombination, site_div, tajimasD)
	with open('%sstats_%d.csv' %(datadir,i), 'w') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerows(csvData)
	csvFile.close()

np.savetxt('%snalleles.txt' %datadir, X, fmt='%d')
