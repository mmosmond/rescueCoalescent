import os #to run SLiM from python
import numpy as np #for vectors etc
import pyslim, tskit, msprime #to recapitate, mutate, and sample tree sequences, and compute statistics
import csv #to save statistics

# parameters for SLiM
K = 799 #carrying capacity
N0 = K #intial population size
d = 0.00 #decline rate of ancestral homozygote
s = 0.13 #selection coefficient
h = 0.5 #dominance coefficient
B = 2 #number of offspring per parent
L = 2e7 #number of sites on chromosome
L0 = round(L/2) #location of selected site (one of the center sites)
u = 0 #mutation rate at selected site
m = 0 #migration rate
rbp = 2e-8 #per basepair recombination rate
k = 10 #initial number of beneficial alleles
datadir = "data/" #location to put output files
nreps = 100 #number of replicates (number of runs where allele fixes and population recovers)
maxt = 1000 #maximum number of generations (safestop that should never be reached)

#parameters for msprime
Ne = round(1e4 * 4/7) #long-term effective population size
U = 6e-9 #per basepair mutation rate at neutral sites
nsamples = 100 #number of chromosomes to sample for stats
nwindows = 100 #number of windows across genome to compute stats for
R = 0.001 #recombination distance between neutral loci to calculate LD for

#genome window calculation
windows = np.linspace(0, L, nwindows+1) #delimit windows
midpoints = windows[:-1]+(windows[1:]-windows[:-1])/2 #midpoint of each window
distances = midpoints-L0 #distance of midpoint from selected site
recombination = [(1-(1-2*rbp)**abs(i))/2*np.sign(i) for i in distances] #recombination rate at midpoint (signed for plotting)

#empty vector for number of unique copies of beneficial allele remaining at fixation in each replicate
X = np.zeros(nreps)

#for each replicate
for i in range(nreps): 
	
	# define SLiM script ...
	script = """
	initialize() {
		
		defineConstant("K", %d); //carrying capacity (integer)
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
		defineConstant("simID", getSeed()); //get the random seed to label temporary file

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

	//discrete generations, hard carrying capacity, census and update fitness
	1:%d early() {

		//initialize population
		if (sim.generation == 1) {
			sim.addSubpop("p1", N0); //initialize population of wildtypes
			target = sample(p1.genomes, k); //choose k chromosomes without replacement...
			for (i in target)
				i.addNewDrawnMutation(m1, L0); //... and give beneficial mutation
			sim.outputFull("/tmp/slim_" + simID + ".txt"); //output this initial state to use for future runs if needed
		}

		//enforce discrete generations
		inds = sim.subpopulations.individuals; //get info on all individuals
		inds[inds.age > 0].fitnessScaling = 0.0; //parents all die at next instance of viability selection

		//hard carrying capacity by random culling
		off = inds[inds.age == 0]; //offspring
		N = length(off); //total number of offspring
		indices = which(inds.age == 0); //indices of offspring
		if (N > K) { //if more than K...
			inds[sample(indices, N-K)].fitnessScaling = 0.0; //...kill a random subset to reduce N to K
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
	        writeFile("data/prescue.csv", "0", append=T); //record extinction
	        writeFile(outfile_dynamics, "t n p"); //erase and restart the output file
			catn("all hope was lost in generation " + sim.generation + " - RESTARTING"); //alert the user
			sim.readFromPopulationFile("/tmp/slim_" + simID + ".txt"); //reinitialize simulation
		}
		else {           
			catn(sim.generation + ": " + N + ", " + freq); //print generation and population size and frequency
			writeFile(outfile_dynamics, sim.generation + " " + N + " " + freq, append=T);
	        if (freq == 1.0 & N == K) { //if mutation fixed and population recovered
	        	writeFile("data/prescue.csv", "1", append=T); //record rescue
				catn("rescue complete in generation " + sim.generation);
	            sim.simulationFinished(); //end simulation
	            sim.treeSeqOutput(outfile_tree); //save tree sequence
	        }
		}

		//fitness scaling (viability selection occurs after early events)
		p1.fitnessScaling = (1.0 - d)/B; //survival probability V = X(1-d)/B, where X is the fitness effect of the selected site (X=1 for wildtype, X=1+s*h for heterozygotes, X=1+s for mutant homozygotes)
	}

	//backup: end simulation if runs too long and print warning to increase maxt
	%d late () {
		catn("times up, make maxt longer!");
		sim.simulationFinished();
	}
	""" %(K,N0,d,s,h,B,L,L0,u,m,rbp,k,datadir,i,datadir,i,maxt,maxt+1)

	# and run it
	os.system("echo '" + script + "' | slim") 

	# then load tree-sequences
	ts = pyslim.load("%stree_%d.trees" %(datadir, i)) 
	
	# calculate the number of unique copies of beneficial allele
	X[i] = len(np.unique([i.genotypes for i in ts.variants()][0])) 
	
	# recapitate the tree sequence and overlay neutral mutations
	ts = msprime.mutate(ts.recapitate(rbp, Ne=Ne).simplify(), rate=U, keep=True)
	
	# take a random sample
	offspring_nodes = (np.where(ts.tables.nodes.time==0.)[0]).astype(np.int32) #chromosomes in offspring (to exclude the parental generation from the samples)
	samples = np.random.choice(offspring_nodes, nsamples, replace=False) #random sample of chromosomes in offspring
	ts = ts.simplify(samples) #simplify to sample only
	ts.dump("%stree_%d_sample.trees" %(datadir, i)) #save sample tree (in case need to calculate any more statistics)
	# ts = tskit.load("%stree_%d_sample.trees" %(datadir,i)) #load dumped version if running new stats below

	# calculate pairwise diversity
	site_div = ts.diversity([ts.samples()], windows=windows) #compute average pairwise diversity in windows
	site_div = [i[0] for i in site_div] #flatten site_div list
	
	# calculate Tajima's D
	tajimasD = ts.Tajimas_D([ts.samples()], windows=windows) #compute average Tajima's D in windows
	tajimasD = [i[0] for i in tajimasD] #flatten tajimasD list
	
	# save pairwise diversity and Tajima's D to file, with associated recombination rates
	csvData = zip(recombination, site_div, tajimasD)
	with open('%sstats_%d.csv' %(datadir,i), 'w') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerows(csvData)
	csvFile.close()

	# calculate site frequency spectrum
	sfs = ts.allele_frequency_spectrum([ts.samples()], windows=windows, polarised=True) #unfolded SFS, averaged within windows
	np.savetxt('%ssfs_%d.txt' %(datadir,i), sfs)

	# calculate linkage disequilibrium
	all_positions = [i.position for i in ts.mutations()] #positions of all mutations
	freqs = [sum(i.genotypes)/nsamples for i in ts.variants()] #frequency of derived alleles at these positions
	seg = [i for i,j in enumerate(freqs) if 0<j and j<1] #indices of segregating mutations
	positions = [all_positions[i] for i in seg] #positions of segregating mutations
	idxs = [np.argmin(np.abs(positions-i)) for i in midpoints] #find mutation indices nearest the window midpoints
	idx_positions = [positions[i] for i in idxs] #positions of mutations nearest midpoints
	distance = np.log(1-R)/np.log(1-2*rbp) #convert the specified recombination rate between mutations to number of sites
	other_positions = idx_positions + distance #this is where we want the other mutation
	other_idxs = [np.argmin(np.abs(positions-i)) for i in other_positions] #mutation indices nearest the desired poition
	lds = [tskit.LdCalculator(ts).r2(seg[idxs[i]], seg[other_idxs[i]]) for i in range(nwindows)] #linkage disequilibrium between the two mutations

	# save ld to file with associated recombination rates
	csvData = zip(recombination, lds)
	with open('%sld_%d.csv' %(datadir,i), 'w') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerows(csvData)
	csvFile.close()

# save number of unique copies of beneficial allele remaining in each replicate
np.savetxt('%snalleles.txt' %datadir, X, fmt='%d')
