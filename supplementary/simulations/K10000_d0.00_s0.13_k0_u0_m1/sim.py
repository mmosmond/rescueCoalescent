import os #to run SLiM from python

# parameters for SLiM
K = 1e4 #carrying capacity
N0 = K #intial population size
d = 0.0 #decline rate of ancestral homozygote
s = 0.13 #selection coefficient
h = 0.5 #dominance coefficient
B = 2 #number of offspring per parent
u = 0 #mutation rate at selected site
m = 1 #migration rate
k = 0 #initial number of beneficial alleles
datadir = "data/" #location to put output files
nreps = 100 #number of replicates (number of runs where allele fixes and population recovers)
maxt = 1000 #maximum number of generations (safestop that should never be reached)

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
		defineConstant("u", %.8f); //mutation rate at beneficial locus [0,1]
		defineConstant("m", %.8f); //migration rate [0,1]
		defineConstant("k", %d); //initial number of mutants
		defineConstant("outfile_dynamics", "%sdynamics_%d.txt"); //where to save dynamics
		defineConstant("simID", getSeed()); //get the random seed to label temporary file

		initializeSLiMModelType("nonWF"); //non Wright Fisher model
		initializeMutationType("m1", h, "f", s); //beneficial mutation characteristics
		m1.mutationStackPolicy = "l"; //keep last mutation
		initializeGenomicElementType("g1", m1, 1.0); //define element g1 to have beneficial mutations
		initializeGenomicElement(g1, 0, 0); //element g1 is just one site
		initializeMutationRate(u, 0); //mutation rate per site 
		initializeRecombinationRate(0); //recombination rate between sites
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
				i.addNewDrawnMutation(m1, 0); //... and give beneficial mutation
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
	 			target.addNewDrawnMutation(m1, 0); //add the migrant allele	
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
				catn("rescue complete in generation " + sim.generation); //alert the user
	            sim.simulationFinished(); //end simulation
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
	""" %(K,N0,d,s,h,B,u,m,k,datadir,i,maxt,maxt+1)

	# and run it
	os.system("echo '" + script + "' | slim") 
