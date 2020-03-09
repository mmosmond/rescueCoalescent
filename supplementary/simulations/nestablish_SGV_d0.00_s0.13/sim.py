import os #to run SLiM from python

N0 = 1e4 #intial population size and carrying capacity
d = 0.0 #decline rate
s = 0.13 #selection coefficient
h = 0.5 #dominance coefficient
B = 2 #number of offspring per parent
u = 0 #mutation rate at selected site
m = 0 #migration rate
datadir = "data" #location to put output files
ks = [1,10,20,40,60,80,100] #k values to explore
nreps = 100 #number of replicates per k value (number of runs where allele fixes and population recovers)
maxt = 1000 #maximum number of generations (safestop that should never be reached)

#for each k value
for k in ks:
	# for each replicate
	for i in range(nreps):

		script = """
		initialize() {
			
			defineConstant("N0", %d); //initial pop size (integer)
			defineConstant("d", %f); // wildtype decline rate [0,1]
			defineConstant("s", %f); //beneficial selection coefficient ([0,1]; s>d for rescue to be possible) 
			defineConstant("h", %f); //beneficial dominance [0,1]
			defineConstant("B", %d); //offspring per parent (positive integer; must be great than 1 for possible persistence) 
			defineConstant("u", %.8f); //mutation rate at beneficial locus [0,1]
			defineConstant("m", %.8f); //migration rate [0,1]
			defineConstant("k", %d); //initial number of mutants
			defineConstant("outfile_nestablish", "%s/nestablish_k%d.txt"); //where to save dynamics
			defineConstant("simID", getSeed()); //get the random seed to label temporary file

			initializeSLiMModelType("nonWF"); //non Wright Fisher model
			initializeMutationType("m1", h, "f", s); //beneficial mutation characteristics
			m1.mutationStackPolicy = "f"; //keep first mutation
			initializeGenomicElementType("g1", m1, 1.0); //define element g1 to have beneficial mutations
			initializeGenomicElement(g1, 0, 0); //element g1 is just one site
			initializeMutationRate(u, 0); //mutation rate per site
			initializeRecombinationRate(0); //per basepair recombination rate  
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
		 			target.addNewDrawnMutation(m1, 0); //add the migrant allele	
				}
			}
			
			// census offspring
			N = length(off); //population size
			freq = sum(asInteger(off.genomes.countOfMutationsOfType(m1)>0))/(2*N); //frequency of beneficial mutation
			n = length(unique(off.genomes.mutations)); //number of unique copies remaining
			if ((u==0 & m==0 & freq == 0) | (N==0)) { //if fail to adapt
				writeFile(outfile_nestablish, "0", append=T); //record number that establish
				catn("all hope was lost in generation " + sim.generation + " - RESTARTING"); //alert the user of restart
				sim.readFromPopulationFile("/tmp/slim_" + simID + ".txt"); //restart simulation
			}
			else {           
		        if (freq > 0.9 & N == N0) { //if mutation (essentially) fixed and population recovered
					writeFile(outfile_nestablish, asString(n), append=T); //record number that establish
					catn("rescue complete in generation " + sim.generation); //alert the uer of rescue
		            sim.simulationFinished(); //end simulation
		        }
			}

			//fitness scaling (viability selection occurs after early events)
			p1.fitnessScaling = (1.0 - d)/B; //survival probability V = X(1-d)/B, where X is the fitness effect of the selected site (X=1 for wildtype, X=1+s*h for heterozygotes, X=1+s for mutant homozygotes)
		}

		//backup: end simulation if runs too long
		%d late () {
			catn("times up, make maxt bigger");
			sim.simulationFinished();
		}
		""" %(N0,d,s,h,B,u,m,k,datadir,k,maxt,maxt+1)

		os.system("echo '" + script + "' | slim") #run script in SLiM

