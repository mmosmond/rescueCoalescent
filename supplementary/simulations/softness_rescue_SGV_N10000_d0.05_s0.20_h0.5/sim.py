import os

N0 = 1e4 #intial population size and carrying capacity
d = 0.05 #decline rate
s = 0.2 #selection coefficient
h = 0.5 #dominance
B = 2 #number of offspring per parent
u = 0 #mutation rate at selected site
m = 0 #migration rate
datadir = "data/" #location to put output files
nreps = 100 #number of replicates per k value
ks = [i for i in range(2,102+1,10)] #kvalues to explore

for j,k in enumerate(ks):
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

			initializeSLiMModelType("nonWF"); //non Wright Fisher model
			initializeMutationType("m1", h, "f", s); //beneficial mutation characteristics
			m1.mutationStackPolicy = "f"; //keep first mutation
			initializeGenomicElementType("g1", m1, 1.0); //define element g1 to have beneficial mutations
			initializeGenomicElement(g1, 0, 0); //element g1 is just one site
			initializeMutationRate(u, 0); //mutation rate per site
			initializeRecombinationRate(0);   
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
				i.addNewDrawnMutation(m1, 0); //... and give beneficial mutation
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
		 			target.addNewDrawnMutation(m1, 0); //add the migrant allele	
				}
			}
			
			// census offspring
			N = length(off); //population size
			freq = sum(asInteger(off.genomes.countOfMutationsOfType(m1)>0))/(2*N); //frequency of beneficial mutation
			if ((u==0 & m==0 & freq == 0) | (N==0)) { //if fail to adapt
				catn("all hope was lost in generation " + (sim.generation-1) + " - RESTARTING"); //alert the user
				sim.readFromPopulationFile("/tmp/slim_" + simID + ".txt"); //reinitialize simulation
				setSeed(getSeed() + 1); //change random seed
			}
			else {           
		        if (freq == 1.0 & N == N0) { //if mutation fixed and population recovered
					catn("rescue complete in generation " + (sim.generation-1.0));
		            n = length(unique(off.genomes.mutations));
		            catn(n + " mutations contributed");
		            writeFile("data/nalleles.csv", paste(c(k, n), ","), append=T);
		            sim.simulationFinished(); //end simulation
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
		""" %(N0,d,s,h,B,u,m,k)

		os.system("echo '" + script + "' | slim") #run script in SLiM

