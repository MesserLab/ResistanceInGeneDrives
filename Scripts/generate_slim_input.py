#!/usr/bin/env python
from argparse import ArgumentParser
from math import floor, sqrt

DEFAULT = {
        "N": 10000,
        "Sww": 0.0,
        "Swr01": 0.0,
        "Swr11": 0.0,
        "SwD": -0.05,
        "Sr01r11": 0.0,
        "Sr01r01": 0.0,
        "Sr01D": -0.05,
        "Sr11r11": 0.0,
        "Sr11D": -0.05,
        "SDD": -0.1,
        "b": 0.35,
        "d": 0.99,
        "c1":0.75,
        "c2":0.95,
        "X": False,
        "two_allele": False,
}

class params:
    def __init__(self, N=DEFAULT["N"], Sww=DEFAULT["Sww"], Swr01=DEFAULT["Swr01"], Swr11=DEFAULT["Swr11"], SwD=DEFAULT["SwD"], Sr01r11=DEFAULT["Sr01r11"],\
            Sr01r01=DEFAULT["Sr01r01"], Sr01D=DEFAULT["Sr01D"], Sr11r11=DEFAULT["Sr11r11"], Sr11D=DEFAULT["Sr11D"], SDD=DEFAULT["SDD"], d=DEFAULT["d"], b=DEFAULT["b"], c1=DEFAULT["c1"], c2=DEFAULT["c2"],\
            X=DEFAULT["X"], two_allele=DEFAULT["two_allele"]):
        self.N = N
        self.Sww = Sww
        self.Swr01=Swr01
        self.Swr11=Swr11
        self.SwD=SwD
        self.Sr01r01=Sr01r01
        self.Sr01r11=Sr01r11
        self.Sr01D=Sr01D
        self.Sr11r11=Sr11r11
        self.Sr11D=Sr11D
        self.SDD=SDD
        self.b=b
        self.d=d
        self.c1=c1
        self.c2=c2
        self.X = X
        self.two_allele = two_allele

def init(args):
    out = """
initialize() {
        initializeMutationRate( 0.0 );
        initializeMutationType("m1", 0.5, "f", -0.05); //D
        initializeMutationType("m2", 0.5, "f", 0.0); //+
        initializeMutationType("m3", 0.5, "f", 0.0); //r01
        initializeMutationType("m4", 0.5, "f", 0.0); //r10
        initializeMutationType("m5", 0.5, "f", 0.0); //r11
        initializeGenomicElementType("g1", m1, 1.0);
        initializeGenomicElement(g1, 0, 1);
        initializeRecombinationRate(0.0);\n"""
    if args.X:
        out += "    initializeSex('X');\n"
    else:
        out += "    initializeSex('A');\n"
    out += """}
1 {
	sim.addSubpop("p1", """ + str(args.N) +  """);
}

1 late() {\n"""

    if args.X:
        out += """
	xGenomes = p1.genomes[!p1.genomes.isNullGenome];
	xGenomes.addNewDrawnMutation(m2, 0); //add wildtype

        femaleGenomes = p1.individuals[p1.individuals.sex == "F"].genomes[0:99];
	femaleGenomes.removeMutations(sim.mutationsOfType(m2));
	femaleGenomes.addNewDrawnMutation(m1, 0); //add driver

        maleGenomes = p1.individuals[p1.individuals.sex == "M"].genomes;
        maleGenomes = maleGenomes[!maleGenomes.isNullGenome][0:49];
	maleGenomes.removeMutations(sim.mutationsOfType(m2));
	maleGenomes.addNewDrawnMutation(m1, 0); //add driver\n"""
    else:
        out += """
	p1.genomes[!p1.genomes.isNullGenome].addNewDrawnMutation(m2, 0); //add wildtype

        femaleGenomes = p1.individuals[p1.individuals.sex == "F"].genomes[0:99];
	femaleGenomes.removeMutations(sim.mutationsOfType(m2));
	femaleGenomes.addNewDrawnMutation(m1, 0); //add driver

        maleGenomes = p1.individuals[p1.individuals.sex == "M"].genomes[0:99];
	maleGenomes.removeMutations(sim.mutationsOfType(m2));
	maleGenomes.addNewDrawnMutation(m1, 0); //add driver\n"""
    out += """}

modifyChild() {
	driver = sim.mutationsOfType(m1);
	wild = sim.mutationsOfType(m2);
	r01 = sim.mutationsOfType(m3);
	r10 = sim.mutationsOfType(m4);
	r11 = sim.mutationsOfType(m5);
	MomGenomes = c(parent1Genome1, parent1Genome2);
	DadGenomes = c(parent2Genome1, parent2Genome2);
	ChildGenomes = c(childGenome1, childGenome2);
        ChildGenomes = ChildGenomes[!ChildGenomes.isNullGenome];
	child.tagF = 0;
	if(any(MomGenomes.containsMutations(driver))) { //Embryo
		c = """ + str(args.c1) + """;
		if(sum(MomGenomes.containsMutations(driver)) > 1){
			c = """ + str(args.c2) + """;
		}\n"""
    if args.two_allele:
        out += """		if(any(ChildGenomes.containsMutations(r10))) {
				if(sum(childGenome1.containsMutations(r10))){
                                        if(runif(1) < c) {		//convert to r11
					childGenome1.removeMutations(r10);
					childGenome1.addNewDrawnMutation(m5, 1);\n"""
        if args.X:
            out += "                }  } if(child.sex == 'F') {\n              if(sum(childGenome2.containsMutations(r10))){\n"
        else:
            out += "                } } if(sum(childGenome2.containsMutations(r10))){\n"
        out += """              if(runif(1) < c) {		//convert to r11
                                        childGenome2.removeMutations(r10);
					childGenome2.addNewDrawnMutation(m5, 1);
                                        }
			}\n"""
        if args.X:
            out += "        }\n"
        out += """		} if(any(ChildGenomes.containsMutations(r01))) {
				if(sum(childGenome1.containsMutations(r01))){
                                    if(runif(1) < c) {		//convert to r11
					childGenome1.removeMutations(r01);
					childGenome1.addNewDrawnMutation(m5, 1);\n"""
        if args.X:
            out += "                } } if(child.sex == 'F') {\n      if(sum(childGenome2.containsMutations(r01))){\n"
        else:
            out += "                } } if(sum(childGenome2.containsMutations(r01))){\n"
        out += """                  if(runif(1) < c) {		//convert to r11
                                    childGenome2.removeMutations(r01);
					childGenome2.addNewDrawnMutation(m5, 1);
                                    }
			}\n"""
        if args.X:
            out += "                }\n"
    if args.two_allele:
        out += """                  } else if(any(ChildGenomes.containsMutations(wild))) {\n"""
    else:
        out += """                  if(any(ChildGenomes.containsMutations(wild))) {\n"""
    out += """                  roll = runif(1);\n
                      if(sum(childGenome1.containsMutations(wild))){\n"""
    if args.two_allele:
        out += """		if(roll < c*(1-c)) {		//convert to r01\n"""
    else:
        out += """		if(roll < c) {		//convert to r01\n"""
    out += """
					childGenome1.removeMutations(wild);
					childGenome1.addNewDrawnMutation(m3, 1);\n"""
    if args.X:
        out += """                }
                                }
                                if(child.sex == 'F'){
                                if((runif(1) < c) & sum(childGenome2.containsMutations(wild))){\n"""
    elif args.two_allele:
        out += """                }
                                }
                                if((runif(1) < c*(1-c)) & sum(childGenome2.containsMutations(wild))){\n""" 
    else:
        out += """                }
                                }
                                if((runif(1) < c) & sum(childGenome2.containsMutations(wild))){\n""" 
    out += """          childGenome2.removeMutations(wild);
					childGenome2.addNewDrawnMutation(m3, 1);\n"""
    if args.X:
        out += """              }\n;"""
    if args.two_allele:
        out += """      }            else if(roll < 2*c*(1-c)) {		//convert to r10
				if(sum(childGenome1.containsMutations(wild))){
					childGenome1.removeMutations(wild);
					childGenome1.addNewDrawnMutation(m4, 1);\n"""
        if args.X:
            out += "                } if(child.sex == 'F') {\n     if((runif(1) < c*(1-c)) & sum(childGenome2.containsMutations(wild))){\n"
        else:
            out += "                } if((runif(1) < c*(1-c)) & sum(childGenome2.containsMutations(wild))){\n"
        out += """              childGenome2.removeMutations(wild);
					childGenome2.addNewDrawnMutation(m4, 1);
				}\n"""
        if args.X:
            out += "                }\n"
        out += """                  } else if(roll <2*c*(1-c)+c*c) {		//convert to r11
				if(sum(childGenome1.containsMutations(wild))){
					childGenome1.removeMutations(wild);
					childGenome1.addNewDrawnMutation(m5, 1);\n"""
        if args.X:
            out += "                } if(child.sex == 'F') {\n          if((runif(1) < c*c) & sum(childGenome2.containsMutations(wild))){\n"
        else:
            out += "                } if((runif(1) < c*c) & sum(childGenome2.containsMutations(wild))){\n"
        out += """                  childGenome2.removeMutations(wild);
					childGenome2.addNewDrawnMutation(m5, 1);
				}
			}
		}
	}\n"""
        if args.X:
            out += "                }\n"
    else:
        out += """ }
        }
    }\n"""
    out += """                  if(any(ChildGenomes.containsMutations(driver))) { //Germline\n"""
    if args.two_allele:
        out += """          	if(any(ChildGenomes.containsMutations(r10))) {
                        child.tagF = """ + str(1+args.Sr01D) + """;
                        roll = runif(1);
                         	if(roll < """ + str(args.b) + """) {		//convert to r11
				if(sum(childGenome1.containsMutations(r10))){
					childGenome1.removeMutations(r10);
					childGenome1.addNewDrawnMutation(m5, 1);\n"""
        if args.X:
            out += "                } else if((child.sex == 'F') & sum(childGenome2.containsMutations(r10))){\n"
        else:
            out += "                } else if(sum(childGenome2.containsMutations(r10))){\n"
        out += """                  childGenome2.removeMutations(r10);
					childGenome2.addNewDrawnMutation(m5, 1);
				}
			} else if(roll < """ + str(args.b) + """+""" +str(args.d) + """*(1-""" + str(args.b) + """)) {		//convert to D
				if(sum(childGenome1.containsMutations(r10))){
					childGenome1.removeMutations(r10);
					childGenome1.addNewDrawnMutation(m1, 1);\n"""
        if args.X:
            out += "                } else if((child.sex == 'F') & sum(childGenome2.containsMutations(r10))){\n"
        else:
            out += "                } else if(sum(childGenome2.containsMutations(r10))){\n"
        out += """                  childGenome2.removeMutations(r10);
					childGenome2.addNewDrawnMutation(m1, 1);
				}
			}
        		} else if(any(ChildGenomes.containsMutations(r01))) {
                  child.tagF = """ + str(1+args.Sr01D) + """;
                        roll = runif(1);
            	if(roll < """ + str(args.b) + """) {		//convert to r11
				if(sum(childGenome1.containsMutations(r01))){
					childGenome1.removeMutations(r01);
					childGenome1.addNewDrawnMutation(m5, 1);\n"""
        if args.X:
            out += "                } else if((child.sex == 'F') & sum(childGenome2.containsMutations(r01))){\n"
        else:
            out += "                } else if(sum(childGenome2.containsMutations(r01))){\n"
        out += """                  childGenome2.removeMutations(r01);
					childGenome2.addNewDrawnMutation(m5, 1);
				}
			} else if(roll <""" + str(args.b) + """+""" + str(args.d) + """*(1- """ + str(args.b) + """)) {		//convert to D
				if(sum(childGenome1.containsMutations(r01))){
					childGenome1.removeMutations(r01);
					childGenome1.addNewDrawnMutation(m1, 1);\n"""
        if args.X:
            out += "                } else if((child.sex == 'F') & sum(childGenome2.containsMutations(r01))){\n"
        else:
            out += "                } else if(sum(childGenome2.containsMutations(r01))){\n"
        out += """                  childGenome2.removeMutations(r01);
					childGenome2.addNewDrawnMutation(m1, 1);
			}
                    }\n"""
    prob = 0
    if args.two_allele:
        prob = args.b*(1-args.d)*(1-args.b)
        out += """		} else if(any(ChildGenomes.containsMutations(wild))) {\n"""
    else:
        prob = args.b
        out += """		if(any(ChildGenomes.containsMutations(wild))) {\n"""
    out += """
                        child.tagF = """ + str(1+args.SwD) + """;
                        roll = runif(1);
			if(roll < """ + str(prob) + """){		//convert to r01
				if(sum(childGenome1.containsMutations(wild))){
					childGenome1.removeMutations(wild);
					childGenome1.addNewDrawnMutation(m3, 1);\n"""
    if args.X:
        out += "                } else if((child.sex == 'F') & sum(childGenome2.containsMutations(wild))){\n"
    else:
        out += "                } else if(sum(childGenome2.containsMutations(wild))){\n"
    out += """                  childGenome2.removeMutations(wild);
					childGenome2.addNewDrawnMutation(m3, 1);
				}\n"""
    if args.two_allele:
        out += """                  	} else if(roll <""" + str(args.b) +  """*(1-""" + str(args.d) + """)*(1-""" + str(args.b) + """)+""" + str(args.b) +  """*(1-""" + str(args.d) + """)*(1-""" + str(args.b) + """)){		//convert to r10
				if(sum(childGenome1.containsMutations(wild))){
					childGenome1.removeMutations(wild);
					childGenome1.addNewDrawnMutation(m4, 1);\n"""
        if args.X:
            out += "                } else if((child.sex == 'F') & sum(childGenome2.containsMutations(wild))){\n"
        else:
            out += "                } else if(sum(childGenome2.containsMutations(wild))){\n"
        out += """                  childGenome2.removeMutations(wild);
					childGenome2.addNewDrawnMutation(m4, 1);
                                    }
			} else if(roll <""" + str(args.b) +  """*(1-""" + str(args.d) + """)*(1-""" + str(args.b) + """)+""" + str(args.b) +  """*(1-""" + str(args.d) + """)*(1-""" + str(args.b) + """+""" + str(args.b) +  """*""" + str(args.b) +  """)){			//convert to r11
				if(sum(childGenome1.containsMutations(wild))){
					childGenome1.removeMutations(wild);
					childGenome1.addNewDrawnMutation(m5, 1);\n"""
        if args.X:
            out += "                } else if((child.sex == 'F') & sum(childGenome2.containsMutations(wild))){\n"
        else:
            out += "                } else if(sum(childGenome2.containsMutations(wild))){\n"
        out += """                  childGenome2.removeMutations(wild);
					childGenome2.addNewDrawnMutation(m5, 1);
				}\n"""
    out += """			} else if(runif(1) < """ +str(args.d) + """) {		//convert to D
                      if(sum(childGenome1.containsMutations(wild))){
					childGenome1.removeMutations(wild);
					childGenome1.addNewDrawnMutation(m1, 1);\n"""
    if args.X:
        out += "                } else if((child.sex == 'F') & sum(childGenome2.containsMutations(wild))){\n"
    else:
        out += "                } else if(sum(childGenome2.containsMutations(wild))){\n"
    out += """                  childGenome2.removeMutations(wild);
					childGenome2.addNewDrawnMutation(m1, 1);
				}
			}
		}
	}
	return T;
}

fitness(NULL) {
	if(individual.tagF > 0)
		return individual.tagF;
	else
		return relFitness;
}

1:40 late() {
	if (sim.countOfMutationsOfType(m1) == 0){
                fixed = (sum(sim.substitutions.mutationType == m1) == 1);
                if(fixed){
			cat(paste('#OUTPUT: 1 ' + sim.generation, '\\n'));
                } else {
			cat(paste('#OUTPUT: 0 ' + sim.generation, '\\n'));
                }
		sim.simulationFinished();
	} else {
                driver = sim.mutationsOfType(m1);
		wild = sim.mutationsOfType(m2);
	     	r01 = sim.mutationsOfType(m3);
                r10 = sim.mutationsOfType(m4);
                r11 = sim.mutationsOfType(m5);
                cat('#OUTPUT: ');
		cat(sum(sim.mutationFrequencies(p1, driver)) + " " + sum(sim.mutationFrequencies(p1, wild)) + " " + sum(sim.mutationFrequencies(p1, r01)) + " " +  sum(sim.mutationFrequencies(p1, r10)) + " " + sum(sim.mutationFrequencies(p1, r11)) + '\\n');
	}
}
"""
    return out


def main(args):
    output = init(args)
    print(output)


if __name__=="__main__":
    """parser = ArgumentParser(description='Process some integers.')
    parser.add_argument('-N', dest='N', type=int, default=DEFAULT["N"], help="Census Population Size, default=" + str(DEFAULT["N"]))
    parser.add_argument('-Ne', dest='Ne', type=int, default=DEFAULT["Ne"], help="Variance effective population size, default=" + str(DEFAULT["N"]))
    parser.add_argument('-sd0', dest='Sd0', type=float, default=DEFAULT["Sd0"], help="Fitness cost of driver/wildtype heterozygotes, default=" + str(DEFAULT["Sd0"]))
    parser.add_argument('-sdr', dest='Sdr', type=float, default=DEFAULT["Sdr"], help="Fitness cost of driver/resistance heterozygotes, default=" + str(DEFAULT["Sdr"]))
    parser.add_argument('-sr0', dest='Sr0', type=float, default=DEFAULT["Sr0"], help="Fitness cost of resistance/wildtype heterozygotes, default=" + str(DEFAULT["Sr0"]))
    parser.add_argument('-sdd', dest='Sdd', type=float, default=DEFAULT["Sdd"], help="Fitness cost of driver homozygotes, default=" + str(DEFAULT["Sdd"]))
    parser.add_argument('-srr', dest='Srr', type=float, default=DEFAULT["Srr"], help="Fitness cost of resistance homozygotes, default=" + str(DEFAULT["Srr"]))
    parser.add_argument('-L', dest='L', type=int, default=DEFAULT["L"], help="Chromosome Length=, default=" + str(DEFAULT["L"]))
    parser.add_argument('-c', dest='C', type=float, default=DEFAULT["C"], help="MCR converseion rate, default=" + str(DEFAULT["C"]))
    parser.add_argument('-r', dest='R', type=float, default=DEFAULT["R"], help="Recombination rate, default=" + str(DEFAULT["R"]))
    parser.add_argument('-d', dest='D', type=float, default=DEFAULT["D"], help="Fraction of cases in which repair generates resistance allele by NHEJ, default=" + str(DEFAULT["D"]))
    parser.add_argument('-mu', dest='MU', type=float, default=DEFAULT["MU"], help="Rate at which wildtype allele mutates into resistance allele, default=" + str(DEFAULT["MU"]))
    parser.add_argument('-G', dest='G', type=int, default=DEFAULT["G"], help="Number of generations, default=" + str(DEFAULT["G"]))
    parser.add_argument('-xd', dest='Xd', type=float, default=DEFAULT["Xd"], help="Introduction frequency of driver allele, default=" + str(DEFAULT["Xd"]))
    parser.add_argument('-xr', dest='Xr', type=float, default=DEFAULT["Xr"], help="Introduction frequency of resistance allele, default=" + str(DEFAULT["Xr"]))
    parser.add_argument('-sb', dest='SB', action='store_true', default=DEFAULT["DriveSexBias"], help="Flag to signify that the drive only occures in one sex")
    parser.add_argument('-hl', dest='HL', action='store_true', default=DEFAULT["HL"], help="Flag to signify that after encountering a driver gamete, the wildtype is homoleathal")
    parser.add_argument('-gamete', dest='Gamete', action='store_true', default=DEFAULT["Gamete"], help="Flag to signify that homing occurs in gametes")
    args = parser.parse_args()"""
    args = params()
    main(args)
