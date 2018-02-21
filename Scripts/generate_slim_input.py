#!/usr/bin/env python
from argparse import ArgumentParser
from math import floor, sqrt

DEFAULT = {
        "N": 10000,
        "SD": -0.05,
        "Sr": 0.0,
        "hD": 0.5,
        "hr": 0.5,
        "b": 0.25,
        "d": 0.93,
        "c1":0.03,
        "c2":0.05,
        "X": False,
        "N_guideRNA": 2,
}

class params:
    def __init__(self, N=DEFAULT["N"], SD=DEFAULT["SD"], Sr=DEFAULT["Sr"],\
            hD=DEFAULT["hD"], hr=DEFAULT["hr"], d=DEFAULT["d"],\
            b=DEFAULT["b"], c1=DEFAULT["c1"], c2=DEFAULT["c2"],\
            X=DEFAULT["X"], N_guideRNA=DEFAULT["N_guideRNA"]):
        self.N = N # population size
        self.SD=SD # selection coefficient of drive allele
        self.Sr=Sr # selection coefficient of resistance allele
        self.hD=hD # dominance coefficient of drive allele
        self.hr=hr # dominance coefficient of resistance allele
        self.b=b   # chance that a resistance allele is formed early at a single site
        self.d=d   # chance that the allele is converted to drive (after resistance alleles form) based on one cut site
        self.c1=c1 # chance of r allele formation during embryo stage if the mother has one D allele
        self.c2=c2 # chance of r allele formation during embryo stage if the mother has two D alleles
        self.X = X # flag set to True if the drive is x-linked
        self.N_guideRNA = N_guideRNA

def init(args):
    out = """
initialize() {
    initializeMutationRate( 0.0 );
    initializeMutationType("m1", """ + str(args.hD) + """, "f", """ + str(args.SD) + """); //D
    initializeMutationType("m3", """ + str(args.hr) + """, "f", """ + str(args.Sr) + """); //r01
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

    femaleGenomes = p1.individuals[p1.individuals.sex == "F"].genomes[0:99];
    femaleGenomes.addNewDrawnMutation(m1, 0); //add driver

    maleGenomes = p1.individuals[p1.individuals.sex == "M"].genomes;
    maleGenomes = maleGenomes[!maleGenomes.isNullGenome][0:49];
    maleGenomes.addNewDrawnMutation(m1, 0); //add driver\n"""
    else:
        out += """
    femaleGenomes = p1.individuals[p1.individuals.sex == "F"].genomes[0:99];
    femaleGenomes.addNewDrawnMutation(m1, 0); //add driver

    maleGenomes = p1.individuals[p1.individuals.sex == "M"].genomes[0:99];
    maleGenomes.addNewDrawnMutation(m1, 0); //add driver\n"""

    out += """}

modifyChild() {
    MomGenomes = c(parent1Genome1, parent1Genome2);
    DadGenomes = c(parent2Genome1, parent2Genome2);
    ChildGenomes = c(childGenome1, childGenome2);
    ChildGenomes = ChildGenomes[!ChildGenomes.isNullGenome];

    // Add resistance allele during embryo stage
    if(any(MomGenomes.containsMutations(sim.mutationsOfType(m1)))) {
            c = """ + str(args.c1) + """;
            if(sum(MomGenomes.containsMutations(sim.mutationsOfType(m1))) > 1){
                c = """ + str(args.c2) + """;
            }

            for(cGenome in ChildGenomes){
                // if child has drive allele
                if(sum(cGenome.containsMutations(sim.mutationsOfType(m1)))==0){

                    // if there are wildtype alleles (not drive and not resistance
                    if(sum(cGenome.containsMutations(sim.mutationsOfType(m3))) < """ + str(args.N_guideRNA) + """) {

                        // for each wildtype allele
                        for (i in 1:(""" + str(args.N_guideRNA) + """-sum(cGenome.containsMutations(sim.mutationsOfType(m3))))) {

                            // convert to resistance allele
                            // with probability c
                            if(runif(1) < c) {		//convert to r11
                                cGenome.addNewDrawnMutation(m3, 1);
                            }
                        }
                    }
                }
            }
    }

    // germline stage
    if(any(ChildGenomes.containsMutations(sim.mutationsOfType(m1)))) {\n"""

    if args.X:
        out += """  if(child.sex == "F") {\n"""
    out += """  for(cGenome in ChildGenomes) {

                    // if child has drive allele
                    if(sum(cGenome.containsMutations(sim.mutationsOfType(m1))) == 0){

                        // if there are wildtype alleles (not drive and not resistance
                        if(sum(cGenome.containsMutations(sim.mutationsOfType(m3))) <""" + str(args.N_guideRNA) + """) {

                            // for each wildtype allele
                            for (i in 1:(""" + str(args.N_guideRNA) + """-sum(cGenome.containsMutations(sim.mutationsOfType(m3))))) {

                                // convert to resistance allele with
                                // probability b
                                if(runif(1) <""" +str(args.b) + """) {
                                    cGenome.addNewDrawnMutation(m3, 1);

                                // convert to drive allele with
                                // probability d
                                } else if(runif(1) <""" +str(args.d) + """) {
                                    cGenome.removeMutations(sim.mutationsOfType(m3));
                                    cGenome.addNewDrawnMutation(m1, 1);
                                    break;
                                }
                            }
                        }
                    }
                } \n"""
    if args.X:
        out += """}\n"""
    out += """        }
        return T;
}\n"""

    out += """
1:41 late() {
    // if drive allele is fixed or lost
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
        r01 = sim.mutationsOfType(m3);
        cat('#OUTPUT: ');
        wild = 1-sum(sim.mutationFrequencies(p1, driver));
        cat(sum(sim.mutationFrequencies(p1, driver)) + " ");
        for(i in 1:2){
            count = 0;
            for(p in p1.genomes[!p1.genomes.isNullGenome]) {
                if(sum(p.containsMutations(r01)) == i) {
                    count = count + 1;
                }
            }\n"""
    if args.X:
        out += """      cat(count/""" + str(1.5*args.N) + """+ " ");
            wild = wild - count/""" + str(1.5*args.N) + """;\n"""
    else:
        out += """      cat(count/""" + str(2*args.N) + """+ " ");
            wild = wild - count/""" + str(2*args.N) + """;\n"""
    out += """}
        cat(wild + \"\\n\");
	}\n"""
    if args.X:
        out += """\n"""
    out += """}
"""
    return out


def main(args):
    output = init(args)
    print(output)


if __name__=="__main__":
    parser = ArgumentParser(description='Process some integers.')
    parser.add_argument('-N', dest='N', type=int, default=DEFAULT["N"], help="Census Population Size, default=" + str(DEFAULT["N"]))
    parser.add_argument('-N_guideRNA', dest='N_guideRNA', type=int, default=DEFAULT["N_guideRNA"], help="Number of guideRNAs, default=" + str(DEFAULT["N_guideRNA"]))
    parser.add_argument('-sd', dest='SD', type=float, default=DEFAULT["SD"], help="Fitness cost of driver homozygote, default=" + str(DEFAULT["SD"]))
    parser.add_argument('-sr', dest='Sr', type=float, default=DEFAULT["Sr"], help="Fitness cost of resistance homozygote, default=" + str(DEFAULT["Sr"]))
    parser.add_argument('-hd', dest='hD', type=float, default=DEFAULT["hD"], help="Dominance coefficient of drive allele, default=" + str(DEFAULT["hD"]))
    parser.add_argument('-hr', dest='hr', type=float, default=DEFAULT["hr"], help="Dominance coefficient of resistance allele, default=" + str(DEFAULT["hr"]))
    parser.add_argument('-c1', dest='c1', type=float, default=DEFAULT["c1"], help="chance of resistance allele formation in embryo if the mother has one D allele, default=" + str(DEFAULT["c1"]))
    parser.add_argument('-c2', dest='c2', type=float, default=DEFAULT["c2"], help="chance of resistance allele formation in embryo if the mother has two D allele, default=" + str(DEFAULT["c1"]))
    parser.add_argument('-b', dest='b', type=float, default=DEFAULT["b"], help="Fraction of cases in which repair generates resistance allele by NHEJ, default=" + str(DEFAULT["b"]))
    parser.add_argument('-d', dest='d', type=float, default=DEFAULT["d"], help="Chance that the allele is converted to drive (after resistance alleles form) based on one cut site, default=" + str(DEFAULT["d"]))
    parser.add_argument('-X', dest='X', action='store_true', default=DEFAULT["X"], help="Flag to signify that the drive is x-linked, default=" + str(DEFAULT["X"]))
    args = parser.parse_args()
    args = params()
    main(args)
