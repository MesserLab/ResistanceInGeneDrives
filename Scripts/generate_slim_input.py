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
        self.N = N
        self.SD=SD
        self.Sr=Sr
        self.hD=hD
        self.hr=hr
        self.b=b
        self.d=d
        self.c1=c1
        self.c2=c2
        self.X = X
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
    child.tagF = 0;
    if(any(MomGenomes.containsMutations(sim.mutationsOfType(m1)))) { //Embryo
            c = """ + str(args.c1) + """;
            if(sum(MomGenomes.containsMutations(sim.mutationsOfType(m1))) > 1){
                c = """ + str(args.c2) + """;
            }

            for(cGenome in ChildGenomes){
                    if(sum(cGenome.containsMutations(sim.mutationsOfType(m1)))==0){
                            if(sum(cGenome.containsMutations(sim.mutationsOfType(m3))) < """ + str(args.N_guideRNA) + """) {
                                    for (i in 1:(""" + str(args.N_guideRNA) + """-sum(cGenome.containsMutations(sim.mutationsOfType(m3))))) {
                                            if(runif(1) < c) {		//convert to r11
                                                    cGenome.addNewDrawnMutation(m3, 1);
                                            }
                                    }
                            }
                    }
            }
    }

    if(any(ChildGenomes.containsMutations(sim.mutationsOfType(m1)))) { //Germline\n"""

    if args.X:
        out += """  if(child.sex == "F") {\n"""
    out += """  for(cGenome in ChildGenomes) {
                    if(sum(cGenome.containsMutations(sim.mutationsOfType(m1))) == 0){
                        if(sum(cGenome.containsMutations(sim.mutationsOfType(m3))) <""" + str(args.N_guideRNA) + """) {
                            for (i in 1:(""" + str(args.N_guideRNA) + """-sum(cGenome.containsMutations(sim.mutationsOfType(m3))))) {
                                if(runif(1) <""" +str(args.b) + """) {		//convert to r11
                                    cGenome.addNewDrawnMutation(m3, 1);
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
