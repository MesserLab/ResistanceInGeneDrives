#!/usr/bin/env python
import subprocess
from random import uniform, randint, seed
import generate_slim_input
from parse import parse
from multiprocessing import Pool
import sys
seed()
ID = randint(1,1000)
RUNS = 10000

def run_slim(filename):
        process = subprocess.Popen(["slim", filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
        out, err = process.communicate()
        return parse(out)


def main(cmd_line):
    args = generate_slim_input.params()
    #medium
    args.c1 = 0.2
    args.c2 = 0.27
    args.b = 0.28
    args.d = 0.94
    output = generate_slim_input.init(args) #HD F
    with open("Data/tmp/X_medium.txt", "w+") as f:
            f.write(output)

    #high
    args.c1 = 0.65
    args.c2 = 0.8
    args.b = 0.3
    args.d = 0.95
    output = generate_slim_input.init(args) #HD F
    with open("Data/tmp/X_high.txt", "w+") as f:
            f.write(output)

    #low
    args.c1 = 0.03
    args.c2 = 0.05
    args.b = 0.25
    args.d = 0.93
    output = generate_slim_input.init(args) #HD F
    with open("Data/tmp/X_low.txt", "w+") as f:
            f.write(output)

    open('Data/Medium_resistance/X_2gRNA.txt', 'w+').close()
    open('Data/High_resistance/X_2gRNA.txt', 'w+').close()
    open('Data/Low_resistance/X_2gRNA.txt', 'w+').close()
    for j in range(1000):
        with Pool(processes=4) as pool:
            medium_pool = pool.apply_async(run_slim, ["Data/tmp/X_medium.txt"])
            high_pool = pool.apply_async(run_slim, ["Data/tmp/X_high.txt"])
            low_pool = pool.apply_async(run_slim, ["Data/tmp/X_low.txt"])

            medium_result = medium_pool.get()
            high_result = high_pool.get()
            low_result = low_pool.get()

            with open('Data/Medium_resistance/X_2gRNA.txt', "a") as f:
                for i in range(0, len(medium_result)):
                    f.write(" ".join(medium_result[i]) + " " + str(i) + "\n")

            with open('Data/High_resistance/X_2gRNA.txt', "a") as f:
                for i in range(0, len(high_result)):
                    f.write(" ".join(high_result[i]) + " " + str(i) + "\n")

            with open('Data/Low_resistance/X_2gRNA.txt', "a") as f:
                for i in range(0, len(low_result)):
                    f.write(" ".join(low_result[i]) + " " + str(i) + "\n")

if __name__=="__main__":
	main(sys.argv)
