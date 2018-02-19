#!/usr/bin/env python
import generate_slim_input
from admin import run_slim
from multiprocessing import Pool
import sys
RUNS = 1000

def main(cmd_line):
    args = generate_slim_input.params()
    args.N_guideRNA = 2
    #medium
    args.c1 = 0.2
    args.c2 = 0.27
    args.b = 0.28
    args.d = 0.94
    output = generate_slim_input.init(args) #HD F
    with open("Data/tmp/medium.txt", "w+") as f:
            f.write(output)

    #high
    args.c1 = 0.65
    args.c2 = 0.8
    args.b = 0.3
    args.d = 0.95
    output = generate_slim_input.init(args) #HD F
    with open("Data/tmp/high.txt", "w+") as f:
            f.write(output)

    #low
    args.c1 = 0.03
    args.c2 = 0.05
    args.b = 0.25
    args.d = 0.93
    output = generate_slim_input.init(args) #HD F
    with open("Data/tmp/low.txt", "w+") as f:
            f.write(output)

    open('Data/Medium_resistance/autosomal_2gRNA.txt', 'w+').close()
    open('Data/High_resistance/autosomal_2gRNA.txt', 'w+').close()
    open('Data/Low_resistance/autosomal_2gRNA.txt', 'w+').close()
    for j in range(RUNS):
        with Pool(processes=4) as pool:
            medium_pool = pool.apply_async(run_slim, ["Data/tmp/medium.txt"])
            high_pool = pool.apply_async(run_slim, ["Data/tmp/high.txt"])
            low_pool = pool.apply_async(run_slim, ["Data/tmp/low.txt"])

            medium_result = medium_pool.get()
            high_result = high_pool.get()
            low_result = low_pool.get()

            with open('Data/Medium_resistance/autosomal_2gRNA.txt', "a") as f:
                for i in range(0, len(medium_result)):
                    f.write(" ".join(medium_result[i]) + " " + str(i) + "\n")

            with open('Data/High_resistance/autosomal_2gRNA.txt', "a") as f:
                for i in range(0, len(high_result)):
                    f.write(" ".join(high_result[i]) + " " + str(i) + "\n")

            with open('Data/Low_resistance/autosomal_2gRNA.txt', "a") as f:
                for i in range(0, len(low_result)):
                    f.write(" ".join(low_result[i]) + " " + str(i) + "\n")

if __name__=="__main__":
	main(sys.argv)
