#!/usr/bin/env python
import generate_slim_input
from admin import run_slim
from multiprocessing import Pool
import sys
RUNS = 1000

def main(cmd_line):
    args = generate_slim_input.params()
    args.N_guideRNA = 1
    #medium
    args.c1 = 0.23
    args.c2 = 0.32
    args.b = 0.3
    args.d = 0.98
    output = generate_slim_input.init(args)
    with open("Data/tmp/medium_one.txt", "w+") as f: #write SliM input file
            f.write(output)

    #high
    args.c1 = 0.75
    args.c2 = 0.95
    args.b = 0.35
    args.d = 0.99
    output = generate_slim_input.init(args)
    with open("Data/tmp/high_one.txt", "w+") as f: #write SliM input file
            f.write(output)

    #low
    args.c1 = 0.03
    args.c2 = 0.05
    args.b = 0.25
    args.d = 0.97
    output = generate_slim_input.init(args)
    with open("Data/tmp/low_one.txt", "w+") as f: #write SliM input file
            f.write(output)

    #erase previous data
    open('Data/Medium_resistance/autosomal_1gRNA.txt', 'w+').close()
    open('Data/High_resistance/autosomal_1gRNA.txt', 'w+').close()
    open('Data/Low_resistance/autosomal_1gRNA.txt', 'w+').close()
    for j in range(RUNS):
        #run simulations in parallel
        with Pool(processes=3) as pool:
            medium_pool = pool.apply_async(run_slim, ["Data/tmp/medium_one.txt"])
            high_pool = pool.apply_async(run_slim, ["Data/tmp/high_one.txt"])
            low_pool = pool.apply_async(run_slim, ["Data/tmp/low_one.txt"])

            medium_result = medium_pool.get()
            high_result = high_pool.get()
            low_result = low_pool.get()

            # write to file
            with open('Data/Medium_resistance/autosomal_1gRNA.txt', "a") as f:
                for i in range(0, len(medium_result)):
                    f.write(" ".join(medium_result[i]) + " " + str(i) + "\n")

            with open('Data/High_resistance/autosomal_1gRNA.txt', "a") as f:
                for i in range(0, len(high_result)):
                    f.write(" ".join(high_result[i]) + " " + str(i) + "\n")

            with open('Data/Low_resistance/autosomal_1gRNA.txt', "a") as f:
                for i in range(0, len(low_result)):
                    f.write(" ".join(low_result[i]) + " " + str(i) + "\n")


if __name__=="__main__":
	main(sys.argv)
