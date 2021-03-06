#!/usr/bin/env python
import generate_slim_input
from admin import run_slim
from multiprocessing import Pool
import sys
RUNS = 1000

def main(cmd_line):
    args = generate_slim_input.params()
    args.X = True
    #medium
    args.c1 = 0.2
    args.c2 = 0.27
    args.b = 0.28
    args.d = 0.94
    args.N_guideRNA = 2
    output = generate_slim_input.init(args)
    with open("Data/tmp/X_medium.txt", "w+") as f: #write SliM input file
            f.write(output)

    #high
    args.c1 = 0.65
    args.c2 = 0.8
    args.b = 0.3
    args.d = 0.95
    output = generate_slim_input.init(args)
    with open("Data/tmp/X_high.txt", "w+") as f: #write SliM input file
            f.write(output)

    #low
    args.c1 = 0.03
    args.c2 = 0.05
    args.b = 0.25
    args.d = 0.93
    output = generate_slim_input.init(args)
    with open("Data/tmp/X_low.txt", "w+") as f: #write SliM input file
            f.write(output)

    #erase previous data
    open('Data/Medium_resistance/X_2gRNA.txt', 'w+').close()
    open('Data/High_resistance/X_2gRNA.txt', 'w+').close()
    open('Data/Low_resistance/X_2gRNA.txt', 'w+').close()
    for j in range(RUNS):
        #run simulations in parallel
        with Pool(processes=3) as pool:
            medium_pool = pool.apply_async(run_slim, ["Data/tmp/X_medium.txt"])
            high_pool = pool.apply_async(run_slim, ["Data/tmp/X_high.txt"])
            low_pool = pool.apply_async(run_slim, ["Data/tmp/X_low.txt"])

            medium_result = medium_pool.get()
            high_result = high_pool.get()
            low_result = low_pool.get()

            # write to file
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
