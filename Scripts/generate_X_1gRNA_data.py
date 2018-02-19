#!/usr/bin/env python
import generate_slim_input
from admin import run_slim
from multiprocessing import Pool
import sys
RUNS = 1000

def main(cmd_line):
    args = generate_slim_input.params()
    args.N_guideRNA = 1
    args.X = True
    #medium
    args.c1 = 0.23
    args.c2 = 0.32
    args.b = 0.3
    args.d = 0.98
    output = generate_slim_input.init(args) #HD F
    with open("Data/tmp/X_medium_one.txt", "w+") as f:
            f.write(output)

    #high
    args.c1 = 0.75
    args.c2 = 0.95
    args.b = 0.35
    args.d = 0.99
    output = generate_slim_input.init(args) #HD F
    with open("Data/tmp/X_high_one.txt", "w+") as f:
            f.write(output)

    #low
    args.c1 = 0.03
    args.c2 = 0.05
    args.b = 0.25
    args.d = 0.97
    output = generate_slim_input.init(args) #HD F
    with open("Data/tmp/X_low_one.txt", "w+") as f:
            f.write(output)

    open('Data/Medium_resistance/X_1gRNA.txt', 'w+').close()
    open('Data/High_resistance/X_1gRNA.txt', 'w+').close()
    open('Data/Low_resistance/X_1gRNA.txt', 'w+').close()
    for j in range():
        with Pool(processes=4) as pool:
            X_medium_pool = pool.apply_async(run_slim, ["Data/tmp/X_medium_one.txt"])
            X_high_pool = pool.apply_async(run_slim, ["Data/tmp/X_high_one.txt"])
            X_low_pool = pool.apply_async(run_slim, ["Data/tmp/X_low_one.txt"])

            X_medium_result = X_medium_pool.get()
            X_high_result = X_high_pool.get()
            X_low_result = X_low_pool.get()

            with open('Data/Medium_resistance/X_1gRNA.txt', "a") as f:
                for i in range(0, len(X_medium_result)):
                    f.write(" ".join(X_medium_result[i]) + " " + str(i) + "\n")

            with open('Data/High_resistance/X_1gRNA.txt', "a") as f:
                for i in range(0, len(X_high_result)):
                    f.write(" ".join(X_high_result[i]) + " " + str(i) + "\n")

            with open('Data/Low_resistance/X_1gRNA.txt', "a") as f:
                for i in range(0, len(X_low_result)):
                    f.write(" ".join(X_low_result[i]) + " " + str(i) + "\n")

if __name__=="__main__":
	main(sys.argv)
