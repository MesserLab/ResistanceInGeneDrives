#! /usr/bin/env python
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import os

def graph(filename, drive, location, two=True):
    df = pd.read_table(filename, header=None,  sep=" ")
    df.loc[df.ix[:,40] == 0,"id"] = "Driver"
    df.loc[df.ix[:,40] == 1,"id"] = "Wild"
    df.loc[df.ix[:,40] == 2,"id"] = "r01"
    df.loc[df.ix[:,40] == 3,"id"] = "r10"
    df.loc[df.ix[:,40] == 4,"id"] = "r11"
    df = df.drop(df.columns[40], axis=1)

    font = {'family' : 'sans-serif', 'size' : 11}
    plt.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rc('font', **font)

    fig, ax1 = plt.subplots(1)
    plt.setp(ax1, xticks=np.arange(0, 41, 10.0))
    fig.set_size_inches(3, 4)
    mean = df.groupby("id").mean().transpose()
    print("max frequency: " + str(mean.Driver.max()))
    print("generation of max frequency: " + str(mean.Driver.idxmax()))
    std_dev = df.groupby("id").std().transpose()
    t = np.arange(40)
    ax1.plot(mean.Driver, label="Driver", color="#ff328b")
    ax1.fill_between(t, mean.Driver-std_dev.Driver, \
            mean.Driver+std_dev.Driver, alpha=0.3,\
            color="#ff328b")

    ax1.plot(mean.Wild, label="+", color="#7384ff")
    ax1.fill_between(t, mean.Wild-std_dev.Wild, \
            mean.Wild+std_dev.Wild, alpha=0.3,\
            facecolor="#7384ff")

    if two:
        ax1.plot(mean.r01 + mean.r10, label="r01", color="#6e3100")
        lower_bound = mean.r01+ mean.r10-\
                (std_dev.r01**2+std_dev.r10**2)**(1/2.0)

        upper_bound = mean.r01+ mean.r10+\
                (std_dev.r01**2+std_dev.r10**2)**(1/2.0)

        ax1.fill_between(t, upper_bound, lower_bound, alpha=0.3,\
                facecolor="#6e3100")

        ax1.plot(mean.r11, label="r11", color="#008454")
        ax1.fill_between(t, mean.r11-std_dev.r11,\
                mean.r11+std_dev.r11, \
                alpha=0.3, facecolor="#008454")
    else:
        ax1.plot(mean.r01 + mean.r10, label="r01", color="#008454")
        lower_bound = mean.r01+ mean.r10-\
                (std_dev.r01**2+std_dev.r10**2)**(1/2.0)

        upper_bound = mean.r01+ mean.r10+\
                (std_dev.r01**2+std_dev.r10**2)**(1/2.0)

        ax1.fill_between(t, lower_bound, upper_bound,\
                alpha=0.3, facecolor="#008454")

        ax1.spines['right'].set_color('none')
    ax1.spines['top'].set_color('none')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    ax1.set_ylim(0, 1)
    ax1.set_xlim(0, 40)
    ax1.set_xlabel('Generations')
    ax1.set_ylabel('Frequency')
    ax1.legend( loc=location, prop={"size":10}, frameon=False)
    filename = filename.split('/')[3]
    print(drive + ' ' + filename)
    plt.tight_layout()
    fig.savefig('Figures/' + drive + '/' + filename.split('.')[0] + '.png', format='png', dpi=1500)

def main(args):
    drives = ["High_resistance", "Medium_resistance", "Low_resistance"]
    legend_location = {("High_resistance", "Autosomal_2gRNA.txt"):"upper right",
            ("High_resistance", "Autosomal_1gRNA.txt"):"upper right",
            ("High_resistance", "X_2gRNA.txt"):"upper right",
            ("High_resistance", "X_1gRNA.txt"):"upper right",
            ("Medium_resistance", "Autosomal_2gRNA.txt"):"center right",
            ("Medium_resistance", "Autosomal_1gRNA.txt"):"upper right",
            ("Medium_resistance", "X_2gRNA.txt"):"center right",
            ("Medium_resistance", "X_1gRNA.txt"):"upper right",
            ("Low_resistance", "Autosomal_2gRNA.txt"):"center right",
            ("Low_resistance", "Autosomal_1gRNA.txt"):"upper right",
            ("Low_resistance", "X_2gRNA.txt"):"center right",
            ("Low_resistance", "X_1gRNA.txt"):"upper right"}

    location = 'Data/'
    if args.example:
        location = 'Data/Example'
    for drive in drives:
        for filename in os.listdir('../Data/' + drive + '/'):
            if filename.split('_')[1].startswith('2'):
                graph(location + drive + '/' + filename, drive,\
                        legend_location[(drive, filename)])
            else:
                graph(location + drive + '/' + filename, drive,\
                        legend_location[(drive, filename)], False)

if __name__=="__main__":
    parser = ArgumentParser(description='Generate Figures.')
    parser.add_argument('-example', dest='example', action='store_true', default=False, help="Flag to generate figures from example data")
    args = parser.parse_args()
    main(args)
