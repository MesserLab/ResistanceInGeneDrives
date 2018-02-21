#! /usr/bin/env python
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import os

def set_ax(filename, ax, two=True):
    ax.clear()
    df = pd.read_table(filename, header=None,  sep=" ")
    df.loc[df.ix[:,41] == 0,"id"] = "Driver"
    df.loc[df.ix[:,41] == 1,"id"] = "Wild"
    df.loc[df.ix[:,41] == 2,"id"] = "r1"
    df.loc[df.ix[:,41] == 3,"id"] = "r2"
    df = df.drop(df.columns[41], axis=1)
    mean = df.groupby("id").mean().transpose()
    print("Driver: " + ",".join(str(i) for i in mean.Driver))
    print("Wild: " + ",".join(str(i) for i in mean.Wild))
    print("r1: " + ",".join(str(i) for i in mean.r1))
    print("r2: " + ",".join(str(i) for i in mean.r2))
    print("\n")
    std_dev = df.groupby("id").std().transpose()

    t = np.arange(41)
    plt.setp(ax, xticks=np.arange(0, 42, 10.0))
    ax.plot(mean.Driver, label="Driver", color="#ff328b", linewidth=2)
    ax.fill_between(t, mean.Driver-std_dev.Driver, \
            mean.Driver+std_dev.Driver, alpha=0.3,\
            color="#ff328b")

    ax.plot(mean.Wild, label="Wild Type", color="#7384ff", linewidth=2)
    ax.fill_between(t, mean.Wild-std_dev.Wild, \
            mean.Wild+std_dev.Wild, alpha=0.3,\
            facecolor="#7384ff")

    if two:
        ax.plot(mean.r1, label="(r,+)", color="#6e3100", linewidth=2)
        lower_bound = mean.r1-\
                (std_dev.r1)

        upper_bound = mean.r1+\
                (std_dev.r1)**(1/2.0)

        ax.fill_between(t, upper_bound, lower_bound, alpha=0.3,\
                facecolor="#6e3100")

        ax.plot(mean.r2, label="Resistance", color="#008454", linewidth=2)
        ax.fill_between(t, mean.r2-std_dev.r2,\
                mean.r2+std_dev.r2, \
                alpha=0.3, facecolor="#008454")
    else:
        ax.plot(mean.r1, label="r01", color="#008454", linewidth=2)
        lower_bound = mean.r1-\
                (std_dev.r1)

        upper_bound = mean.r1+\
                (std_dev.r1)

        ax.fill_between(t, lower_bound, upper_bound,\
                alpha=0.3, facecolor="#008454")

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 40)
    ax.tick_params(length=20, width=2)

def main(args):
    font = {'family' : 'sans-serif', 'size' : 11}
    plt.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rc('font', **font)

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)

    location = 'Data/'
    if args.example:
        location = 'Data/Example/'
    print("X_1gRNA_high")
    set_ax(location + "High_resistance/X_1gRNA.txt", ax, False)
    fig.tight_layout()
    fig.savefig('Figures/X_1_H.pdf', format='pdf', dpi=1500)

    print("X_2gRNA_high")
    set_ax(location + "High_resistance/X_2gRNA.txt", ax)
    fig.tight_layout()
    fig.savefig('Figures/X_2_H.pdf', format='pdf', dpi=1500)

    print("autosomal_1gRNA_high")
    set_ax(location + "High_resistance/autosomal_1gRNA.txt", ax, False)
    fig.tight_layout()
    fig.savefig('Figures/A_1_H.pdf', format='pdf', dpi=1500)

    print("autosomal_2gRNA_high")
    set_ax(location + "High_resistance/autosomal_2gRNA.txt", ax)
    fig.tight_layout()
    fig.savefig('Figures/A_2_H.pdf', format='pdf', dpi=1500)

    print("X_1gRNA_medium")
    set_ax(location + "Medium_resistance/X_1gRNA.txt", ax, False)
    fig.tight_layout()
    fig.savefig('Figures/X_1_M.pdf', format='pdf', dpi=1500)

    print("X_2gRNA_medium")
    set_ax(location + "Medium_resistance/X_2gRNA.txt", ax)
    fig.tight_layout()
    fig.savefig('Figures/X_2_M.pdf', format='pdf', dpi=1500)

    print("autosomal_1gRNA_medium")
    set_ax(location + "Medium_resistance/autosomal_1gRNA.txt", ax, False)
    fig.tight_layout()
    fig.savefig('Figures/A_1_M.pdf', format='pdf', dpi=1500)

    print("autosomal_2gRNA_medium")
    set_ax(location + "Medium_resistance/autosomal_2gRNA.txt", ax)
    fig.tight_layout()
    fig.savefig('Figures/A_2_M.pdf', format='pdf', dpi=1500)

    print("X_1gRNA_low")
    set_ax(location + "Low_resistance/X_1gRNA.txt", ax, False)
    fig.tight_layout()
    fig.savefig('Figures/X_1_L.pdf', format='pdf', dpi=1500)

    print("X_2gRNA_low")
    set_ax(location + "Low_resistance/X_2gRNA.txt", ax)
    fig.tight_layout()
    fig.savefig('Figures/X_2_L.pdf', format='pdf', dpi=1500)

    print("autosomal_1gRNA_low")
    set_ax(location + "Low_resistance/autosomal_1gRNA.txt", ax, False)
    fig.tight_layout()
    fig.savefig('Figures/A_1_L.pdf', format='pdf', dpi=1500)

    print("autosomal_2gRNA_low")
    set_ax(location + "Low_resistance/autosomal_2gRNA.txt", ax)
    fig.tight_layout()
    fig.savefig('Figures/A_2_L.pdf', format='pdf', dpi=1500)

if __name__=="__main__":
    parser = ArgumentParser(description='Generate Figures.')
    parser.add_argument('-example', dest='example', action='store_true', default=False, help="Flag to generate figures from example data")
    args = parser.parse_args()
    main(args)
