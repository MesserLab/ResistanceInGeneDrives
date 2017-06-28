#! /usr/bin/env python
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import os

def set_ax(filename, ax, two=True):
    df = pd.read_table(filename, header=None,  sep=" ")
    df.loc[df.ix[:,40] == 0,"id"] = "Driver"
    df.loc[df.ix[:,40] == 1,"id"] = "Wild"
    df.loc[df.ix[:,40] == 2,"id"] = "r01"
    df.loc[df.ix[:,40] == 3,"id"] = "r10"
    df.loc[df.ix[:,40] == 4,"id"] = "r11"
    df = df.drop(df.columns[40], axis=1)
    mean = df.groupby("id").mean().transpose()
    std_dev = df.groupby("id").std().transpose()

    t = np.arange(40)
    plt.setp(ax, xticks=np.arange(0, 41, 10.0))
    ax.plot(mean.Driver, label="Driver", color="#ff328b")
    ax.fill_between(t, mean.Driver-std_dev.Driver, \
            mean.Driver+std_dev.Driver, alpha=0.3,\
            color="#ff328b")

    ax.plot(mean.Wild, label="Wild Type", color="#7384ff")
    ax.fill_between(t, mean.Wild-std_dev.Wild, \
            mean.Wild+std_dev.Wild, alpha=0.3,\
            facecolor="#7384ff")

    if two:
        ax.plot(mean.r01 + mean.r10, label="(r,+)", color="#6e3100")
        lower_bound = mean.r01+ mean.r10-\
                (std_dev.r01**2+std_dev.r10**2)**(1/2.0)

        upper_bound = mean.r01+ mean.r10+\
                (std_dev.r01**2+std_dev.r10**2)**(1/2.0)

        ax.fill_between(t, upper_bound, lower_bound, alpha=0.3,\
                facecolor="#6e3100")

        ax.plot(mean.r11, label="Resistance", color="#008454")
        ax.fill_between(t, mean.r11-std_dev.r11,\
                mean.r11+std_dev.r11, \
                alpha=0.3, facecolor="#008454")
    else:
        ax.plot(mean.r01 + mean.r10, label="r01", color="#008454")
        lower_bound = mean.r01+ mean.r10-\
                (std_dev.r01**2+std_dev.r10**2)**(1/2.0)

        upper_bound = mean.r01+ mean.r10+\
                (std_dev.r01**2+std_dev.r10**2)**(1/2.0)

        ax.fill_between(t, lower_bound, upper_bound,\
                alpha=0.3, facecolor="#008454")

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 40)





def main(args):
    font = {'family' : 'sans-serif', 'size' : 11}
    plt.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rc('font', **font)

    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4)
    fig.set_size_inches(10, 7)

    axes_font_size = 11
    axes_title_size = 11

    location = 'Data/'
    if args.example:
        location = 'Data/Example/'
    set_ax(location + "High_resistance/X_1gRNA.txt", ax1, False)
    ax1.set_ylabel('Allele Frequency', fontsize=axes_font_size)
    ax1.set_title("One gRNA", fontsize=axes_title_size)

    set_ax(location + "High_resistance/X_2gRNA.txt", ax2)
    ax2.set_title("Two gRNA", fontsize=axes_title_size)

    set_ax(location + "High_resistance/Autosomal_1gRNA.txt", ax3, False)
    ax3.set_title("One gRNA", fontsize=axes_title_size)


    set_ax(location + "High_resistance/Autosomal_2gRNA.txt", ax4)
    ax4.set_title("Two gRNA", fontsize=axes_title_size)


    set_ax(location + "Medium_resistance/X_1gRNA.txt", ax5, False)
    ax5.set_ylabel('Allele Frequency', fontsize=axes_font_size)

    set_ax(location + "Medium_resistance/X_2gRNA.txt", ax6)

    set_ax(location + "Medium_resistance/Autosomal_1gRNA.txt", ax7, False)

    set_ax(location + "Medium_resistance/Autosomal_2gRNA.txt", ax8)
    legend_properties = {'weight':'bold', "size":10.5}


    set_ax(location + "Low_resistance/X_1gRNA.txt", ax9, False)
    ax9.set_ylabel('Allele Frequency', fontsize=axes_font_size)
    ax9.set_xlabel('Generation', fontsize=axes_font_size)

    set_ax(location + "Low_resistance/X_2gRNA.txt", ax10)
    ax10.set_xlabel('Generation', fontsize=axes_font_size)

    set_ax(location + "Low_resistance/Autosomal_1gRNA.txt", ax11, False)
    ax11.set_xlabel('Generation', fontsize=axes_font_size)

    set_ax(location + "Low_resistance/Autosomal_2gRNA.txt", ax12)
    ax12.set_xlabel('Generation', fontsize=axes_font_size)
    ax12.legend(loc="center right", prop=legend_properties, bbox_to_anchor=[1.02, 0.5])

    plt.tight_layout()
    fig.savefig('Figures/tiles.png', format='png', dpi=1500)




if __name__=="__main__":
    parser = ArgumentParser(description='Generate Figures.')
    parser.add_argument('-example', dest='example', action='store_true', default=False, help="Flag to generate figures from example data")
    args = parser.parse_args()
    main(args)
