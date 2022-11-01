from pdb import pm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    draw_time_line()

def draw_time_line():
    pMap = pd.read_csv('./Results/multi_solver.csv')
    print(pMap)
    plt.scatter(pMap['SCP'],pMap['COP'], marker='o')
    plt.rcParams["font.size"] = 14

    plt.xlabel('SCP  W')
    plt.ylabel('COP a.u.')
    plt.xlim(0,500)
    plt.ylim(0,4)
    plt.savefig('./Fig/pMap.png')
    plt.show()


if __name__ == '__main__':
    main()