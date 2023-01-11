from pdb import pm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    draw_time_line()

def draw_time_line():
    pMap = pd.read_csv('./Results/303_2_30.csv')
    print(pMap)
    plt.rcParams["font.size"] = 17
    plt.scatter(pMap['SCP'],pMap['COP'], marker='o')

    plt.xlabel('SCP  W')
    plt.ylabel('COP a.u.')
    plt.savefig('./Fig/pMap.png')
    plt.show()


if __name__ == '__main__':
    main()