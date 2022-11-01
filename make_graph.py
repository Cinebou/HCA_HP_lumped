import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    draw_time_line()

def draw_time_line():
    pMap = pd.read_csv('./Results/time_len_solver.csv')
    plt.plot(pMap.iloc[:,1],pMap.iloc[:,2], marker='o')
    plt.show()


if __name__ == '__main__':
    main()