#
#numpyをフル活用した結果解析
#
import re
import time
import numpy as np
import matplotlib
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from euler_high_init import dt, TIME


def cell_movement(posdyr):
    pass


if __name__ == '__main__':
    start = time.time()
    fin_position = open('results_posi_np', 'rt')
    fin_state = open('results_stat_np', 'rt')

    cells_list = []
    phi = []
    ind, col = 0, 0
    for i, (position, state) in enumerate(zip(fin_position, fin_state)):
        pos = re.split('/', position)
        pos.remove('\n')
        pos1 =[]
        for p in pos:
            p = re.split(',', p)
            pos1 += [p]
        pos = pos1
        sta = re.split(',', state)
        sta.remove(' \n')

        pos_list = []
        for j, (p, s) in enumerate(zip(pos, sta)):
            pos_list += [[float(p[0]), float(p[1]), float(p[2]), float(s)]]
            col = j
        cells_list.append(pos_list)
        ind = i
    fin_position.close()
    cells = np.zeros((col+1, ind+1,4))
    for i, cell in enumerate(cells_list):
        if i == 103:
            print(cell)
        for j, ce in enumerate(cell):
            cells[j,i] = ce

    print('実行時間TIME...', TIME, ', 時間幅dt...', dt, ', STEP数(TIME/dt)...', int(TIME/dt))
    print('細胞数...', len(cells))

    cells = np.array(cells) # cells[細胞番号,時刻,(0;x,1;y,2;z,3;phi)]

    print(len(np.where(cells[:,-1,-2]<100)[0]))