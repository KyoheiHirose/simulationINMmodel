# coding=utf-8
"""
euler法による数値解析

"""

import numpy as np
import matplotlib.pyplot as plt
import random
import math
from mpl_toolkits.mplot3d import Axes3D
import copy

# constant
R = 10  # radius of cell
dt = 0.01
N = 100
DOM_X = 20
DOM_Y = 20
DOM_Z = 100
BETA = 15
GAMMA = 4
TIME = 50
V_INIT = 20.0


class Cell(object):

    def __init__(self, phi):
        self.phi = phi
        # domain内にランダムに配置
        self.r = np.random.rand(3) * np.array([DOM_X, DOM_Y, DOM_Z])  # 細胞の座標[0,dom)でランダム
        # z方向でV_INITに2割の幅をもたせて初期値を設定 <- 改良の余地あり
        self.v = np.array([0., 0., V_INIT]) * random.uniform(0.9, 1.1)  # 細胞の速度をV_INITの90~110%でランダムに与える
        self.timer = 0
        self.timer2 = 0
        self.DEND = copy.deepcopy(self.r) * np.array([1, 1, 0])

    def calc_next(self, cell):
        """
        細胞の状態に関するアクションを行う関数
        :param cell: 現在存在する細胞のリスト
        :return:
        """
        if self.phi == 1:
            self.phi1()
        elif self.phi == 2:
            self.phi2()
        elif self.phi == 3:
            self.phi3()
        elif self.phi == 4:
            self.phi4()
        elif self.phi == 5:
            self.phi5()
        elif self.phi == 6:
            self.phi6(cell)
        elif self.phi == 7:
            self.phi7()
        else:
            print('---Error happen in calc_next')
        self.border_change()

    def phi1(self):
        """ G2 phase """
        if self.r[2] < R / 2:  # apical面に到着したら、M phaseに入る
            self.phi = 2  # change to M phase
            self.timer = 1 + 2 * random.random()  # G2-M phaseは[0,2)の乱数で表記

    def phi2(self):
        """ M phase """
        self.r[2] = R / 2  # M期の間はapical面に固定
        if self.timer > 0:
            self.timer -= dt
        else:
            self.phi = 6  # change to the cell division

    def phi3(self):
        """ G1 phase (突起あり)"""
        if self.r[2] < 0:
            self.r[2] = R / 2
        if self.timer > 0:
            self.timer -= dt
        else:
            num = random.random()
            if num < 0.67:  # 67%は S phase に再突入
                self.phi = 5
                self.timer = 4
            else:  # the others は differentiated に
                self.phi = 7

    def phi4(self):
        """ G1 phase (突起なし) """
        if self.timer2 > 0:
            self.r[2] = R / 2  # timer2がゼロになるまでapical面に固定
            self.timer2 -= dt
        else:
            if self.timer > 0:
                self.timer -= dt
            else:
                num = random.random()
                if num < 0.67:
                    self.phi = 5
                    self.timer = 4  # S期は時間
                else:
                    self.phi = 7

    def phi5(self):
        """ S phase """
        if self.timer > 0:
            self.timer -= dt
        else:
            self.phi = 1  # change to G2 phase

    def phi6(self, cell):
        theta = 2 * math.pi * random.random()  # 0~2piでランダムな分裂角度を生成
        rad = np.array([math.cos(theta), math.sin(theta), 0])
        cell += [copy.deepcopy(cell[i])]
        # in xy plane で位置を変更し突起なしの細胞を作る
        cell[-1].phi = 4
        cell[-1].r -= R / 4 * rad
        cell[-1].DEND = copy.deepcopy(cell[-1].r)
        cell[-1].timer = 9 + random.uniform(-1, 1)
        cell[-1].timer2 = random.uniform(0, 3)  # apical面を離れない時間
        # 自分自身の座標を変更して娘細胞となる
        cell[i].phi = 3
        cell[i].r += R / 4 * rad
        cell[i].timer = 9 + random.uniform(-1, 1)
        print('cell was divided!!!')

    def phi7(self):
        pass

    def border_change(self):
        if self.r[0] < 0:
            self.r[0] = DOM_X - self.r[0]
        elif self.r[0] > DOM_X:
            self.r[0] = self.r[0] - DOM_X
        if self.r[1] < 0:
            self.r[1] = DOM_Y - self.r[1]
        elif self.r[1] > DOM_Y:
            self.r[1] = self.r[1] - DOM_Y


def f_intaract(cell):
    """
    cells[i]に対する相互作用による力を計算して和を取る
    但し、phi == 2(M phase), 4(G1 phase 突起なし) は apical面(z == 0)にあるためz方向には拘束
    :param cell:
    :return: -BETA * f
    """
    f = np.zeros(3)
    for j in range(len(cell)):
        r_ji = cell[j].r - cell[i].r
        norm = np.linalg.norm(r_ji)
        if 0 < norm < R:
            f += (R - norm) / norm * r_ji
    if cell[i].phi == 2 or (cell[i].phi == 4 and cell[i].timer2 > 0):
        f = f * np.array([1, 1, 0])
    # 境界条件領域外に細胞
    # if R/(-2) < cell[i].r[0] < R/2:
    #     r_norm = R/2 + cell[i].r[0]
    #     f += (R-r_norm)*np.array([-1, 0, 0])
    # elif DOM_X-R/2 < cell[i].r[0] < DOM_X+R/2:
    #     r_norm = DOM_X + R/2 - cell[i].r[0]
    #     f += (R-r_norm)*np.array([1, 0, 0])
    # if R/(-2) < cell[i].r[1] < R/2:
    #     r_norm = R/2 + cell[i].r[1]
    #     f += (R-r_norm)*np.array([0, -1, 0])
    # elif DOM_Y-R/2 < cell[i].r[1] < DOM_Y+R/2:
    #     r_norm = DOM_Y + R/2 - cell[i].r[1]
    #     f += (R-r_norm)*np.array([0, 1, 0])
    # 境界条件2領域外R/2に壁
    # if cell[i].r[0] < 0:
    #     cell[i].r[0] = 0
    # elif DOM_X < cell[i].r[0]:
    #     cell[i].r[0] = DOM_X
    # if cell[i].r[1] < 0:
    #     cell[i].r[1] = 0
    # elif DOM_Y < cell[i].r[1]:
    #     cell[i].r[1] = DOM_Y

    return -1 * BETA * f


def vphi(r, v, phi, timer2):
    """
    :param r: 現在の細胞の位置
    :param v: 細胞の速度
    :param phi:
    :param timer2:
    :return: v_phi * v * np.array([0, 0, 1])
    """
    if phi == 1:
        v_phi = 8.5
    elif phi == 3:
        if r[2] < 10 + R / 2:
            v_phi = -1.6
        else:
            v_phi = -0.1
    elif phi == 4:
        if timer2 > 0:
            v_phi = 0
        else:
            if r[2] < 10 + R / 2:
                v_phi = -1.6
            else:
                v_phi = 0
    else:
        v_phi = 0
    return v_phi * v * np.array([0, 0, -1])


def h(r, DEND):
    return GAMMA * (DEND - r) * np.array([1, 1, 0])


if __name__ == '__main__':
    cells = []
    for i in range(N):
        number = random.random()
        if number <= 0.1:
            cells += [Cell(1)]
            cells[i].timer = 0
        elif number <= 0.55:
            cells += [Cell(3)]
            z = cells[i].r[2]
            cells[i].timer = 9 * (DOM_Z - z) / DOM_Z
        else:
            cells += [Cell(4)]
            z = cells[i].r[2]
            cells[i].timer = 9 * (DOM_Z - z) / DOM_Z

    fig = plt.figure()
    ax = Axes3D(fig)
    t = 0
    while t < TIME:
        for i in range(len(cells)):
            # euler method
            ft = f_intaract(cells) \
                 + vphi(cells[i].r, cells[i].v, cells[i].phi, cells[i].timer2) \
                 + h(cells[i].r, cells[i].DEND)
            cells[i].r += ft * dt
            # 細胞の状態を更新
            cells[i].calc_next(cells)
        # 以下plotに関する
        list1 = []
        list2 = []
        list3 = []
        list4 = []
        list5 = []
        list6 = []
        list7 = []

        for i in range(len(cells)):
            if cells[i].phi == 1:
                list1 += [cells[i].r]
            elif cells[i].phi == 2:
                list2 += [cells[i].r]
            elif cells[i].phi == 3:
                list3 += [cells[i].r]
            elif cells[i].phi == 4:
                list4 += [cells[i].r]
            elif cells[i].phi == 5:
                list5 += [cells[i].r]
            elif cells[i].phi == 6:
                list6 += [cells[i].r]
            elif cells[i].phi == 7:
                list7 += [cells[i].r]
            else:
                print('Unexpected error')

        list1 = np.array(list1)
        list2 = np.array(list2)
        list3 = np.array(list3)
        list4 = np.array(list4)
        list5 = np.array(list5)
        list6 = np.array(list6)
        list7 = np.array(list7)
        print('=========================t=', t, '=================================', '\n', '\n')
        # print('the number of cells is .....................................', len(cells))
        # print("the number of G1 phase cells is ............................", list1.shape[0])
        # print("the number of M phase cells is .............................", list2.shape[0])
        # print("the number of early G1 phase cells is ......................", list3.shape[0])
        # print("the number of late G1 phase cells is .......................", list4.shape[0])
        # print("the number of S phase cells is .............................", list5.shape[0])
        # print("the number of cell division phase cells is .................", list6.shape[0])
        # print("the number of differentiated cells is ......................", list7.shape[0])
        t += dt

        ax.cla()
        if list1.shape[0] > 0:
            ax.plot(list1[:, 0], list1[:, 1], list1[:, 2], "o", color = "orange", ms = 8, mew = 0.5,
                    label = "G2")  # G2 phase
        if list2.shape[0] > 0:
            ax.plot(list2[:, 0], list2[:, 1], list2[:, 2], "o", color = "red", ms = 8, mew = 0.5,
                    label = "M")  # M phase
        if list3.shape[0] > 0:
            ax.plot(list3[:, 0], list3[:, 1], list3[:, 2], "^", color = "magenta", ms = 8, mew = 0.5,
                    label = "early G1")  # early G1
        if list4.shape[0] > 0:
            ax.plot(list4[:, 0], list4[:, 1], list4[:, 2], "^", color = "pink", ms = 8, mew = 0.5,
                    label = "late G1")  # late G1
        if list5.shape[0] > 0:
            ax.plot(list5[:, 0], list5[:, 1], list5[:, 2], "o", color = "deepskyblue", ms = 8, mew = 0.5, label = "S")
        # if list6.shape[0] > 0:
        #     ax.plot(list6[:, 0], list6[:, 1], list6[:, 2], "o", color="red", ms=8, mew=0.5, label="S")
        if list7.shape[0] > 0:
            ax.plot(list7[:, 0], list7[:, 1], list7[:, 2], "o", color = "grey", ms = 8, mew = 0.5,
                    label = "differentiated")
        ax.set_xlim(0, DOM_X)
        ax.set_ylim(0, DOM_Y)
        ax.set_zlim(0, DOM_Z + 100)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.legend()
        plt.pause(0.001)
        for i in range(3):
            print(" ")

    plt.show()
