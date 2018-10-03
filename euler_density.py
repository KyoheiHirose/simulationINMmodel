# coding=utf-8
"""
euler法による数値解析
細胞の各領域内の密度
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
BETA = 15  # 相互作用に用いる定数
GAMMA = 4  # 平面内拘束に用いる
DELTA = 2  # apical面近傍で平面内拘束を強めるために用いる　packageと名前がかぶるらしいので注意
TIME = 50
V_INIT = 20.0


class Cell(object):

    def __init__(self, phi):
        self.phi = phi
        # domain内にランダムに配置
        self.r = np.random.rand(3) * np.array([DOM_X, DOM_Y, DOM_Z])  # 細胞の座標[0,dom)でランダム
        # z方向でV_INITに2割の幅をもたせて初期値を設定 <- 改良の余地あり
        self.v = np.array([0., 0., V_INIT])*random.uniform(0.9, 1.1) # 細胞の速度をV_INITの90~110%でランダムに与える
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
        if self.r[2] < R/2:  # apical面に到着したら、M phaseに入る
            self.phi = 2  # change to M phase
            self.timer = 1 + 2*random.random()  # G2-M phaseは[0,2)の乱数で表記

    def phi2(self):
        """ M phase """
        self.r[2] = R/2  # M期の間はapical面に固定
        if self.timer > 0:
            self.timer -= dt
        else:
            self.phi = 6  # change to the cell division

    def phi3(self):
        """ G1 phase (突起あり)"""
        if self.r[2] < 0:
            self.r[2] = R/2
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
            self.r[2] = R/2  # timer2がゼロになるまでapical面に固定
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
        """ cell division """
        theta = 2 * math.pi * random.random()  # 0~2piでランダムな分裂角度を生成
        rad = np.array([math.cos(theta), math.sin(theta), 0])
        cell += [copy.deepcopy(cell[i])]
        # in xy plane で位置を変更し突起なしの細胞を作る
        cell[-1].phi = 4
        cell[-1].r -= R/4 * rad
        cell[-1].DEND = copy.deepcopy(cell[-1].r)
        cell[-1].timer = 9 + random.uniform(-1, 1)
        cell[-1].timer2 = random.uniform(0, 3)  # apical面を離れない時間
        # 自分自身の座標を変更して娘細胞となる
        cell[i].phi = 3
        cell[i].r += R/4 * rad
        cell[i].timer = 9 + random.uniform(-1, 1)
        print('cell was divided!!!')

    def phi7(self):
        """ differentiated """
        # do nothing for now
        # if self.r[2] > DOM_Z + R/2:
        #     self.r[2] = DOM_Z + DOM_Z/2  # 150umで除外
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
            f += (R - norm)/norm * r_ji
    if cell[i].phi == 2 or (cell[i].phi == 4 and cell[i].timer2 > 0):
        f = f*np.array([1, 1, 0])
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
    if cell[i].r[0] < 0:
        cell[i].r[0] = 0
    elif DOM_X < cell[i].r[0]:
        cell[i].r[0] = DOM_X
    if cell[i].r[1] < 0:
        cell[i].r[1] = 0
    elif DOM_Y < cell[i].r[1]:
        cell[i].r[1] = DOM_Y


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
        if r[2] < 10 + R/2:
            v_phi = -1.6
        else:
            v_phi = -0.3
    elif phi == 4:
        if timer2 > 0:
            v_phi = 0
        else:
            if r[2] < 10 + R/2:
                v_phi = -1.6
            else:
                v_phi = -0.3
    else:
        v_phi = 0
    return v_phi * v * np.array([0, 0, -1])


def h(r, DEND):
    # delta = DELTA * abs((1-r[2]/DOM_Z))
    return GAMMA * (DEND - r) * np.array([1, 1, 0])


if __name__ == '__main__':

    # initiating
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

    # main loop

    t = 0
    while t < TIME:

        # listの初期化
        list0to10 = [0, 0, 0, 0, 0, 0, 0]  # bin [0, 10)の細胞をカウントphiとindexが一致
        list10to20 = [0, 0, 0, 0, 0, 0, 0]
        list20to30 = [0, 0, 0, 0, 0, 0, 0]
        list30to40 = [0, 0, 0, 0, 0, 0, 0]
        list40to50 = [0, 0, 0, 0, 0, 0, 0]
        list50to60 = [0, 0, 0, 0, 0, 0, 0]
        list60to70 = [0, 0, 0, 0, 0, 0, 0]
        list70to80 = [0, 0, 0, 0, 0, 0, 0]
        list80to90 = [0, 0, 0, 0, 0, 0, 0]
        list90to100 = [0, 0, 0, 0, 0, 0, 0]
        listover100 = [0, 0, 0, 0, 0, 0, 0]

        # 初期値を表示したいから更新前に情報を取り出す
        for i in range(len(cells)):
            if 0 < cells[i].r[2] <= 10:
                if cells[i].phi == 1:
                    list0to10[0] += 1
                elif cells[i].phi == 2:
                    list0to10[1] += 1
                elif cells[i].phi == 3:
                    list0to10[2] += 1
                elif cells[i].phi == 4:
                    list0to10[3] += 1
                elif cells[i].phi == 5:
                    list0to10[4] += 1
                elif cells[i].phi == 6:
                    list0to10[5] += 1
                else:
                    list0to10[6] += 1
            elif 10 < cells[i].r[2] <= 20:
                if cells[i].phi == 1:
                    list10to20[0] += 1
                elif cells[i].phi == 2:
                    list10to20[1] += 1
                elif cells[i].phi == 3:
                    list10to20[2] += 1
                elif cells[i].phi == 4:
                    list10to20[3] += 1
                elif cells[i].phi == 5:
                    list10to20[4] += 1
                elif cells[i].phi == 6:
                    list10to20[5] += 1
                else:
                    list20to30[6] += 1
            elif 20 < cells[i].r[2] <= 30:
                if cells[i].phi == 1:
                    list20to30[0] += 1
                elif cells[i].phi == 2:
                    list20to30[1] += 1
                elif cells[i].phi == 3:
                    list20to30[2] += 1
                elif cells[i].phi == 4:
                    list20to30[3] += 1
                elif cells[i].phi == 5:
                    list20to30[4] += 1
                elif cells[i].phi == 6:
                    list20to30[5] += 1
                else:
                    list20to30[6] += 1
            elif 30 < cells[i].r[2] <= 40:
                if cells[i].phi == 1:
                    list30to40[0] += 1
                elif cells[i].phi == 2:
                    list30to40[1] += 1
                elif cells[i].phi == 3:
                    list30to40[2] += 1
                elif cells[i].phi == 4:
                    list30to40[3] += 1
                elif cells[i].phi == 5:
                    list30to40[4] += 1
                elif cells[i].phi == 6:
                    list30to40[5] += 1
                else:
                    list30to40[6] += 1
            elif 40 < cells[i].r[2] <= 50:
                if cells[i].phi == 1:
                    list40to50[0] += 1
                elif cells[i].phi == 2:
                    list40to50[1] += 1
                elif cells[i].phi == 3:
                    list40to50[2] += 1
                elif cells[i].phi == 4:
                    list40to50[3] += 1
                elif cells[i].phi == 5:
                    list40to50[4] += 1
                elif cells[i].phi == 6:
                    list40to50[5] += 1
                else:
                    list40to50[6] += 1
            elif 50 < cells[i].r[2] <= 60:
                if cells[i].phi == 1:
                    list50to60[0] += 1
                elif cells[i].phi == 2:
                    list50to60[1] += 1
                elif cells[i].phi == 3:
                    list50to60[2] += 1
                elif cells[i].phi == 4:
                    list50to60[3] += 1
                elif cells[i].phi == 5:
                    list50to60[4] += 1
                elif cells[i].phi == 6:
                    list50to60[5] += 1
                else:
                    list50to60[6] += 1
            elif 60 < cells[i].r[2] <= 70:
                if cells[i].phi == 1:
                    list60to70[0] += 1
                elif cells[i].phi == 2:
                    list60to70[1] += 1
                elif cells[i].phi == 3:
                    list60to70[2] += 1
                elif cells[i].phi == 4:
                    list60to70[3] += 1
                elif cells[i].phi == 5:
                    list60to70[4] += 1
                elif cells[i].phi == 6:
                    list60to70[5] += 1
                else:
                    list60to70[6] += 1
            elif 70 < cells[i].r[2] <= 80:
                if cells[i].phi == 1:
                    list70to80[0] += 1
                elif cells[i].phi == 2:
                    list70to80[1] += 1
                elif cells[i].phi == 3:
                    list70to80[2] += 1
                elif cells[i].phi == 4:
                    list70to80[3] += 1
                elif cells[i].phi == 5:
                    list70to80[4] += 1
                elif cells[i].phi == 6:
                    list70to80[5] += 1
                else:
                    list70to80[6] += 1
            elif 80 < cells[i].r[2] <= 90:
                if cells[i].phi == 1:
                    list80to90[0] += 1
                elif cells[i].phi == 2:
                    list80to90[1] += 1
                elif cells[i].phi == 3:
                    list80to90[2] += 1
                elif cells[i].phi == 4:
                    list80to90[3] += 1
                elif cells[i].phi == 5:
                    list80to90[4] += 1
                elif cells[i].phi == 6:
                    list80to90[5] += 1
                else:
                    list80to90[6] += 1
            elif 90 < cells[i].r[2] <= 100:
                if cells[i].phi == 1:
                    list90to100[0] += 1
                elif cells[i].phi == 2:
                    list90to100[1] += 1
                elif cells[i].phi == 3:
                    list90to100[2] += 1
                elif cells[i].phi == 4:
                    list90to100[3] += 1
                elif cells[i].phi == 5:
                    list90to100[4] += 1
                elif cells[i].phi == 6:
                    list90to100[5] += 1
                else:
                    list90to100[6] += 1
            else:
                if cells[i].phi == 1:
                    listover100[0] += 1
                elif cells[i].phi == 2:
                    listover100[1] += 1
                elif cells[i].phi == 3:
                    listover100[2] += 1
                elif cells[i].phi == 4:
                    listover100[3] += 1
                elif cells[i].phi == 5:
                    listover100[4] += 1
                elif cells[i].phi == 6:
                    listover100[5] += 1
                else:
                    listover100[6] += 1
            # euler method
            ft = f_intaract(cells) \
                 + vphi(cells[i].r, cells[i].v, cells[i].phi, cells[i].timer2) \
                 + h(cells[i].r, cells[i].DEND)
            cells[i].r += ft*dt
            # 細胞の状態を更新
            cells[i].calc_next(cells)

        print('=========================t=', t,'=================================','\n','\n')
        t += dt


