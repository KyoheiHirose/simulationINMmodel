# coding=utf-8
"""
誤差判定用の出力設定
テスト用にランダムな要因を削除
dt = 0.05  では小数点以下2位まで信頼できる
"""
import numpy as np
import matplotlib.pyplot as plt
import random
import math
from mpl_toolkits.mplot3d import Axes3D
import copy
import re
from prettytable import PrettyTable


# 定数
R = 10  # 細胞の直径
dt = 0.1  # 時間幅
N = 100  # 初期細胞の数
domX = 20  # 領域の設定
domY = 20
domZ = 100
beta = 15  # 方程式の定数の指定
gumma = 4
TIME = 100  # シミュレーションを行う時の長さ
V_INIT = 20.0  # 速度の初期値


class Cell(object):

    def __init__(self):
        self.phi = 0  # 細胞周期の状態
        self.r = np.array([0,0,0])  # 細胞の座標[0,dom)でランダム
        # self.v = np.array([0,0,0])  # 細胞の速度をV_INITの90~110%でランダムに与える
        self.timer = 0  # 細胞がその状態でどれだけ時間が経過したか
        self.timer2 = 0  # phi==4の時のみ用いる遅れを表現するタイマー
        self.dend = np.array([0,0,0])  # 神経突起の位置を定数として保存

    def calc_next(self, cell):
        """
        各細胞の状態に対応した振る舞いを実行させる関数
        """
        if self.phi == 1:  # G2
            self.phi1()
        elif self.phi == 2:  # M
            self.phi2()
        elif self.phi == 3:  # G1 突起ありの方
            self.phi3()
        elif self.phi == 4:  # G1 突起のない方
            self.phi4()
        elif self.phi == 5:  # S
            self.phi5()
        elif self.phi == 6:  # 細胞分裂
            self.phi6(cell)
        elif self.phi == 7:  # 領域外にでる
            self.phi7()

        self.border_change()

    def phi1(self):
        """G2期の細胞"""
        if self.r[2] < R/2:
            # print("change G2 -> M !!")
            self.phi = 2  # G2からMへの移行
            self.timer = 0.9  # M期の１時間を[0,1)のランダムで指定

    def phi2(self):
        """M期の細胞"""
        self.r[2] = R/2
        if self.timer > 0:
            self.timer -= dt
        else:
            # print("change M -> cell division !!")
            self.phi = 6  # 細胞分裂を行う状態に移行

    def phi3(self):
        """G1期(突起あり)の細胞"""
        if self.timer > 0:
            self.timer -= dt
        else:
            self.timer = 4
            self.phi = 5

    def phi4(self):
        """G1期(突起なし)の細胞"""
        if self.timer2 > 0:
            self.r[2] = R/2
            self.timer2 -= dt
        else:
            if self.timer > 0:
                self.timer -= dt
            else:
                self.timer = 4
                self.phi = 5

    def phi5(self):
        """S期"""
        if self.timer > 0:
            self.timer -= dt
        else:
            self.phi = 1  # G2期に移行

    def phi6(self, cell):
        theta = random.random() * 2 * math.pi  # 0~2piでランダムな角度生成
        rad = np.array([math.cos(theta), math.sin(theta), 0])
        # 現在iの細胞のコピーを生成し、cellsの末尾に追加
        divided_cell = copy.deepcopy(cell[i])
        cell = np.append(cell, [divided_cell], axis=0)
        # xy平面内での位置を変更して娘細胞とする
        cell[-1].phi = 4
        cell[-1].r -= R / 4 * rad
        cell[-1].v = cell[-1].v * np.array([-1, -1, 1])  # xy方向の速度を逆転
        cell[-1].dend = copy.deepcopy(cell[-1].r)  # 突起の位置を記録
        cell[-1].timer = 9
        cell[-1].timer2 = 1.5
        # 自分自身の座標を変更して娘細胞となる
        cell[i].phi = 3
        cell[i].r += R / 4 * rad
        cell[i].timer = 9
        # print("cell was divided")

    def phi7(self):
        if self.r[2] < 0:
            self.r[2] = R/2

    def border_change(self):
        if self.r[0] < 0:
            self.r[0] = domX - self.r[0]
        elif self.r[0] > domX:
            self.r[0] = self.r[0] - domX
        if self.r[1] < 0:
            self.r[1] = domY - self.r[1]
        elif self.r[1] > domY:
            self.r[1] = self.r[1] - domY

def f_int(cell, vect_r):
    """
    cells[i]に対して相互作用による力をすべて計算
    """
    f = np.zeros(3)
    for j in range(N):
        r_ij = cells[j].r - vect_r
        r_norm = np.linalg.norm(r_ij)
        e_r_ij = r_ij/r_norm
        if 0 < r_norm < R and i != j:  # 計算は細胞間距離がRより小さい場合のみ
            f += (-1)*beta * (R-r_norm) * e_r_ij
    if R/(-2) < vect_r[2] < R/2:
        r_norm = R/2 + vect_r[2]
        f += (-1)*beta * (R-r_norm)*np.array([0, 0, 1])
    return f


def V(r, phi, timer2):
    """timer2はphi==4用"""
    if phi == 1:
        v_phi = 8.5
    elif phi == 3:
        if r[2] < 15:
            v_phi = -1.6
        else:
            v_phi = 0
    elif phi == 4:
        if timer2 > 0:  # timer2>0の時は推進力は押さえ込まれてる
            v_phi = 0
        else:
            if r[2] < 15:
                v_phi = -1.6
            else:
                v_phi = 0
    else:
        v_phi = 0
    e_z = np.array([0, 0, 1])
    return v_phi * V_INIT * e_z


def H(r, dend):
    return gumma * (dend - r) * np.array([1, 1, 0])


def num_to_phase(phi):
    if phi == 1:
        return "G2"
    elif phi == 2:
        return "M"
    elif phi == 3:
        return "early G1"
    elif phi == 5:
        return "S"
    elif phi == 7:
        return "differentiated"
    else:
        return "late G1"


# main
# 初期値の読み込み
init = []
fin = open('inits_set', 'rt')
for line in fin:
    init += [re.split(',',line)]
fin.close()

cells = [Cell() for i in range(N)]

for i in range(N):
    r0 = float(init[i][0])
    r1 = float(init[i][1])
    r2 = float(init[i][2])
    cells[i].r = np.array([r0, r1, r2])
    # v0 = float(init[i][3])
    # v1 = float(init[i][4])
    # v2 = float(init[i][5])
    # cells[i].v = np.array([v0, v1, v2])
    phi = float(init[i][6])
    cells[i].phi = phi
    timer = float(init[i][7])
    cells[i].timer = timer
    timer2 = float(init[i][8])
    cells[i].timer2 = timer2
# fig = plt.figure()
# ax = Axes3D(fig)
TIME = 1  # シミュレーションを行う時の長さ
I = 10
time = [0]
dt = 0.1
t = dt
while t < TIME+0.0001:
    time += [t]
    t += dt

cells1 = copy.deepcopy(cells)
list1 = [cells1[I].r[2]]
dt = 0.02
t = dt
count = 1
countlist = [5*i for i in range(1,11)]
while t < TIME+0.001:
    for i in range(len(cells1)):
        # print("the size of cells1 is ", cells1.shape[0])
        # runge-kutta法によって次の時刻のr,vを取得
        ft = f_int(cells1, cells1[i].r) - V(cells1[i].r, cells1[i].phi, cells1[i].timer2) + H(cells1[i].r, cells1[i].dend)
        cells1[i].r += dt*ft

        # 細胞の状態を更新
        cells1[i].calc_next(cells1)
    if count in countlist:
        list1 += [cells1[I].r[2]]
    count += 1
    t += dt

    # list1 += [frunge]



cells2 = copy.deepcopy(cells)
list2 = [cells2[I].r[2]]

list2phi = [num_to_phase(cells2[I].phi)]
dt = 0.01
t = dt
count = 1
countlist = [10*i for i in range(1,11)]
while t < TIME+0.001:
    for i in range(len(cells2)):
        ft = f_int(cells2, cells2[i].r) - V(cells2[i].r, cells2[i].phi, cells2[i].timer2) + H(cells2[i].r, cells2[i].dend)
        cells2[i].r += dt*ft

        # 細胞の状態を更新
        cells2[i].calc_next(cells2)
    if count in countlist:
        list2 += [cells2[I].r[2]]
    t += dt
    count += 1
    # list1 += [frunge]
#
cells3 = copy.deepcopy(cells)
list3 = [cells3[I].r[2]]
dt = 0.005
t = dt
count = 1
countlist = [20*i for i in range(1,11)]
while t < TIME+0.001:
    for i in range(len(cells3)):
        ft = f_int(cells3, cells3[i].r) - V(cells3[i].r, cells3[i].phi, cells3[i].timer2) + H(cells3[i].r, cells3[i].dend)
        cells3[i].r += dt*ft
        cells3[i].calc_next(cells3)
    if count in countlist:
        list3 += [cells3[I].r[2]]
    t += dt
    count += 1

list4 = []
list5 = []
for i in range(len(list1)):
    list4 += [list1[i] - list2[i]]
    list5 += [list2[i] - list3[i]]

table = PrettyTable()
table.add_column('time', time)
table.add_column('euler dt = 0.02', list1)
table.add_column('euler dt = 0.01', list2)
table.add_column('euler dt = 0.005', list3)
table.add_column('dt=0.01の誤差', list4)
table.add_column('dt=0.005の誤差', list5)
print(table)
