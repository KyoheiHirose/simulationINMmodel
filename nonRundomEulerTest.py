# coding=utf-8
"""
速度の初期条件を追加
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
dt = 0.05  # 時間幅
N = 100  # 初期細胞の数
domX = 20  # 領域の設定
domY = 20
domZ = 100
beta = 15  # 方程式の定数の指定
gumma = 4
TIME = 100  # シミュレーションを行う時の長さ
V_INIT = 20.0  # 速度の初期値
randNum = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99]  # 左の番号のついたものだけをグラフ上に表示


class Cell(object):

    def __init__(self):
        self.phi = 0  # 細胞周期の状態
        self.r = np.array([0,0,0])  # 細胞の座標[0,dom)でランダム
        self.v = np.array([0,0,0])  # 細胞の速度をV_INITの90~110%でランダムに与える
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

    def phi1(self):
        """G2期の細胞"""
        if self.r[2] <= R/2:
            # print("change G2 -> M !!")
            self.phi = 2  # G2からMへの移行
            self.timer = 0.9  # M期の１時間を[0,1)のランダムで指定

    def phi2(self):
        """
        M期の細胞
        apical面に固定
        """
        self.r[2] = R/2
        if self.timer > 0:
            self.timer -= dt
        else:
            # print("change M -> cell division !!")
            self.phi = 6  # 細胞分裂を行う状態に移行

    def phi3(self):
        """
        G1期(突起あり)の細胞
        """
        if self.timer > 0:
            self.timer -= dt
        else:
            self.timer = 4
            self.phi = 5

    def phi4(self):
        """
        G1期(突起なし)の細胞
        apical面に固定

        """
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
        """
        細胞分裂させる関数
        :param cell:
        :return:
        """
        theta = random.random() * 2 * math.pi  # 0~2piでランダムな角度生成
        rad = np.array([math.cos(theta), math.sin(theta), 0])
        # 現在iの細胞のコピーを生成し、cellsの末尾に追加
        divided_cell = copy.deepcopy(cell[i])
        cell = np.append(cell, [divided_cell], axis=0)
        # xy平面内での位置を変更して娘細胞とする
        cell[-1].phi = 4
        cell[-1].r -= R / 4 * rad
        cell[-1].dend = copy.deepcopy(cell[-1].r)  # 新たに分裂した細胞の突起の位置を記録　変更したい
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


def f_int(cell, vect_r):
    """
    cells[i]に対して相互作用による力をすべて計算
    """
    f = np.zeros(3)
    for j in range(N):
        r_ij = cell[j].r - vect_r
        r_norm = np.linalg.norm(r_ij)
        e_r_ij = r_ij/r_norm
        if 0 < r_norm < R and i != j:  # 計算は細胞間距離がRより小さい場合のみ
            f += (-1)*beta * (R-r_norm) * e_r_ij
    # if R/(-2) < vect_r[2] < R/2:
    #     r_norm = R/2 + vect_r[2]
    #     f += (-1)*beta * (R-r_norm)*np.array([0, 0, 1])
    # 境界条件...箱の中にあるとして境界に接すると反発を受ける係数は細胞と同じ
    if R/(-2) < vect_r[0] < R/2:
        r_norm = R/2 + vect_r[0]
        f += (-1)*beta * (R-r_norm)*np.array([-1, 0, 0])
    if domX-R/2 < vect_r[0] < domX+R/2:
        r_norm = domX + R/2 - vect_r[0]
        f += (-1)*beta * (R-r_norm)*np.array([1, 0, 0])
    if R/(-2) < vect_r[1] < R/2:
        r_norm = R/2 + vect_r[1]
        f += (-1)*beta * (R-r_norm)*np.array([0, -1, 0])
    if domY-R/2 < vect_r[1] < domY+R/2:
        r_norm = domY + R/2 - vect_r[1]
        f += (-1)*beta * (R-r_norm)*np.array([0, 1, 0])

    return f


def V(r, v, phi, timer2):
    """timer2はphi==4用"""
    if phi == 1:
        v_phi = 8.5
    elif phi == 3:
        if r[2] < 15:
            v_phi = -1.6
        else:
            v_phi = -0.5
    elif phi == 4:
        if timer2 > 0:  # timer2>0の時は推進力は押さえ込まれてる
            v_phi = 0
        else:
            if r[2] < 15:
                v_phi = -1.6
            else:
                v_phi = -0.5
    else:
        v_phi = 0
    e_z = np.array([0, 0, 1])
    return v_phi * v * e_z


def H(r, dend):
    return gumma * (dend - r) * np.array([1, 1, 0])


def num_to_phase(phi):
    if phi == 1:
        return "G2"
    elif phi == 2:
        return "M"
    elif phi == 3:
        return "early G1"
    elif phi == 4:
        return "late G1"
    elif phi == 5:
        return "S"
    elif phi == 7:
        return "differentiated"
    else:
        return "it's ERROR!!!"

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

fig = plt.figure()
ax = Axes3D(fig)

t = 0
while t < TIME+0.001:
    for i in range(len(cells)):
        ft = f_int(cells, cells[i].r) - V(cells[i].r, cells[i].v, cells[i].phi, cells[i].timer2) + H(cells[i].r, cells[i].dend)
        cells[i].r += dt*ft
        if cells[i].phi == 1:
            print("I'm G2 cell")
            print(ft[2])
        # 細胞の状態を更新
        cells[i].calc_next(cells)

    print("全細胞数...", len(cells), "分裂回数...")
    list1 = np.array([[0, 0, 0]])
    list2 = np.array([[0, 0, 0]])
    list3 = np.array([[0, 0, 0]])
    list4 = np.array([[0, 0, 0]])
    list5 = np.array([[0, 0, 0]])
    listother = np.array([[0, 0, 0]])
    # 代表の10個を表示
    for i in range(len(cells)):  # range(cells.shape[0])で全部randNumで10個
        if cells[i].phi == 2:  # M期
            list2 = np.append(list2, [cells[i].r], axis=0)
        elif cells[i].phi == 1:  # G2
            list1 = np.append(list1, [cells[i].r], axis=0)
        elif cells[i].phi == 3:  # early G1
            list3 = np.append(list3, [cells[i].r], axis=0)
        elif cells[i].phi == 4:  # late G1
            list4 = np.append(list4, [cells[i].r], axis=0)
        elif cells[i].phi == 5:  # S
            list5 = np.append(list5, [cells[i].r], axis=0)
        else:
            listother = np.append(listother, [cells[i].r], axis=0)
    # 各リストには初期化した分の要素が存在するから要素数はー1
    print("the number of G1 phase cells is ................", list1.shape[0]-1)
    print("the number of M phase cells is .................", list2.shape[0]-1)
    print("the number of early G1 phase cells is ..........", list3.shape[0]-1)
    print("the number of late G1 phase cells is ...........", list4.shape[0]-1)
    print("the number of S phase cells is .................", list5.shape[0]-1)
    print("the number of differentiated cells is ..........", listother.shape[0]-1)

    ax.cla()
    ax.plot(list1[1:, 0], list1[1:, 1], list1[1:, 2], "o", color="orange", ms=8, mew=0.5, label="G2")  # G2 phase
    ax.plot(list2[1:, 0], list2[1:, 1], list2[1:, 2], "o", color="red", ms=8, mew=0.5, label="M")  # M phase
    ax.plot(list3[1:, 0], list3[1:, 1], list3[1:, 2], "^", color="magenta", ms=8, mew=0.5, label="early G1")  # early G1
    ax.plot(list4[1:, 0], list4[1:, 1], list4[1:, 2], "^", color="pink", ms=8, mew=0.5, label="late G1")  # late G1
    ax.plot(list5[1:, 0], list5[1:, 1], list5[1:, 2], "o", color="deepskyblue", ms=8, mew=0.5, label="S")
    ax.plot(listother[1:, 0], listother[1:, 1], listother[1:, 2], "o", color="grey", ms=8, mew=0.5, label="differentiated")
    ax.set_xlim(0, domX)
    ax.set_ylim(0, domY)
    ax.set_zlim(0, domZ)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.pause(0.0001)
    for i in range(3):
        print(" ")
    t += dt

plt.show()