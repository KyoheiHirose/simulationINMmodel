# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import random
import math
from mpl_toolkits.mplot3d import Axes3D
import copy

# 定数
R = 10  # 細胞の直径
dt = 0.1  # 時間幅
N = 100  # 初期細胞の数
domX = 20
domY = 20
domZ = 100
beta = 15
gumma = 4
TIME = 100
V_INIT = 20.0
randNum = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99]


class Cell(object):

    def __init__(self, phi):
        self.phi = phi  # 細胞周期の状態
        self.r = np.random.rand(3) * np.array([domX-1, domY-1, domZ])  # 細胞の座標[0,dom)でランダム
        self.v = np.array([0., 0., V_INIT*random.uniform(0.9, 1.1)])  # 細胞の速度[0,1)でランダム Vを作用させるためにz方向のみ小さな値をあたえる
        self.timer = 0  # 細胞がその状態でどれだけ時間が経過したか
        self.timer2 = 0  # phi==4の時のみ用いる遅れを表現するタイマー
        self.dend = copy.deepcopy(self.r)  # 軸の位置

    def calc_next(self):
        """
        細胞を次の状態に移行させる
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
        elif self.phi == 7:  # 領域外にでる
            self.phi7()

    def phi1(self):
        """G2期の細胞"""
        if self.r[2] < R/2:
            print("change G2 -> M !!")
            if self.r[2] < 0:
                self.r[2] = R/2
            self.phi = 2  # G2からMへの移行
            self.timer = random.random()

    def phi2(self):
        """M期の細胞"""
        if self.r[2] < 0:
            self.r[2] = 0
        if self.timer > 0:
            self.timer -= dt
        else:
            # print("change M -> cell division !!")
            self.phi = 6  # 細胞分裂を行う状態に移行

    def phi3(self):
        """G1期(突起あり)の細胞"""
        if self.r[2] < 0:
            self.r[2] = 0
        if self.timer > 0:
            self.timer -= dt
        else:
            number = random.random()
            if number < 0.67:
                # print("change G1 -> S !!")
                self.phi = 5
                self.timer = 4  # S期は時間
            else:
                # print("change G1 -> I'm differentiated !!")
                self.phi = 7

    def phi4(self):
        """G1期(突起なし)の細胞"""
        if self.r[2] < 0:
            self.r[2] = 0
        if self.timer2 > 0:
            self.timer2 -= dt
        else:
            if self.timer > 0:
                self.timer -= dt
            else:
                number = random.random()
                if number < 0.67:
                    print("change G1 -> S !!")
                    self.phi = 5
                    self.timer = 4  # S期は時間
                else:
                    print("change G1 -> I'm differentiated !!")
                    self.phi = 7

    def phi5(self):
        if self.r[2] < 0:
            self.r[2] = 0
        if self.timer > 0:
            self.timer -= dt
        else:
            self.phi = 1  # G2期に移行

    def phi7(self):
        if self.r[2] < 0:
            self.r[2] = 0



def f_int(vect_r):
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
    if R/(-2) < vect_r[2] < R/2:
        r_norm = R/2 + vect_r[2]
        f += (-1)*beta * (R-r_norm)*np.array([0, 0, 0])

    return f


def V(r, v, phi, timer2):
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
                v_phi = -1.6  # 半分にした方が自然ではない????????????????????????????????
            else:
                v_phi = 0
    else:
        v_phi = 0
    e_z = np.array([0, 0, 1])
    return v_phi * v * e_z


def H(r, dend):
    return gumma * (dend - r) * np.array([1, 1, 0])


def numToPhase(phi):
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
cells = np.array([])  # 細胞を生成
for i in range(N):
    number = random.random()
    if number <= 0.1:
        cells = np.append(cells, Cell(1))
        cells[i].timer = 0
    elif number <= 5.5:
        cells = np.append(cells, Cell(3))
    else:
        cells = np.append(cells, Cell(4))

        cells[i].timer = 9
fig = plt.figure()
ax = Axes3D(fig)

counter = 0
t = 0
while t < TIME:
    print("t = ", t)
    print("=======================START============================")
    for i in range(cells.shape[0]):
        # print("the size of cells is ", cells.shape[0])
        # runnge-kutta法によって次の時刻のr,vを取得

        kb1 = f_int(cells[i].r) - V(cells[i].r, cells[i].v, cells[i].phi, cells[i].timer2) + H(cells[i].r, cells[i].dend)
        kb2 = f_int(cells[i].r+dt/2*kb1) - V(cells[i].r+dt/2*kb1, cells[i].v+dt/2*kb1, cells[i].phi, cells[i].timer2) + H(cells[i].r+dt/2*kb1, cells[i].dend)
        kb3 = f_int(cells[i].r+dt/2*kb2) - V(cells[i].r+dt/2*kb2, cells[i].v+dt/2*kb2, cells[i].phi, cells[i].timer2) + H(cells[i].r+dt/2*kb2, cells[i].dend)
        kb4 = f_int(cells[i].r+dt*kb3) - V(cells[i].r+dt*kb3, cells[i].v+dt*kb3, cells[i].phi, cells[i].timer2) + H(cells[i].r+dt*kb3, cells[i].dend)
        """
        if i in range(cells.shape[0]) and cells[i].phi == 4:
            print("I'm NO.", i, "My phase is ", numToPhase(cells[i].phi))
            print("  f_int = ", f_int(cells[i].r))
            print("  v = ", V(cells[i].r, cells[i].v, cells[i].phi, cells[i].timer2))
            print("  kb1 = ", kb1)
            print("  h = ", H(cells[i].r+dt/2*kb1, cells[i].dend))
            print("dr = ", dt*(kb1+2*kb2+2*kb3+kb4)/6)
        """
        cells[i].r += dt*(kb1+2*kb2+2*kb3+kb4)/6

        # 細胞の状態を更新
        cells[i].calc_next()

        # 細胞分裂によってあたらしいcellを生成
        if cells[i].phi == 6:
            thita = random.random() * 2 * math.pi  # 0~2piでランダムな角度生成
            rad = np.array([math.cos(thita), math.sin(thita), 0])
            # 現在iの細胞のコピーを生成し、cellsの末尾に追加
            dividedCell = copy.deepcopy(cells[i])
            cells = np.append(cells, [dividedCell], axis=0)
            # xy平面内での位置を変更して娘細胞とする
            cells[-1].phi = 4
            cells[-1].r -= R / 4 * rad
            cells[-1].v = cells[-1].v * np.array([-1, -1, 1])  # xy方向の速度を逆転
            cells[-1].dend = copy.deepcopy(cells[-1].r)  # 突起の位置を記録
            cells[-1].timer = 9 + random.uniform(-1, 1)
            cells[-1].timer2 = random.uniform(0, 3)
            # 自分自身の座標を変更して娘細胞となる
            cells[i].phi = 3
            cells[i].r += R / 4 * rad
            cells[i].timer = 9 + random.uniform(-1, 1)
            counter += 1
            # print("cell was divided")
        """
        if i in randNum:
            if cells[i].phi != 7:
                print("I'm No.", i, "My phase is ", numToPhase(cells[i].phi))
        """
    print("=======================END============================")

    # plot
    print("全細胞数...", cells.shape[0], "分裂回数...", counter)
    list1 = np.array([[0, 0, 0]])
    list2 = np.array([[0, 0, 0]])
    list3 = np.array([[0, 0, 0]])
    list4 = np.array([[0, 0, 0]])
    list5 = np.array([[0, 0, 0]])
    listother = np.array([[0, 0, 0]])
    # 代表の10個を表示
    for i in range(cells.shape[0]):  # range(cells.shape[0])で全部randNumで10個
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