#########################################
# ##高速化用
# 各細胞の位置変化と細胞周期をresultsファイルに保存
#########################################

import numpy as np
import re
import matplotlib.pyplot as plt
import random
import math
from mpl_toolkits.mplot3d import Axes3D
import copy
import calculate
import time

# CONSTANTS
R = 10
BETA = 15/6
GAMMA = 4
V_INIT = 20.0
N = 100
TIME = 5 # 20stepで0.25secほど 400stepで17sec
dt = 0.01
phi3up = 0
phi4up = 0
v_phi_ini3 = -1.6/6 + phi3up
v_phi_ini4 = -1.6/6 + phi4up


DOM_X = 20
DOM_Y = 20
DOM_Z = 100

class Cell(object):

	def __init__(self, x, y, z, xd, yd, phi, ct1, ct2, v):
		self.phi = phi
		# domain内にランダムに配置
		self.r = np.array([x,y,z])
		# z方向でV_INITに2割の幅をもたせて初期値を設定 <- 改良の余地あり
		self.v = np.array([0., 0., v])
		self.timer = ct1
		self.timer2 = ct2
		self.DEND = np.array([xd, yd, 0])

	def calc_next(self, cell, dt):
		"""
		細胞の状態に関するアクションを行う関数
		:param cell: 現在存在する細胞のリスト
		:return:
		"""
		if self.phi == 1:
			self.phi1()
		elif self.phi == 2:
			self.phi2(dt)
		elif self.phi == 3:
			self.phi3(dt)
		elif self.phi == 4:
			self.phi4(dt)  # 実験的にphi4を封印
		elif self.phi == 5:
			self.phi5(dt)
		elif self.phi == 6:
			self.phi6(cell)
		elif self.phi == 7:
			self.phi7()
		else:
			print('phi = ',self.phi,'---Error happen in calc_next')
		self.border_change()

	def phi1(self):
		""" G2 phase """
		if self.r[2] < R/2:  # apical面に到着したら、M phaseに入る
			self.phi = 2  # change to M phase
			self.timer = 1 + 2*random.random()  # G2-M phaseは[0,2)の乱数で表記

	def phi2(self, dt):
		""" M phase """
		self.r[2] = R/2  # M期の間はapical面に固定
		if self.timer > 0:
			self.timer -= dt
		else:
			self.phi = 6  # change to the cell division

	def phi3(self, dt):
		""" G1 phase (突起あり)"""
		if self.r[2] < R/2:
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

	def phi4(self, dt):
		""" G1 phase (突起なし) """
		if self.r[2] < R/2:
			self.r[2] = R/2
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

	def phi5(self, dt):
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
		if self.r[2] < 0:
			self.r[2] = 0

# @profile
# def f_intaract(cell, i):
#     """
#     cells[i]に対する相互作用による力を計算して和を取る
#     但し、phi == 2(M phase), 4(G1 phase 突起なし) は apical面(z == 0)にあるためz方向には拘束
#     :param cell:
#     :return: -BETA * f
#     """
#     f = [0, 0, 0]
#     for j in range(len(cell)):
#         r_ji = cell[j].r - cell[i].r
#         norm2 = r_ji[0]*r_ji[0] + r_ji[1]*r_ji[1] + r_ji[2]*r_ji[2]
#         if norm2 < R**2 and j != i:
#             norm = math.sqrt(norm2)
#             a = (R - norm) / norm
#             f[0] += a * r_ji[0]
#             f[1] += a * r_ji[1]
#             f[2] += a * r_ji[2]
#     if cell[i].phi == 2 or (cell[i].phi == 4 and cell[i].timer2 > 0):
#         f[2] *= 0
#     return f
#

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
			v_phi = v_phi_ini3
		else:
			v_phi = phi3up
	elif phi == 4:
		if timer2 > 0:
			v_phi = 0
		else:
			if r[2] < 10 + R/2:
				v_phi = v_phi_ini4
			else:
				v_phi = phi4up
	else:
		v_phi = 0
	return v_phi * v * np.array([0, 0, -1])


def h(r, DEND):
	return GAMMA *(1 - r[2]/DOM_Z)* (DEND - r) * np.array([1, 1, 0])


if __name__=="__main__":
	print('euler_high.py started!!')
	start = time.time()
	# 細胞の初期状態を作成
	f_position = open('final.dat', 'rt')  # final.dat...細胞体の位置x,y,z,軸の位置x,y,z,phi,timer,v
	cells_init = []
	for i, position in enumerate(f_position):
		pos = re.split(' ', position)
		for i, p in enumerate(pos):
			if i != 6:
				try:
					pos[i] = float(p)
				except:
					p = p[:-1]
					pos[i] = float(p)
			else:
				p = int(p)
				if p == 1 or p == 10:
					pos[i] = int(1)
				elif p == 3 or p == 30:
					pos[i] = int(3)
				elif p == 5:
					pos[i] = int(1)
				elif p == 6:
					pos[i] = int(5)
				elif p == 7 or p == 70:
					n = random.random()
					if n < 0.5:
						pos[i] = int(3)
					else:
						pos[i] = int(4)
				elif p == 8 or p == 9:
					pos[i] = int(7)
				else:
					pos[i] = int(p)
		cells_init += [pos]
	f_position.close()
	cells = []
	for i, cell in enumerate(cells_init):
		if cell[6] == 4:
			timer1 = 9*random.random()
			cells += [Cell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[6], timer1,cell[7], cell[8])]
		else:
			cells += [Cell(cell[0],cell[1],cell[2],cell[3],cell[4],cell[6],cell[7],0,cell[8])]

	t = 0
	fout_posi = open('results_posi_np_V', 'wt')
	fout_stat = open('results_stat_np_V', 'wt')
	while t < TIME:
		position = ''
		stats = ''
		# オイラー法で計算し、細胞の状態を更新
		for i in range(len(cells)):
			# euler method
			f = calculate.f_intaract(cells, i)
			ft = -1*BETA*np.array(f) \
				 + vphi(cells[i].r, cells[i].v, cells[i].phi, cells[i].timer2) \
				 + h(cells[i].r, cells[i].DEND)
			cells[i].r += ft * dt

			# 細胞の状態を更新
			cells[i].calc_next(cells, dt)

			# 細胞の状態を外部ファイルに出力
			position += str(cells[i].r[0]) + ',' + str(cells[i].r[1]) + ',' + str(cells[i].r[2]) + '/'
			stats += str(cells[i].phi) + ', '
		position += '\n'
		stats += '\n'
		fout_posi.write(position)
		fout_stat.write(stats)

		t += dt

	fout_posi.close()
	fout_stat.close()

	end = time.time() - start
	print('took ', end, 'sec')

