
import re
import numpy as np
from euler_high_init import dt, TIME, phi3up, phi4up
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

t = np.arange(0, TIME, dt)  # 時間軸を用いる際の横軸
samples = [1,5,10,15,35,40,45,50,60,65,70,75]  # 左の番号のついたものだけをグラフ上に表示


def dens_bins(pos, sta, start, step):
	"""
	各binの細胞数を計算し出力する
	:param pos: positions
	:param sta: states
	:param start: 記録開始時間　これ以前の記録は表示されない
	:param step: STEP=TIME/dt
	:return: グラフの出力
	"""
	d_count = start
	bin0 = [[], [[],[],[],[],[],[],[]]]
	bin1 = [[], [[],[],[],[],[],[],[]]]
	bin2 = [[], [[],[],[],[],[],[],[]]]
	bin3 = [[], [[],[],[],[],[],[],[]]]
	bin4 = [[], [[],[],[],[],[],[],[]]]
	bin5 = [[], [[],[],[],[],[],[],[]]]
	bin6 = [[], [[],[],[],[],[],[],[]]]
	bin7 = [[], [[],[],[],[],[],[],[]]]
	bin8 = [[], [[],[],[],[],[],[],[]]]
	bin9 = [[], [[],[],[],[],[],[],[]]]
	while d_count < step:
		bin0[0] += [np.sum((0. < pos[:, d_count])&(pos[:,d_count] <= 10.))]
		bin1[0] += [np.sum((10. < pos[:, d_count])&(pos[:,d_count] <= 20.))]
		bin2[0] += [np.sum((20. < pos[:, d_count])&(pos[:,d_count] <= 30.))]
		bin3[0] += [np.sum((30. < pos[:, d_count])&(pos[:,d_count] <= 40.))]
		bin4[0] += [np.sum((40. < pos[:, d_count])&(pos[:,d_count] <= 50.))]
		bin5[0] += [np.sum((50. < pos[:, d_count])&(pos[:,d_count] <= 60.))]
		bin6[0] += [np.sum((60. < pos[:, d_count])&(pos[:,d_count] <= 70.))]
		bin7[0] += [np.sum((70. < pos[:, d_count])&(pos[:,d_count] <= 80.))]
		bin8[0] += [np.sum((80. < pos[:, d_count])&(pos[:,d_count] <= 90.))]
		bin9[0] += [np.sum((90. < pos[:, d_count])&(pos[:,d_count] <= 100.))]

		# phi == 1
		bin0[1][0] += [np.sum((0. < pos[:, d_count]) & (pos[:, d_count] <= 10.) & (sta[:, d_count] == 1))]
		bin1[1][0] += [np.sum((10. < pos[:, d_count]) & (pos[:, d_count] <= 20.) & (sta[:, d_count] == 1))]
		bin2[1][0] += [np.sum((20. < pos[:, d_count]) & (pos[:, d_count] <= 30.) & (sta[:, d_count] == 1))]
		bin3[1][0] += [np.sum((30. < pos[:, d_count]) & (pos[:, d_count] <= 40.) & (sta[:, d_count] == 1))]
		bin4[1][0] += [np.sum((40. < pos[:, d_count]) & (pos[:, d_count] <= 50.) & (sta[:, d_count] == 1))]
		bin5[1][0] += [np.sum((50. < pos[:, d_count]) & (pos[:, d_count] <= 60.) & (sta[:, d_count] == 1))]
		bin6[1][0] += [np.sum((60. < pos[:, d_count]) & (pos[:, d_count] <= 70.) & (sta[:, d_count] == 1))]
		bin7[1][0] += [np.sum((70. < pos[:, d_count]) & (pos[:, d_count] <= 80.) & (sta[:, d_count] == 1))]
		bin8[1][0] += [np.sum((80. < pos[:, d_count]) & (pos[:, d_count] <= 90.) & (sta[:, d_count] == 1))]
		bin9[1][0] += [np.sum((90. < pos[:, d_count]) & (pos[:, d_count] <= 100.) & (sta[:, d_count] == 1))]
		# phi == 2
		bin0[1][1] += [np.sum((0. < pos[:, d_count]) & (pos[:, d_count] <= 10.) & (sta[:, d_count] == 2))]
		bin1[1][1] += [np.sum((10. < pos[:, d_count]) & (pos[:, d_count] <= 20.) & (sta[:, d_count] == 2))]
		bin2[1][1] += [np.sum((20. < pos[:, d_count]) & (pos[:, d_count] <= 30.) & (sta[:, d_count] == 2))]
		bin3[1][1] += [np.sum((30. < pos[:, d_count]) & (pos[:, d_count] <= 40.) & (sta[:, d_count] == 2))]
		bin4[1][1] += [np.sum((40. < pos[:, d_count]) & (pos[:, d_count] <= 50.) & (sta[:, d_count] == 2))]
		bin5[1][1] += [np.sum((50. < pos[:, d_count]) & (pos[:, d_count] <= 60.) & (sta[:, d_count] == 2))]
		bin6[1][1] += [np.sum((60. < pos[:, d_count]) & (pos[:, d_count] <= 70.) & (sta[:, d_count] == 2))]
		bin7[1][1] += [np.sum((70. < pos[:, d_count]) & (pos[:, d_count] <= 80.) & (sta[:, d_count] == 2))]
		bin8[1][1] += [np.sum((80. < pos[:, d_count]) & (pos[:, d_count] <= 90.) & (sta[:, d_count] == 2))]
		bin9[1][1] += [np.sum((90. < pos[:, d_count]) & (pos[:, d_count] <= 100.) & (sta[:, d_count] == 2))]
		# phi == 3
		bin0[1][2] += [np.sum((0. < pos[:, d_count]) & (pos[:, d_count] <= 10.) & (sta[:, d_count] == 3))]
		bin1[1][2] += [np.sum((10. < pos[:, d_count]) & (pos[:, d_count] <= 20.) & (sta[:, d_count] == 3))]
		bin2[1][2] += [np.sum((20. < pos[:, d_count]) & (pos[:, d_count] <= 30.) & (sta[:, d_count] == 3))]
		bin3[1][2] += [np.sum((30. < pos[:, d_count]) & (pos[:, d_count] <= 40.) & (sta[:, d_count] == 3))]
		bin4[1][2] += [np.sum((40. < pos[:, d_count]) & (pos[:, d_count] <= 50.) & (sta[:, d_count] == 3))]
		bin5[1][2] += [np.sum((50. < pos[:, d_count]) & (pos[:, d_count] <= 60.) & (sta[:, d_count] == 3))]
		bin6[1][2] += [np.sum((60. < pos[:, d_count]) & (pos[:, d_count] <= 70.) & (sta[:, d_count] == 3))]
		bin7[1][2] += [np.sum((70. < pos[:, d_count]) & (pos[:, d_count] <= 80.) & (sta[:, d_count] == 3))]
		bin8[1][2] += [np.sum((80. < pos[:, d_count]) & (pos[:, d_count] <= 90.) & (sta[:, d_count] == 3))]
		bin9[1][2] += [np.sum((90. < pos[:, d_count]) & (pos[:, d_count] <= 100.) & (sta[:, d_count] == 3))]
		# phi == 4
		bin0[1][3] += [np.sum((0. < pos[:, d_count]) & (pos[:, d_count] <= 10.) & (sta[:, d_count] == 4))]
		bin1[1][3] += [np.sum((10. < pos[:, d_count]) & (pos[:, d_count] <= 20.) & (sta[:, d_count] == 4))]
		bin2[1][3] += [np.sum((20. < pos[:, d_count]) & (pos[:, d_count] <= 30.) & (sta[:, d_count] == 4))]
		bin3[1][3] += [np.sum((30. < pos[:, d_count]) & (pos[:, d_count] <= 40.) & (sta[:, d_count] == 4))]
		bin4[1][3] += [np.sum((40. < pos[:, d_count]) & (pos[:, d_count] <= 50.) & (sta[:, d_count] == 4))]
		bin5[1][3] += [np.sum((50. < pos[:, d_count]) & (pos[:, d_count] <= 60.) & (sta[:, d_count] == 4))]
		bin6[1][3] += [np.sum((60. < pos[:, d_count]) & (pos[:, d_count] <= 70.) & (sta[:, d_count] == 4))]
		bin7[1][3] += [np.sum((70. < pos[:, d_count]) & (pos[:, d_count] <= 80.) & (sta[:, d_count] == 4))]
		bin8[1][3] += [np.sum((80. < pos[:, d_count]) & (pos[:, d_count] <= 90.) & (sta[:, d_count] == 4))]
		bin9[1][3] += [np.sum((90. < pos[:, d_count]) & (pos[:, d_count] <= 100.) & (sta[:, d_count] == 4))]
		# phi == 5
		bin0[1][4] += [np.sum((0. < pos[:, d_count]) & (pos[:, d_count] <= 10.) & (sta[:, d_count] == 5))]
		bin1[1][4] += [np.sum((10. < pos[:, d_count]) & (pos[:, d_count] <= 20.) & (sta[:, d_count] == 5))]
		bin2[1][4] += [np.sum((20. < pos[:, d_count]) & (pos[:, d_count] <= 30.) & (sta[:, d_count] == 5))]
		bin3[1][4] += [np.sum((30. < pos[:, d_count]) & (pos[:, d_count] <= 40.) & (sta[:, d_count] == 5))]
		bin4[1][4] += [np.sum((40. < pos[:, d_count]) & (pos[:, d_count] <= 50.) & (sta[:, d_count] == 5))]
		bin5[1][4] += [np.sum((50. < pos[:, d_count]) & (pos[:, d_count] <= 60.) & (sta[:, d_count] == 5))]
		bin6[1][4] += [np.sum((60. < pos[:, d_count]) & (pos[:, d_count] <= 70.) & (sta[:, d_count] == 5))]
		bin7[1][4] += [np.sum((70. < pos[:, d_count]) & (pos[:, d_count] <= 80.) & (sta[:, d_count] == 5))]
		bin8[1][4] += [np.sum((80. < pos[:, d_count]) & (pos[:, d_count] <= 90.) & (sta[:, d_count] == 5))]
		bin9[1][4] += [np.sum((90. < pos[:, d_count]) & (pos[:, d_count] <= 100.) & (sta[:, d_count] == 5))]
		# phi == 6
		bin0[1][5] += [np.sum((0. < pos[:, d_count]) & (pos[:, d_count] <= 10.) & (sta[:, d_count] == 6))]
		bin1[1][5] += [np.sum((10. < pos[:, d_count]) & (pos[:, d_count] <= 20.) & (sta[:, d_count] == 6))]
		bin2[1][5] += [np.sum((20. < pos[:, d_count]) & (pos[:, d_count] <= 30.) & (sta[:, d_count] == 6))]
		bin3[1][5] += [np.sum((30. < pos[:, d_count]) & (pos[:, d_count] <= 40.) & (sta[:, d_count] == 6))]
		bin4[1][5] += [np.sum((40. < pos[:, d_count]) & (pos[:, d_count] <= 50.) & (sta[:, d_count] == 6))]
		bin5[1][5] += [np.sum((50. < pos[:, d_count]) & (pos[:, d_count] <= 60.) & (sta[:, d_count] == 6))]
		bin6[1][5] += [np.sum((60. < pos[:, d_count]) & (pos[:, d_count] <= 70.) & (sta[:, d_count] == 6))]
		bin7[1][5] += [np.sum((70. < pos[:, d_count]) & (pos[:, d_count] <= 80.) & (sta[:, d_count] == 6))]
		bin8[1][5] += [np.sum((80. < pos[:, d_count]) & (pos[:, d_count] <= 90.) & (sta[:, d_count] == 6))]
		bin9[1][5] += [np.sum((90. < pos[:, d_count]) & (pos[:, d_count] <= 100.) & (sta[:, d_count] == 6))]
		# phi ==27
		bin0[1][6] += [np.sum((0. < pos[:, d_count]) & (pos[:, d_count] <= 10.) & (sta[:, d_count] == 7))]
		bin1[1][6] += [np.sum((10. < pos[:, d_count]) & (pos[:, d_count] <= 20.) & (sta[:, d_count] == 7))]
		bin2[1][6] += [np.sum((20. < pos[:, d_count]) & (pos[:, d_count] <= 30.) & (sta[:, d_count] == 7))]
		bin3[1][6] += [np.sum((30. < pos[:, d_count]) & (pos[:, d_count] <= 40.) & (sta[:, d_count] == 7))]
		bin4[1][6] += [np.sum((40. < pos[:, d_count]) & (pos[:, d_count] <= 50.) & (sta[:, d_count] == 7))]
		bin5[1][6] += [np.sum((50. < pos[:, d_count]) & (pos[:, d_count] <= 60.) & (sta[:, d_count] == 7))]
		bin6[1][6] += [np.sum((60. < pos[:, d_count]) & (pos[:, d_count] <= 70.) & (sta[:, d_count] == 7))]
		bin7[1][6] += [np.sum((70. < pos[:, d_count]) & (pos[:, d_count] <= 80.) & (sta[:, d_count] == 7))]
		bin8[1][6] += [np.sum((80. < pos[:, d_count]) & (pos[:, d_count] <= 90.) & (sta[:, d_count] == 7))]
		bin9[1][6] += [np.sum((90. < pos[:, d_count]) & (pos[:, d_count] <= 100.) & (sta[:, d_count] == 7))]
		d_count += 1

	plt.plot(t, bin0[0], label='bin0', lw='3')
	plt.plot(t, bin1[0], label='bin1')
	plt.plot(t, bin2[0], label='bin2')
	plt.plot(t, bin3[0], label='bin3')
	plt.plot(t, bin4[0], label='bin4')
	plt.plot(t, bin5[0], label='bin5')
	plt.plot(t, bin6[0], label='bin6')
	plt.plot(t, bin7[0], label='bin7')
	plt.plot(t, bin8[0], label='bin8')
	plt.plot(t, bin9[0], label='bin9')
	plt.legend()
	plt.title('The density of each bins')
	plt.xlabel('time(hour)')

	# return plt.savefig('outputs/figbin_vphi=3')
	return plt.show()


# def g2tog1_movement(pos, sta):
# 	g2g1s_pos = pos[samples]
# 	g2g1s_sta = sta[samples]
# 	g2g1_ind = set(np.where(g2g1s_sta == 2)[0])
# 	for i,p in enumerate(pos):
# 		pass
# 		if i in samples:
		# 	plt.plot(t, p)
	# plt.title('CELL MOVEMENT FROM G2 to G1')
	# plt.ylabel('(um)')
	# plt.xlabel('time(hour)')
	#
	# return plt.show()


def cell_displacement(pos):
	"""
	各細胞の初期値からの変位を表示
	:param pos:
	:return:
	"""
	pos0 = np.reshape(pos[:, 0], (len(pos), 1))
	dplace = pos - pos0
	for i, dpla in enumerate(dplace):
		if i in samples:
			plt.plot(t, dpla, label='cell_' + str(i))
	plt.title('DISPLACEMENT')
	plt.ylabel('displacement(um)')
	plt.xlabel('time(hour)')
	plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
	# return plt.savefig('outputs/figdisp_init')
	return plt.show()


def cell_movement(pos):
	for i, p in enumerate(pos):
		if i in samples:
			plt.plot(t, p, label='cell_' + str(i))
	plt.title('MOVEMENT')
	plt.ylabel('(um)')
	plt.xlabel('time(hour)')
	plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
	# return plt.savefig('outputs/figmov_init')
	return plt.show()


def read_file(name, type):
	"""
	ファイルの読み込み　入力は行:時間　列:各細胞
	:param name: file name
	:param type: reading type; int or float in str
	:return name_list, num: num;合計細胞数
	"""
	name_list =[]
	num = 0  # 細胞数のカウンタ
	fin_name = open(name, 'rt')
	for i, name in enumerate(fin_name):
		nam = re.split(',', name)
		nam.remove(' \n')
		nam_list = []
		for n in nam:
			if type == 'int':
				nam_list += [int(n)]
			elif type == 'float':
				nam_list += [float(n)]
		name_list += [nam_list]
		num = i
	fin_name.close()
	return name_list, num



if __name__ == '__main__':
	STEP = int(TIME/dt)

	positions_list, num_cell = read_file('results_posi', 'float')
	states_list, num_cell = read_file('results_stat','int')
	# num_cell = len(positions_list[])
	print('STEP = ', STEP, '細胞数 = ',num_cell)
	# 値が0になっているところは細胞が存在しないことと同値とする
	# 行:細胞番号　列:step数
	positions = np.zeros((num_cell, STEP),float)
	states = np.zeros((num_cell, STEP), int)
	for i, (pos, sta) in enumerate(zip(positions_list, states_list)):
		for j, (p, s) in enumerate(zip(pos,sta)):
			positions[j,i] = p
			states[j,i] = s

	# dens_bins(positions,states, 0, STEP)
	# g2tog1_movement(positions, states)
	cell_displacement(positions)
	cell_movement(positions)