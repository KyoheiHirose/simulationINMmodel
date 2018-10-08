import re
import numpy as np
from euler_high import dt, TIME
import matplotlib.pyplot as plt


def dens_bins(pos, sta, step):
	d_count = 0
	bin0 = []
	bin1 = []
	bin2 = []
	bin3 = []
	bin4 = []
	bin5 = []
	bin6 = []
	bin7 = []
	bin8 = []
	bin9 = []
	bin10 = []
	while d_count < step:
		bin0 += [np.sum((0. < pos[:, d_count])&(pos[:,d_count] <= 10.))]
		bin1 += [np.sum((10. < pos[:, d_count])&(pos[:,d_count] <= 20.))]
		bin2 += [np.sum((20. < pos[:, d_count])&(pos[:,d_count] <= 30.))]
		bin3 += [np.sum((30. < pos[:, d_count])&(pos[:,d_count] <= 40.))]
		bin4 += [np.sum((40. < pos[:, d_count])&(pos[:,d_count] <= 50.))]
		bin5 += [np.sum((50. < pos[:, d_count])&(pos[:,d_count] <= 60.))]
		bin6 += [np.sum((60. < pos[:, d_count])&(pos[:,d_count] <= 70.))]
		bin7 += [np.sum((70. < pos[:, d_count])&(pos[:,d_count] <= 80.))]
		bin8 += [np.sum((80. < pos[:, d_count])&(pos[:,d_count] <= 90.))]
		bin9 += [np.sum((90. < pos[:, d_count])&(pos[:,d_count] <= 100.))]
		bin10 += [np.sum((100. < pos[:, d_count]))]
		d_count += 1
	x = np.arange(0, TIME, dt)
	plt.plot(x,bin0)
	plt.plot(x,bin1)
	plt.plot(x,bin2)
	plt.plot(x,bin3)
	plt.plot(x,bin4)
	plt.plot(x,bin5)
	plt.plot(x,bin6)
	plt.plot(x,bin7)
	plt.plot(x,bin8)
	plt.plot(x,bin9)
	plt.plot(x,bin10)
	return plt.show()


if __name__ == '__main__':
	STEP = int(TIME/dt)
	positions_list = []
	states_list = []
	fin_positions = open('results_posi', 'rt')
	fin_states = open('results_stat', 'rt')
	for pos, sta in zip(fin_positions, fin_states):
		pos = re.split(',', pos)
		sta = re.split(',', sta)
		pos.remove(' \n')
		sta.remove(' \n')
		num = len(pos)  # posとstaは同じサイズだから一方のみ
		pos_list = []
		sta_list = []
		for p, s in zip(pos, sta):
			pos_list += [float(p)]
			sta_list += [int(s)]
		positions_list += [pos_list]
		states_list += [sta_list]
	fin_positions.close()
	fin_states.close()

	num_cell = len(positions_list[0])
	# 値が0になっているところは細胞が存在しないことと同値とする
	# 行:細胞番号　列:step数
	positions = np.zeros((num_cell, STEP),float)
	states = np.zeros((num_cell, STEP), int)
	for i, (pos, sta) in enumerate(zip(positions_list, states_list)):
		for j, (p, s) in enumerate(zip(pos,sta)):
			positions[j,i] = p
			states[j,i] = s

	dens_bins(positions,states,STEP)
