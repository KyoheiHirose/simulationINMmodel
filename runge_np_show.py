#
#numpyをフル活用した結果解析
#
import re
import time
import numpy as np
import matplotlib
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from euler_high import dt, TIME


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
	print(cells[103])
	colors ={'S':'blue', 'G1e':'yellow', 'G2':'green', 'G1l':'orange', 'M':'deepskyblue', 'diff':'gray', 'div':'red'}
	x = np.arange(0,int(TIME/dt)+1)
	for i, cell in enumerate(cells):
		# print('cell',i,' = ', cell)
		t1=np.where(cells[i,:,3] != 0)[0].min()
		for t, phi in enumerate(cells[i,:,3]):
			if phi == 1:
				try:
					if cells[i,t+1,3] == 1:
						continue
					else:
						plt.plot(x[t1:t+2], cells[i,t1:t+2,2], color=colors['G2'])
						t1 = t+1
				except:
					plt.plot(x[t1:t+1], cells[i,t1:t+1,2], color=colors['G2'])
			elif phi == 2:
				try:
					if cells[i,t+1,3] == 3:
						continue
					else:
						plt.plot(x[t1:t+2], cells[i,t1:t+2,2], color=colors['M'])
						t1 = t+1
				except:
					plt.plot(x[t1:t+1], cells[i,t1:t+1,2], color=colors['M'])
			elif phi == 3:
				try:
					if cells[i,t+1,3] == 3:
						continue
					else:
						plt.plot(x[t1:t+2], cells[i,t1:t+2,2], color=colors['G1e'])
						t1 = t+1
				except:
					plt.plot(x[t1:t+1], cells[i,t1:t+1,2], color=colors['G1e'])
			elif phi == 4:
				try:
					if cells[i,t+1,3] == 4:
						continue
					else:
						plt.plot(x[t1:t+2], cells[i,t1:t+2,2], color=colors['G1l'])
						t1 = t+1
				except:
					plt.plot(x[t1:t+1], cells[i,t1:t+1,2], color=colors['G1l'])
			elif phi == 5:
				try:
					if cells[i,t+1,3] == 5:
						continue
					else:
						plt.plot(x[t1:t+2], cells[i,t1:t+2,2], color=colors['S'])
						t1 = t+1
				except:
					plt.plot(x[t1:t+1], cells[i,t1:t+1,2], color=colors['S'])
			elif phi == 6:
				try:
					if cells[i,t+1,3] == 6:
						continue
					else:
						plt.plot(x[t1:t+2], cells[i,t1:t+2,2], color=colors['div'])
						t1 = t+1
				except:
					plt.plot(x[t1:t+1], cells[i,t1:t+1,2], color=colors['div'])
			elif phi == 7:
				try:
					if cells[i,t+1,3] == 7:
						continue
					else:
						plt.plot(x[t1:t+1], cells[i,t1:t+1,2], color=colors['diff'])
						t1 = t+1
				except:
					plt.plot(x[t1:t+1], cells[i, t1:t+1,2], color=colors['diff'])
	plt.grid()
	plt.show()
	print('this program took ', time.time()-start, 'sec')