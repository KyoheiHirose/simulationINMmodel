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
import matplotlib
matplotlib.use('Agg') # -----(1)
from pylab import rcParams
rcParams['figure.figsize'] = 30,20

x = np.arange(0, int(TIME / dt) + 1)


def cell_movement(posdyr):
	for i, pdyr in enumerate(posdyr):
		if i in [i for i in range(posdyr.shape[0])]:
			# print('now i am pdyr',i, )
			t1=np.where(posdyr[i,:,3] != 0)[0].min()
			for t, phi in enumerate(posdyr[i,:,3]):
				if phi == 1:
					try:
						if posdyr[i,t+1,3] == 1:
							continue
						else:
							plt.plot(x[t1:t+2], posdyr[i,t1:t+2,2], color=colors['G2'])
							t1 = t+1
					except:
						plt.plot(x[t1:t+1], posdyr[i,t1:t+1,2], color=colors['G2'])
				elif phi == 2:
					try:
						if posdyr[i,t+1,3] == 3:
							continue
						else:
							plt.plot(x[t1:t+2], posdyr[i,t1:t+2,2], color=colors['M'])
							t1 = t+1
					except:
						plt.plot(x[t1:t+1], posdyr[i,t1:t+1,2], color=colors['M'])
				elif phi == 3:
					try:
						if posdyr[i,t+1,3] == 3:
							continue
						else:
							plt.plot(x[t1:t+2], posdyr[i,t1:t+2,2], color=colors['G1e'])
							t1 = t+1
					except:
						plt.plot(x[t1:t+1], posdyr[i,t1:t+1,2], color=colors['G1e'])
				elif phi == 4:
					try:
						if posdyr[i,t+1,3] == 4:
							continue
						else:
							plt.plot(x[t1:t+2], posdyr[i,t1:t+2,2], color=colors['G1l'])
							t1 = t+1
					except:
						plt.plot(x[t1:t+1], posdyr[i,t1:t+1,2], color=colors['G1l'])
				elif phi == 5:
					try:
						if posdyr[i,t+1,3] == 5:
							continue
						else:
							plt.plot(x[t1:t+2], posdyr[i,t1:t+2,2], color=colors['S'])
							t1 = t+1
					except:
						plt.plot(x[t1:t+1], posdyr[i,t1:t+1,2], color=colors['S'])
				elif phi == 6:
					try:
						if posdyr[i,t+1,3] == 6:
							continue
						else:
							plt.plot(x[t1:t+2], posdyr[i,t1:t+2,2], color=colors['div'])
							t1 = t+1
					except:
						plt.plot(x[t1:t+1], posdyr[i,t1:t+1,2], color=colors['div'])
				elif phi == 7:
					try:
						if posdyr[i,t+1,3] == 7:
							continue
						else:
							plt.plot(x[t1:t+1], posdyr[i,t1:t+1,2], color=colors['diff'])
							t1 = t+1
					except:
						plt.plot(x[t1:t+1], posdyr[i, t1:t+1,2], color=colors['diff'])
	plt.grid()
	return plt.show()


def cell_displacement(posdyr):
	pdyr = posdyr
	pdyr[:,:,3].dot(0)
	pass


def cell_displacement_1to34(cells):
	num6=0
	num0=0
	endphi3=0
	endphi4=0
	for i in range(cells.shape[0]):
		cell = cells[i]
		if np.where(cell[:,3]==6.)[0].size>0:  # 分裂した細胞のみを選別
			if np.where(cell[:,3]==6.)[0].size>0:
				num6 = np.where(cell[:,3]==6.)[0][0]  # 最初の分裂した時期を調べる
				ind = num6+1
				while 1:
					try:
						if cell[:,3][ind] == cell[:,3][ind+1]:
							ind += 1
						else:
							break
					except IndexError:
						break
				endphi3 = ind
				ind = num6-1
				startphi1 = 0
				endphi1 = 0
				while 1:
					if cell[:,3][ind] == 2.:
						ind -= 1
					elif cell[:,3][ind]==1.:
						try:
							if cell[:,3][ind+1]==2.:
								endphi1 = ind
								ind -= 1
							elif cell[:,3][ind-1]==5.:
								startphi1 = ind
								break
							elif ind == 0:
								break
							else:
								ind -= 1
						except IndexError:
							break
					else:
						print('something is wrong')
						break
			# plt.plot(cell[num6:endphi3,2])

		elif np.where(cell[:,3]==0.)[0].size>0:
			if np.where(cell[:,2]>5.)[0].size>0:
				num = np.where(cell[:,2]>5.)[0].min()
				ind = num
				while 1:
					try:
						if cell[:, 3][ind] == cell[:, 3][ind + 1]:
							ind += 1
						else:
							break
					except IndexError:
						break
				endphi4 = ind
				plt.plot(cell[num:endphi4,2])
			else:
				break
			#plt.plot(np.hstack((cell[startphi1:endphi1,2],cell[num6:endphi3,2])))
	plt.xlim(0,200)
	plt.ylim(0,35)
	# return plt.savefig('myo0_3.png') # -----(2)
	return plt.show()


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
		for j, ce in enumerate(cell):
			cells[j,i] = ce

	cells = np.array(cells) # cells[細胞番号,時刻,(0;x,1;y,2;z,3;phi)]

	print('実行時間TIME...', TIME, ', 時間幅dt...', dt, ', STEP数(TIME/dt)...', int(TIME/dt))
	print('全細胞数...', cells.shape[0])
	print('領域VZ内細胞数',np.where(cells[:,int(TIME/dt),2]<100)[0].shape[0])

	#############################経過をグラフ表示###############################
	colors ={'S':'deepskyblue', 'G1e':'magenta', 'G2':'green', 'G1l':'orange', 'M':'yellow', 'diff':'gray', 'div':'red'}

	cell_movement(cells)
	# cell_movement(cells[0,:,:3] - cells[0,0,:3])
	# for i, cell in enumerate(cells):
	# cellhist = np.where(cells[:,190:210,3]==6)
	# plt.hist(cellhist)
	# plt.show()
	print(np.where(cells[:,190:210,3]==6))

	# 各binの存在する細胞期の個数を計算

	# print(cells[:,100])
	# cell_displacement_1to34(cells)
	# plt.show()

	print('this program took ', time.time()-start, 'sec')