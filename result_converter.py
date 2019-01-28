import time
import re
import numpy as np
import random
from euler_high_init import dt, TIME

if __name__ == '__main__':

	start = time.time()
	f_position = open('final.dat', 'rt')
	#final.dat...細胞体の位置x,y,z,軸の位置x,y,z,phi,timer,v
	cells = []
	ind, col = 0, 0
	for i, position in enumerate(f_position):
		print(position)
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
					if n <0.5:
						pos[i] = int(3)
					else:
						pos[i] = int(4)
				elif p == 8 or p == 9:
					pos[i] = int(7)
		cells += [pos]
	f_position.close()
	print(len(cells))
