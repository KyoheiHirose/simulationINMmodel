import sys
import matplotlib.pyplot as plt
import random
import numpy as np

N = 10
MAXT = 50
P = 0.3

def rule(ca, i):
	"""
	次の時刻でのi番目の状態を返す
	11* -> 1
	01* -> 0
	10* -> *
	00* -> *
	apical 1** -> 01*
	VZ *11 -> pass
	"""
	if i == 0:
		p = random.random()
		if p < P:
			return 1
		else:
			return 0
	elif i == 1:
		if ca[0] == 1:
			return 1
		else:
			return ca[i]
	elif i == 2:
		if ca[1] == 1 and ca[0] == 1:
			return 1
		elif ca[1] == 1 and ca[0]==0:
			return 0
		else:
			return ca[i]
	elif ca[i-2] == 1 and ca[i-1] == 1:
		return 1
	elif ca[i-2] == 0 and ca[i-1] == 1:
		return 0
	elif ca[i-2] == 1 and ca[i-1] == 0:
		return ca[i]
	elif ca[i-2] == 0 and ca[i-1] == 0:
		return ca[i]
	else:
		print(ca[i-2],ca[i-1],ca[i])
		print('something wrong!!')


def next_state(ca):
	nextca = [0 for i in range(N)]
	for i in range(N):
		nextca[i] = rule(ca, i)
	for i in range(N):
		ca[i] = nextca[i]


if __name__=='__main__':
	outputdata = [[0 for i in range(N)] for j in range(MAXT+1)]
	ca = [int(random.random()*2) for i in range(N)]
	print(ca)
	for i in range(N):
		outputdata[0][i] = ca[i]
	for t in range(MAXT):
		next_state(ca)
		for i in range(N):
			outputdata[t+1][i]=ca[i]
	for i, output in enumerate(outputdata):
		output = output[::-1]
		outputdata[i] = output
	outputdata = np.array(outputdata)
	plt.imshow(outputdata.T,cmap='binary_r')
	plt.show()