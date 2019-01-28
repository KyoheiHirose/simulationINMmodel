import sys
import matplotlib.pyplot as plt
import random
import numpy as np

X = 10
Y = 10
Z = 10
MAXT = 50
P = 0.3


def rule(ca, i,j,k):
	"""
	次の時刻でのi番目の状態を返す
	11* -> 1
	01* -> 0
	10* -> *
	00* -> *
	apical 1** -> 01*
	VZ *11 -> pass
	"""
	if k == 0:
		p = random.random()
		if p < P:
			return 1
		else:
			return 0
	elif k == 1:
		if ca[0] == 1:
			return 1
		else:
			if ca[i-1,j+1,k-1]+ca[i,j+1,k-1]+ca[i+1,j+1,k-1]+ca[i-1,j,k-1]+ca[i+1,j,k-1]+ca[i-1,j-1,k-1]+ca[i,j-1,k-1]+ca[i+1,j-1,k-1] < 5:
				return ca[i]
			else:
				0
	elif k > 1:
		if ca[i,j,k-1] == 1:
			if ca[i,j,k-2] == 1:
				return 1
			elif ca[i-1,j+1,k-2]+ca[i,j+1,k-2]+ca[i+1,j+1,k-2]+ca[i-1,j,k-2]+ca[i+1,j,k-2]+ca[i-1,j-1,k-2]+ca[i,j-1,k-2]+ca[i+1,j-1,k-2] > 4:
				return 1
			else:
				return 0
		else:
			if ca[i-1,j+1,k-1]+ca[i,j+1,k-1]+ca[i+1,j+1,k-1]+ca[i-1,j,k-1]+ca[i+1,j,k-1]+ca[i-1,j-1,k-1]+ca[i,j-1,k-1]+ca[i+1,j-1,k-1] < 5:
				return ca[i,j,k]
			else:
				return 0


def next_state(ca):
	nextca = [[0 for i in range(X)] for j in range(Y)]
	for k in range(Z):
		for j in range(Y):
			for i in range(X):
				nextca[i,j] = rule(ca, i, j, k)
	for i in range(N):
		ca[i] = nextca[i]


if __name__=='__main__':
	outputdata = [[0 for i in range(X)] for j in range(MAXT+1)]
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