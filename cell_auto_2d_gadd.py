import sys
import matplotlib.pyplot as plt
import random
import numpy as np

X = 20
Y = 20
TMAX = 10

def rule(ca,i,j):
	"""
	次の時刻でのi番目の状態を返す
	1/8で分化してしまう
	"""
	p = random.random()
	if p < (ca[i+1][j+1]+ca[i][j+1]+ca[i-1][j+1]+ca[i+1][j]+ca[i-1][j]+ca[i+1][j-1]+ca[i][j-1]+ca[i-1][j-1])/8:
		print(i,',',j,' is ','differentiated!!')
		return 1
	else:
		if p < 0.2:
			return 1
		else:
			return 0


def next_state(ca):
	nextca = np.array([[0 for i in range(X+2)] for j in range(Y+2)])
	for i in range(1,X):
		for j in range(1,Y):
			new = rule(ca, i, j)
			print(new)
			nextca[i][j]=new
	print('nextca')
	print(nextca)
	for i in range(X+2):
		for j in range(Y+2):
			ca[i][j] = nextca[i][j]
	print('ca')
	print(ca)


if __name__=='__main__':
	#caの境界はすべて0で固定
	ca = np.array([[int(random.random()>0.7) for i in range(X+2)] for j in range(Y+2)])
	ca[0] = 0
	ca[Y+1] = 0
	ca[:,0] = 0
	ca[:,X+1] = 0
	print('init state is ')
	print(ca)
	plt.imshow(ca[1:X, 1:Y], cmap='binary_r')
	title = 't = ' + str(0)
	plt.title(title)
	plt.show()

	for t in range(TMAX):
		print('t = ', t+1)
		next_state(ca)
		print(ca)
		plt.imshow(ca[1:X, 1:Y], cmap='binary_r')
		title = 't = ' + str(t+1)
		plt.title(title)
		plt.show()
		# for i in range(N):
		# 	outputdata[t+1][i]=ca[i]
	# for i, output in enumerate(outputdata):
	# 	output = output[::-1]
	# 	outputdata[i] = output
	# outputdata = np.array(outputdata)
