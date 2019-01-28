import random

mt = [0 for i in range(10)]  # 微小管
wa = 1
wd = 0.1
wb = 1125
ws = 145
wf = 55
wh = 125
c = 0.5

def rule(ca, i):
	if c == 0:
		return
