# coding=utf-8
"""
クラスを用いずに配列のみで計算する
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import math
from mpl_toolkits.mplot3d import Axes3D
import copy

# constant
R = 10  # radius of cell
dt = 0.05
N = 100
DOM_X = 20
DOM_Y = 20
DOM_Z = 100
BETA = 15
GAMMA = 4
TIME = 20
V_INIT = 20.0

def calc_next(dyr):
    """

    :param dyr:
    :return:
    """
    if dyr