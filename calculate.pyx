import math
cdef unsigned int R = 10
def f_intaract(cell, unsigned int i):
    """
    cells[i]に対する相互作用による力を計算して和を取る
    但し、phi == 2(M phase), 4(G1 phase 突起なし) は apical面(z == 0)にあるためz方向には拘束
    :param cell:
    :return: -BETA * f
    """
    cdef float norm, norm2, a
    f = [0., 0., 0.]
    for j in range(len(cell)):
        r_ji = cell[j].r - cell[i].r
        norm2 = r_ji[0]*r_ji[0] + r_ji[1]*r_ji[1] + r_ji[2]*r_ji[2]
        if norm2 < R**2 and j != i:
            norm = math.sqrt(norm2)
            a = (R - norm) / norm
            f[0] += a * r_ji[0]
            f[1] += a * r_ji[1]
            f[2] += a * r_ji[2]
    if cell[i].phi == 2 or (cell[i].phi == 4 and cell[i].timer2 > 0):
        f[2] *= 0
    return f