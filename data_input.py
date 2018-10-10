"""
各細胞ごとの変位をデータとして取得する
"""
import re

if __name__ == '__main__':
	# read file
	fposit_x = open('results_posi_x', 'rt')
	fposit_y = open('results_posi_y', 'rt')
	fposit_z = open('results_posi_z', 'rt')


	posit_list_x = []
	posit_list_y = []
	posit_list_z = []

	for positx, posity, positz in zip(fposit_x, fposit_y, fposit_z):

		posx = re.split(',', positx)
		posy = re.split(',', posity)
		posz = re.split(',', positz)

		posx.remove(' \n')
		posy.remove(' \n')
		posz.remove(' \n')

		pos_list_x = []
		pos_list_y = []
		pos_list_z = []

		for px, py, pz in zip(posx, posy, posz):

			pos_list_x += [float(px)]
			pos_list_y += [float(py)]
			pos_list_z += [float(pz)]

		posit_list_x += [pos_list_x]
		posit_list_y += [pos_list_y]
		posit_list_z += [pos_list_z]

	fposit_x.close()
	fposit_y.close()
	fposit_z.close()

