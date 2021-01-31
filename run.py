import sys, getopt, os

if __name__ == '__main__':
	file = sys.argv[1]
	stem = os.path.splitext(file)[0]
	try:
		opts, args = getopt.getopt(sys.argv[2:], 'fcrs', ['coeff', 'colnames', 'rownames', 'single'])
	except getopt.GetoptError:
		print('run.py <data-directory> --coeff --colnames --rownames --single')
		sys.exit(1)

	coeff = '0'
	colnames = '0'
	rownames = '0'
	single = '0'
	for opt, arg in opts:
		if opt == '-h':
			print('run.py <data-directory> --coeff --colnames --rownames --single')
			sys.exit()
		elif opt in ('-f', '--coeff'):
			coeff = '1'
		elif opt in ('-c', '--colnames'):
			colnames = '1'
		elif opt in ('-r', '--rownames'):
			rownames = '1'
		elif opt in ('-s', '--single'):
			single = '1'

	pol = open('polaratio.cpp')
	content = pol.read()
	pol.close()
	data = open(file)
	if colnames: data.readline()
	n = len(data.readline().split()) - (1 if rownames else 0)
	m = len(data.read().strip().split('\n')) + 1
	pol = open('polaratio.cpp', 'w')
	pol.write('const int m = ' + str(m) + ', n = ' + str(n) + ';\n')
	pol.write(content.split('\n', 1)[1])
	pol.close()

	os.system('g++ -std=c++17 -o polaratio polaratio.cpp' + ('' if single else ' -fopenmp'))
	os.system(('polaratio.exe ' if 'win' in sys.platform else './polaratio ') + file + ' ' + stem + ' ' + coeff + ' ' + colnames + ' ' + rownames)
