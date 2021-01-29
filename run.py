import sys, getopt, os

if __name__ == "__main__":
	file = sys.argv[1]
	stem = os.path.splitext(file)[0]
	try:
		opts, args = getopt.getopt(sys.argv[2:], "fcr", ["coeff" "colnames", "rownames"])
	except getopt.GetoptError:
		print('run.py <data-directory> --coeff --colnames --rownames')
		sys.exit(1)

	coeff = False
	colnames = False
	rownames = False
	for opt, arg in opts:
		if opt == '-h':
			print('run.py <data-directory> --coeff --colnames --rownames')
			sys.exit()
		elif opt in ("-f", "--coeff"):
			coeff = True
		elif opt in ("-c", "--colnames"):
			colnames = True
		elif opt in ("-r", "--rownames"):
			rownames = True

	pol = open('polaratio.cpp')
	content = pol.read()
	pol.close()
	data = open(file)
	if colnames: data.readline()
	n = len(data.readline().split()) - (1 if rownames else 0)
	m = len(data.read().split('\n'))
	pol = open('polaratio.cpp', 'w')
	pol.write('const int m = ' + str(m) + ', n = ' + str(n) + ';\n')
	pol.write(content.split('\n', 1)[1])
	pol.close()

	os.system('g++ -std=c++17 -o "polaratio" "polaratio".cpp -fopenmp')
	os.system('./polaratio ' + file + ' ' + stem + ' ' + str(1 if coeff else 0) + ' ' + str(1 if colnames else 0) + ' ' + str(1 if rownames else 0))