import sys, getopt, os

if __name__ == "__main__":
	file = sys.argv[1]
	stem = os.path.splitext(file)[0]
	try:
		opts, args = getopt.getopt(sys.argv[2:], "b:cr", ["base=" "colnames", "rownames"])
	except getopt.GetoptError:
		print('run.py <data-directory> -b <base-for-coefficient> --colnames --rownames')
		sys.exit(1)

	base = 0
	colnames = False
	rownames = False
	for opt, arg in opts:
		if opt == '-h':
			print('run.py <data-directory> -b <base-for-coefficient> --colnames --rownames')
			sys.exit()
		elif opt in ("-b", "--base"):
			if arg <= 1:
				print('Base must be > 1')
				sys.exit()
			base = arg
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
	os.system('./polaratio ' + file + ' ' + stem + ' ' + str(base) + ' ' + str(1 if colnames else 0) + ' ' + str(1 if rownames else 0))