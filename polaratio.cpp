const int m = 6323, n = 49;

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

double polar[n][n],mat[n][m], sum[n];
int invRank[n][m], ranks[n][m];
double one[m + 1], x[m + 1], y[m + 1], xy[m + 1];

int main(int argc, char **argv) {
	string dir = argv[1], stem = argv[2];
	ifstream in(dir);
	ios_base::sync_with_stdio(0); cin.tie(0);

	double zero = 1e9;
	bool coeff = stoi(argv[3]), colnames = stoi(argv[4]), rownames = stoi(argv[5]);
	if (colnames) { string line; getline(in, line); }
	for (int i = 0; i < m; i++) {
		string line; getline(in, line);
		istringstream iss(line);
		vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
		for (int j = 0; j < n; j++) {
			mat[j][i] = stod(tokens[j + rownames]);
			sum[j] += mat[j][i];
			zero = min(zero, mat[j][i]);
		}
	}
	for (int i = 0; i < n; i++) {
		iota(invRank[i], invRank[i] + m, 0);
		sort(invRank[i], invRank[i] + m, [=](int a, int b) { return mat[i][a] < mat[i][b]; });
		for (int j = 0; j < m; j++) {
			ranks[i][invRank[i][j]] = j;
		}
	}

	int i, j, num = 0, nC2 = n * (n-1) / 20, dec = 1;
	#pragma omp parallel private(i, j, one, x, y, xy)
	#pragma omp for
	for (j = 1; j < n; j++) {
		double *c2 = mat[j], avg2 = sum[j] / m;
		for (i = 0; i < j; i++) {
			double *c1 = mat[i], avg1 = sum[i] / m, cov = 0;
			for (int k = 0; k < m; k++) {
				cov += (avg1 - c1[k]) * (avg2 - c2[k]);
			}
			
			memset(one, 0, sizeof one);
			memset(x, 0, sizeof x);
			memset(y, 0, sizeof y);
			memset(xy, 0, sizeof xy);
			double neg = 0;

			for (int k : invRank[i]) {
				double ik = c1[k], jk = c2[k], ijk = ik * jk;
				if (ik > zero) {
					double t1 = 0, t2 = 0, t3 = 0, t4 = 0;
					int idx = m - ranks[j][k];
					while (idx) {
						t1 += one[idx];
						t2 += y[idx];
						t3 += x[idx];
						t4 += xy[idx];
						idx -= idx & -idx;
					}
					neg += -t1*ijk + t2*ik + t3*jk - t4;
				}

				if (jk > zero) {
					int idx = m - ranks[j][k];
					while (idx <= m) {
						one[idx]++;
						x[idx] += ik;
						y[idx] += jk;
						xy[idx] += ijk;
						idx += idx & -idx;
					}
				}
			}

			polar[i][j] = polar[j][i] = neg / (m * cov + neg + 1e-9);
			if (++num == nC2 * dec) cout << dec++ << "0%" << endl;
		}
	}

	ofstream out(stem + "PD.txt");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			out << polar[i][j] << " \n"[j == n - 1];
		}
	}
	if (coeff) {
		ofstream out(stem + "PC.txt");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				out << 2 / (polar[i][j] + 1) - 1 << " \n"[j == n - 1];
			}
		}
	}

	return 0;
}