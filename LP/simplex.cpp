#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <vector>
using namespace std;
typedef long long LL;
typedef pair<int, int> PII;
#define MP(a, b) make_pair(a, b)
#define FOREACH(e,x) for(__typeof(x.begin()) e=x.begin();e!=x.end();++e)

const int maxn = 5, maxm = 5;
double INF = 1e30, eps = 1e-6;
int sig(double x) { return abs(x) < eps ? 0 : (x > 0 ? 1 : -1); }

struct Simplex {
	// 非基本变量, xj in N =0
	set<int> N;
	// 基本变量, xi in B = bi
	set<int> B;
	double b[maxn + maxm];
	double A[maxn + maxm][maxn + maxm];
	double c[maxn + maxm];
	double v;

	int n, m;

	double c_backup[maxn + maxm];
	int init(double _A[maxm][maxn], double _b[maxm], double _c[maxn], int _n, int _m) {
		N.clear(); B.clear();
		n = _n; m = _m;

		for (int j = 0; j < n; j++) N.insert(j);
		for (int i = n; i < n + m; i++) B.insert(i);
		for (int i = n; i < n + m; i++) b[i] = _b[i - n];
		memset(A, 0, sizeof(A));
		for (int i = n; i < n + m; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = _A[i - n][j];
		for (int j = 0; j < n; j++) c[j] = _c[j];
		v = 0;

		// find feasible solution
		int l = n;
		for (int i = n + 1; i < n + m; i++)
			if (b[i] < b[l]) l = i;
		if (sig(b[l]) < 0) {
			// add one dimension
			N.insert(n + m);
			memcpy(c_backup, c, sizeof(c));
			memset(c, 0, sizeof(c));
			c[n + m] = -1;
			for (int i = n; i < n + m; i++) A[i][n + m] = -1;
			pivot(n + m, l);
			double r = work();
			if (sig(r) != 0)
				return -1; // no solution
			
			// x[n + m] <= 0
			B.insert(n + m + 1);
			if (N.find(n + m) != N.end()) {
				A[n + m + 1][n + m] = 1.0;
			} else {
				FOREACH(j, N) A[n + m + 1][*j] = -A[n + m][*j];	
				b[n + m + 1] = -b[n + m];
			}
			// restore c and make it valid
			memcpy(c, c_backup, sizeof(c));
			for (int i = 0; i < n; i++)
			if (B.find(i) != B.end()) {
				// 非基本变量i 已经变成 基本变量了, 那么ci需change
				v += c[i] * b[i];
				FOREACH(j, N) c[*j] -= c[i] * A[i][*j];
				c[i] = 0.0;
			}
		}
		return 0;
	}

	void pivot(int e, int l) {
		static int que[maxn], qr;
		A[e][l] = 1.0 / A[l][e];
		qr = 0;
		// opt.3 A[l][*j] != 0, 623ms -> 46ms
		FOREACH(j, N) if (*j != e && sig(A[l][*j])) A[e][*j] = A[l][*j] / A[l][e], que[qr++] = *j;
		b[e] = b[l] / A[l][e];
		memset(A[l], 0, sizeof(A[l]));
		b[l] = 0;

		B.erase(l); B.insert(e);
		N.erase(e); N.insert(l);
		que[qr++] = l;
		// opt.2 A[i][e] != 0, 920ms -> 623ms
		FOREACH(_i, B) if (*_i != e && sig(A[*_i][e])) {
			int i = *_i;
			b[i] = b[i] - A[i][e] * b[e];
			for (int qi = 0; qi < qr; qi++) {
				int j = que[qi];
				A[i][j] -= A[i][e] * A[e][j];
			}
			A[i][e] = 0.0;
		}

		v += c[e] * b[e];
		FOREACH(j, N) c[*j] -= c[e] * A[e][*j];
		c[e] = 0.0;
	}

	double work() {
		while (1) {
			int e = -1; // 非基本变量 -> 基本变量
			double maxc = -INF;
			// opt.1 find max c, 951ms->920ms
			FOREACH(j, N) if (c[*j] > maxc) {
				maxc = c[*j];
				e = *j;
			}
			if (sig(maxc) <= 0)
				break;
			double delta = INF;
			int l; // 基本变量 -> 非基本变量
			FOREACH(_i, B) {
				int i = *_i;
				if (sig(A[i][e]) > 0) {
					double t_delta = b[i] / A[i][e];
					if (delta > t_delta) {
						delta = t_delta;
						l = i;
					}
				}
			}
			if (delta == INF) 
				return INF;
			pivot(e, l);
//FOREACH(i, N) printf("%d ", *i);
//puts("");
		}
		return v;
	}
} simplex;

double A[maxm][maxn], b[maxm], c[maxn];
int main() {
	// x1+x2<=20, 2*x1+x2<=25
	// max{2*x1 + 3*x2}
	A[0][0] = 1; A[0][1] = 1; b[0] = 20;
	A[1][0] = 2; A[1][1] = 1; b[1] = 25;
	// x2>=16
	A[2][0] = 0; A[2][1] = -1; b[2] = -16;
	c[0] = 3;	c[1] = 2;
	simplex.init(A, b, c, 2, 3);
	simplex.work();
	for (int j = 0; j < 2; j++) printf("%.2f\n", simplex.b[j]);
	printf("min = %.2f\n", simplex.v);

	return 0;
}
