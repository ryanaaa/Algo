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

	double x[maxn + maxm];
	int n, m;

	void init(double _A[maxm][maxn], double _b[maxm], double _c[maxn], int _n, int _m) {
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
		memset(x, 0, sizeof(x));
		//for (int i = n; i < n + m; i++) x[i] = b[i - n];
	}

	void pivot(int e, int l) {
		A[e][l] = 1.0 / A[l][e];
		FOREACH(j, N) if (*j != e) A[e][*j] = A[l][*j] / A[l][e];
		b[e] = b[l] / A[l][e];
		memset(A[l], 0, sizeof(A[l]));
		b[l] = 0;

		B.erase(l); B.insert(e);
		N.erase(e); N.insert(l);
		FOREACH(_i, B) if (*_i != e) {
			int i = *_i;
			b[i] = b[i] - A[i][e] * b[e];
			FOREACH(j, N) 
				A[i][*j] -= A[i][e] * A[e][*j];
			A[i][e] = 0.0;
		}

		FOREACH(j, N) c[*j] -= c[e] * A[e][*j];
		v += c[e] * b[e];
		c[e] = 0.0;
	}

	void work() {
		int infinite = 0;
		while (1) {
			int e = -1; // 基本变量 -> 非基本变量
			FOREACH(j, N) if (sig(c[*j]) > 0) {
				e = *j;
				break;
			}
			if (e == -1) break;
			double delta = INF;
			int l; // 非基本变量 -> 基本变量
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
			if (delta == INF) {
				infinite = 1;
				break;
			}
			pivot(e, l);
		}
		for (int j = 0; j < n; j++) x[j] = b[j];
	}
} simplex;

double A[maxm][maxn], b[maxm], c[maxn];
int main() {
	A[0][0] = 1; A[0][1] = 1; b[0] = 20;
	A[1][0] = 2; A[1][1] = 1; b[1] = 25;
	c[0] = 3;	c[1] = 2;
	simplex.init(A, b, c, 2, 2);
	simplex.work();
	for (int j = 0; j < 2; j++) printf("%.2f\n", simplex.x[j]);
	printf("min = %.2f\n", simplex.v);

	return 0;
}
