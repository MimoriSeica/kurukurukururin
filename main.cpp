#include <bits/stdc++.h>

using namespace std;

#define REP(i,n) for(ll (i) = (0);(i) < (n);++i)
#define REV(i,n) for(ll (i) = (n) - 1;(i) >= 0;--i)
#define PB push_back
#define EB emplace_back
#define MP make_pair
#define FI first
#define SE second
#define SHOW1d(v,n) {REP(WW,n)cerr << v[WW] << ' ';cerr << endl << endl;}
#define SHOW2d(v,WW,HH) {REP(W_,WW){REP(H_,HH)cerr << v[W_][H_] << ' ';cerr << endl;}cerr << endl;}
#define ALL(v) v.begin(),v.end()
#define Decimal fixed<<setprecision(20)
#define INF 1000000000
#define LLINF 1000000000000000000LL
#define MOD 998244353

typedef long long ll;
typedef pair<ll,ll> P;
//--------geometry original ------------------
#define curr(PP, i) PP[i]
#define next(PP, i) PP[(i+1)%PP.size()]
#define diff(PP, i) (next(PP, i) - curr(PP, i))

typedef long long ll;
typedef pair<ll, ll> P;

const double EPS = 1e-7;
const double EPS_GIG = 1e-3;
const double PI = acos(-1.0);
typedef complex<double> point;
namespace std {
	bool operator < (const point& a, const point& b) {
		return real(a) != real(b) ? real(a) < real(b) : imag(a) < imag(b);
	}

	bool operator == (const point& a,const point& b) {
		return (abs(a.real() - b.real()) < EPS && abs(a.imag() - b.imag()) < EPS);
	}
}
struct circle {
	point p; double r;
	circle(){}
	circle(const point &p, double r) : p(p), r(r) { }
};
struct segment : public vector<point> {
	segment(const point &a, const point &b) {
		push_back(a); push_back(b);
	}
};

double cross(const point& a, const point& b) {
	return imag(conj(a)*b);
}

double dot(const point& a, const point& b) {
	return real(conj(a)*b);
}

/*
a → b で時計方向に折れて b → c
a → b で半時計方向に折れて b → c
a → b で逆を向いて a を通り越して b → c
a → b でそのまま b → c
a → b で逆を向いて b → c ( または b == c )
*/

int ccw(point a, point b, point c) {
	b -= a; c -= a;
	if (cross(b, c) > EPS)   return +1;       // counter clockwise
	if (cross(b, c) + EPS < 0)   return -1;       // clockwise
	if (dot(b, c) + EPS < 0)     return +2;       // c--a--b on line
	if (norm(b) < norm(c)) return -2;       // a--b--c on line
	return 0;
}

bool intersectLL(const segment &l, const segment &m) {
	return abs(cross(l[1] - l[0], m[1] - m[0])) > EPS || // non-parallel
		abs(cross(l[1] - l[0], m[0] - l[0])) < EPS;   // same line
}
bool intersectLS(const segment &l, const segment &s) {
	return cross(l[1] - l[0], s[0] - l[0])*       // s[0] is left of l
		cross(l[1] - l[0], s[1] - l[0]) < EPS; // s[1] is right of l
}
bool intersectLP(const segment &l, const point &p) {
	return abs(cross(l[1] - p, l[0] - p)) < EPS;
}
bool intersectSP(const segment &s, const point &p) {
	return ccw(s[0], p, s[1]) == -2; // triangle inequality
}
//端点の交差も考える
bool intersectSS(const segment &s, const segment &t) {
	if(intersectSP(s,t[0]) || intersectSP(s,t[1]) || intersectSP(t,s[0]) || intersectSP(t,s[1]))return true;
	return ccw(s[0], s[1], t[0]) * ccw(s[0], s[1], t[1]) + EPS <= 0 &&
		ccw(t[0], t[1], s[0]) * ccw(t[0], t[1], s[1]) + EPS <= 0;
}

point projection(const segment &l, const point &p) {
	double t = dot(p - l[0], l[0] - l[1]) / norm(l[0] - l[1]);
	return l[0] + t*(l[0] - l[1]);
}
point reflection(const segment &l, const point &p) {
	return p + 2. * (projection(l, p) - p);
}
double distanceLP(const segment &l, const point &p) {
	return abs(p - projection(l, p));
}
double distanceLL(const segment &l, const segment &m) {
	return intersectLL(l, m) ? 0 : distanceLP(l, m[0]);
}
double distanceLS(const segment &l, const segment &s) {
	if (intersectLS(l, s)) return 0;
	return min(distanceLP(l, s[0]), distanceLP(l, s[1]));
}
double distanceSP(const segment &s, const point &p) {
	const point r = projection(s, p);
	if (intersectSP(s, r)) return abs(r - p);
	return min(abs(s[0] - p), abs(s[1] - p));
}
double distanceSS(const segment &s, const segment &t) {
	if (intersectSS(s, t)) return 0;
	return min(min(distanceSP(s, t[0]), distanceSP(s, t[1])),
		min(distanceSP(t, s[0]), distanceSP(t, s[1])));
}
double distancePP(const point &a,const point &b){
	return abs(a-b);
}

//交点
point crosspoint(const segment &l, const segment &m) {
	double A = cross(l[1] - l[0], m[1] - m[0]);
	double B = cross(l[1] - l[0], l[1] - m[0]);
	if (abs(A) < EPS && abs(B) < EPS) return m[0]; // same line
	if (abs(A) < EPS) return point(INF,INF); // !!!PRECONDITION NOT SATISFIED!!!
	return m[0] + B / A * (m[1] - m[0]);
}

//凸包
vector<point> convex_hull(vector<point> ps) {
	int n = ps.size(), k = 0;
	sort(ps.begin(), ps.end());
	vector<point> ch(2*n);
	for (int i = 0; i < n; ch[k++] = ps[i++]) // lower-hull
		while (k >= 2 && ccw(ch[k-2], ch[k-1], ps[i]) == -1) --k;
	for (int i = n-2, t = k+1;i >= 0; ch[k++] = ps[i--]) // upper-hull
		while (k >= t && ccw(ch[k-2], ch[k-1], ps[i]) == -1) --k;
	ch.resize(k - 1);
	return ch;
}

/*多角形内包判定
OUT:0
ON:1
IN:2
*/

int contains(const vector<point>& Poly, const point& p) {
	bool in = false;
	for (int i = 0; i < Poly.size(); ++i) {
		point a = curr(Poly, i) - p, b = next(Poly, i) - p;
		if (imag(a) > imag(b)) swap(a, b);
		if (imag(a) + EPS <= 0 && EPS < imag(b))
			if (cross(a, b) + EPS < 0) in = !in;
		if (abs(cross(a, b)) < EPS && dot(a, b) + EPS <= 0) return 1;
	}
	return in ? 2 : 0;
}

//見えるか(可視グラフ用)
bool block_off(const point &a, const point &b, const vector<point> &obj) {
  point m = (a+b)/2.0;
  bool on = false, in = false;
  for (int j = 0; j < obj.size(); ++j) {
    point c = curr(obj,j), d = next(obj,j);
    if (imag(d) < imag(c)) swap(c, d);
    if (cross(a-c,b-c) * cross(a-d,b-d) < 0 &&    // strictly intersect.
        cross(c-a,d-a) * cross(c-b,d-b) < 0) return true;
    if (cross(a-c,b-c) == 0 && dot(a-c,b-c) < 0) return true;
    if (imag(c) <= imag(m) && imag(m) < imag(d))  // strictly contain.
      if (cross(c-m,d-m) < 0) in = !in;
    if (cross(c-m,d-m) == 0 && dot(c-m,d-m) <= 0) on = true;
  }
  return !on && in;
}

//面積
double area(const vector<point>& p) {
	double A = 0;
	for (int i = 0; i < p.size(); ++i)
		A += cross(curr(p, i), next(p, i));
	return A / 2.;
}

//角度足し算
double add_rad(double a,double b){
	double ret = a + b;
	if(ret > PI)ret -= 2 * PI;
	return ret;
}

//なす角(vector)
double angle(const point &a,const point &b) {
	auto tmp = abs(arg(a) - arg(b));
	return min(tmp, 2 * PI - tmp);
}

double angle(const segment &s1,const segment &s2) {
	return angle(s1[1] - s1[0], s2[1] - s2[0]);
}

//点の回転
point rotate(const point &p, double rad) {
	double x = p.real() * cos(rad) - p.imag() * sin(rad);
	double y = p.imag() * cos(rad) + p.real() * sin(rad);
	return point(x, y);
}

//並行
bool isParallel(const point &a, const point &b){
    return abs(cross(a,b)) < EPS;
}
bool isParallel(const segment &a, const segment &b){
    return isParallel(a[1]-a[0], b[1]-b[0]);
}
//直行
bool isOrthogonal(const point &a,const point &b){
	return abs(angle(a,b) - PI / 2) < EPS;
}
bool isOrthogonal(const segment &a,const segment &b){
	return isOrthogonal(a[1]-a[0],b[1]-b[0]);
}
//凸包判定
bool isConvex(vector<point> poly){
	int sz = poly.size();
	REP(i,sz){
		int tmp = ccw(poly[i],poly[(i+1)%sz],poly[(i+2)%sz]);
		if(tmp == -1){
			return false;
		}
	}
	return true;
}

double convex_diameter(const vector<point> &pt) {
  const int n = pt.size();
  int is = 0, js = 0;
  for (int i = 1; i < n; ++i) {
    if (imag(pt[i]) > imag(pt[is])) is = i;
    if (imag(pt[i]) < imag(pt[js])) js = i;
  }
  double maxd = norm(pt[is]-pt[js]);

  int i, maxi, j, maxj;
  i = maxi = is;
  j = maxj = js;
  do {
    if (cross(diff(pt,i), diff(pt,j)) >= 0) j = (j+1) % n;
    else i = (i+1) % n;
    if (norm(pt[i]-pt[j]) > maxd) {
      maxd = norm(pt[i]-pt[j]);
      maxi = i; maxj = j;
    }
  } while (i != is || j != js);
  return sqrt(maxd); /* farthest pair is (maxi, maxj). */
}

vector<point> convex_cut(const vector<point> P, const segment& l) {
  vector<point> Q;
  for (int i = 0; i < P.size(); ++i) {
    point A = curr(P, i), B = next(P, i);
    if (ccw(l[0], l[1], A) != -1) Q.push_back(A);
    if (ccw(l[0], l[1], A)*ccw(l[0], l[1], B) < 0)
      Q.push_back(crosspoint(segment(A, B), l));
  }
  return Q;
}

//扇型、中心と二つの端点
struct sector {
	point o;
	point a, b;
	sector(){}
	sector(point O, point A, point B) :o(O), a(A), b(B) {}
};

double L, sx, sy, gx, gy;
int point_size = 0, n, r;
vector<int> dist;
vector<int> goal_list;
vector<segment> seg_list;
vector<vector<P>> v;
// 各回転時点での点候補(点, id)
vector<vector<pair<point, int>>> r_point_list(20);
// 各回転時点での四角形
vector<vector<vector<point>>> r_square_list(20);
// 各回転時点での線分
vector<vector<segment>> r_segment_list(20);
// 各回転時点での扇
vector<vector<sector>> r_sector_list(20);

bool can_rotate(point p, int rad_num);

bool add_point_list(int rad_num, point p){
	REP(i, r_square_list[rad_num].size()){
		if(contains(r_square_list[rad_num][i], p) == 2)return false;
	}
	r_point_list[rad_num].EB(p, point_size);
	point_size++;
	return true;
}

void make_r_segment_list_square_list(int rad_num) {
	double rad = (double)rad_num * PI / r;
	double next_rad = (double)((rad_num + 1) % r) * PI / r;
	point vec = rotate(point(L, 0), rad);
	point next_vec = rotate(point(L, 0), next_rad);

	REP(i, n) {
		segment now = seg_list[i];
		point A = now[0];
		point B = now[1];
		point now_v = B - A;
		vector<point> tmp;
		if (isParallel(now_v,vec)) {
			point rotated_vec = rotate(point(EPS_GIG,0),add_rad(rad,PI/2));
			point C, D;
			if(norm(A - (B + vec)) > norm(A - (B - vec)))C = B + vec;
			else C = B - vec;
			if(norm(B - (A + vec)) > norm(B - (A - vec)))D = A + vec;
			else D = A - vec;

			tmp.PB(C + rotated_vec);
			tmp.PB(C - rotated_vec);
			tmp.PB(D + rotated_vec);
			tmp.PB(D - rotated_vec);
		}
		else {
			tmp.PB(A + vec);
			tmp.PB(A - vec);
			tmp.PB(B + vec);
			tmp.PB(B - vec);
		}
		tmp = convex_hull(tmp);
		r_square_list[rad_num].EB(tmp);
		REP(j, 4) {
			r_segment_list[rad_num].EB(tmp[j], tmp[(j + 1) % 4]);
		}
	}
}

void make_r_sector_list(int rad_num){
	double rad = (double)rad_num * PI / r;
	double next_rad = (double)((rad_num + 1) % r) * PI / r;
	point vec = rotate(point(L, 0), rad);
	point next_vec = rotate(point(L, 0), next_rad);

	REP(i, seg_list.size()){
		REP(j, 2){
			auto p = seg_list[i][j];
			r_sector_list[rad_num].EB(p, p + vec, p + next_vec);
			r_sector_list[rad_num].EB(p, p - vec, p - next_vec);
		}
	}
}

void make_r_point_list(int rad_num) {
	add_point_list(rad_num, point(sx, sy));
	if(add_point_list(rad_num, point(gx, gy))){
		goal_list.EB(point_size-1);
	}

	REP(i, r_segment_list[rad_num].size()){
		auto seg_a = r_segment_list[rad_num][i];
		REP(j, i){
			auto seg_b = r_segment_list[rad_num][j];
			if(!intersectSS(seg_a, seg_b))continue;
			auto p = crosspoint(seg_a, seg_b);
			add_point_list(rad_num, p);
			if(can_rotate(p, rad_num)){
				add_point_list((rad_num + r - 1) % r, p);
			}
			if(can_rotate(p, (rad_num + 1) % r)){
				add_point_list((rad_num + 1) % r, p);
			}
		}
	}
	REP(i, r_square_list[rad_num].size()){
		REP(j, 4){
			add_point_list(rad_num, r_square_list[rad_num][i][j]);
			if(can_rotate(r_square_list[rad_num][i][j], rad_num)){
				add_point_list((rad_num + r - 1) % r, r_square_list[rad_num][i][j]);
			}
			if(can_rotate(r_square_list[rad_num][i][j], (rad_num + 1) % r)){
				add_point_list((rad_num + 1) % r, r_square_list[rad_num][i][j]);
			}
		}
	}
}

bool check_visible(point p_a, point p_b, int rad_num) {
	segment seg = segment(p_a, p_b);
	point mid = (p_a + p_b) * 0.5;
	REP(i, r_square_list[rad_num].size()){
		int cou = 0;
		bool flag = false;
		auto sq = r_square_list[rad_num][i];
		if(contains(sq, mid) == 2)return false;
		REP(j, 4){
			if(ccw(sq[j], sq[(j+1)%4], seg[0]) *
			   ccw(sq[j], sq[(j+1)%4], seg[1]) + EPS <= 0 &&
			   ccw(seg[0], seg[1], sq[j]) *
			   ccw(seg[0], seg[1], sq[(j+1)%4]) + EPS <= 0)return false;
		}
	}
	return true;
}

void make_visible_graph(int rad_num) {
	REP(i, r_point_list[rad_num].size()){
		auto p_a = r_point_list[rad_num][i].FI;
		auto id_a = r_point_list[rad_num][i].SE;
		REP(j, i){
			auto p_b = r_point_list[rad_num][j].FI;
			auto id_b = r_point_list[rad_num][j].SE;
			if(check_visible(p_a, p_b, rad_num)){
				v[id_a].EB(id_b, 0);
				v[id_b].EB(id_a, 0);
			}
		}
	}
}

bool contain_sector(sector &sec,point &p){
	point o = sec.o;
	point a = sec.a;
	point b = sec.b;
	if(distancePP(p,o) > L + EPS)return false;
	if(intersectSP(segment(o,a), p))return true;
	if(intersectSP(segment(o,b), p))return true;
	point vec = p - o;
	point vecA = a - o;
	point vecB = b - o;
	if(angle(vec,vecA) < angle(vecA,vecB) && angle(vec,vecB) < angle(vecA,vecB))return true;
	return false;
}

bool can_rotate(point p, int rad_num) {
	int next_rad = (rad_num + 1) % r;
	REP(i, r_square_list[rad_num].size()){
		if(contains(r_square_list[rad_num][i], p) == 2)return false;
	}
	REP(i, r_square_list[next_rad].size()){
		if(contains(r_square_list[next_rad][i], p) == 2)return false;
	}
	REP(i, r_sector_list[rad_num].size()){
		if(contain_sector(r_sector_list[rad_num][i], p))return false;
	}
	return true;
}

void make_rotate_graph(int rad_num) {
	int next_rad = (rad_num + 1) % r;
	REP(i, r_point_list[rad_num].size()){
		point p_a = r_point_list[rad_num][i].FI;
		int id_a = r_point_list[rad_num][i].SE;
		REP(j, r_point_list[next_rad].size()){
			point p_b = r_point_list[next_rad][j].FI;
			int id_b = r_point_list[next_rad][j].SE;

			if(norm(p_a - p_b) > EPS)continue;
			if(can_rotate(p_a, rad_num)){
				v[id_a].EB(id_b, 1);
			}
		}
	}
}

int Dijkstra() {
	priority_queue<P, vector<P>, greater<>> pq;
	pq.push(MP(0, 0));

	while(!pq.empty()){
		int d = pq.top().FI;
		int id = pq.top().SE;
		pq.pop();
		REP(i, v[id].size()){
			int aite = v[id][i].FI;
			int cost = v[id][i].SE;
			if(dist[aite] > d + cost){
				dist[aite] = d + cost;
				pq.push(MP(dist[aite], aite));
			}
		}
	}

	int ret = INF;
	REP(i, goal_list.size()){
		ret = min(ret, dist[goal_list[i]]);
	}
	return ret;
}

void output_visible_graph() {
	cout << r << endl;
	REP(rad_num, r){
		int m = r_point_list[rad_num].size();
		cout << m << endl;
		REP(i, m){
			auto p = r_point_list[rad_num][i].FI;
			bool flag = can_rotate(p, rad_num);
			cout << p.real() << " " << p.imag() << " " << flag << endl;
		}
		m = r_segment_list[rad_num].size();
		cout << m << endl;
		REP(i, m){
			segment seg = r_segment_list[rad_num][i];
			cout << seg[0].real() << " " << seg[0].imag() << " " << seg[1].real() << " " << seg[1].imag() << endl;
		}
		vector<segment> g;
		REP(i, r_point_list[rad_num].size()){
			auto p_a = r_point_list[rad_num][i].FI;
			auto id_a = r_point_list[rad_num][i].SE;
			REP(j, i){
				auto p_b = r_point_list[rad_num][j].FI;
				auto id_b = r_point_list[rad_num][j].SE;
				if(check_visible(p_a, p_b, rad_num)){
					g.EB(p_a, p_b);
				}
			}
		}
		m = g.size();
		cout << m << endl;
		REP(i, m){
			segment seg = g[i];
			cout << seg[0].real() << " " << seg[0].imag() << " " << seg[1].real() << " " << seg[1].imag() << endl;
		}
	}
}

int main(){
	cin.tie(0);cout.tie(0);ios::sync_with_stdio(false);
	cin >> L >> r;
	cin >> sx >> sy;
	cin >> gx >> gy;

	cin >> n;
	REP(i, n){
		double x1, y1, x2, y2;cin >> x1 >> y1 >> x2 >> y2;
		seg_list.EB(point(x1, y1), point(x2, y2));
	}
	clock_t start = clock();

	REP(i, r){
		make_r_segment_list_square_list(i);
		make_r_sector_list(i);
	}
	REP(i, r){
		make_r_point_list(i);
	}
	v.resize(point_size);
	dist.resize(point_size);
	REP(i, r){
		make_visible_graph(i);
		make_rotate_graph(i);
	}
/*
	output_visible_graph();
	return 0;
*/
	cout << "point_size " << point_size << endl;

	REP(i, point_size)dist[i] = INF;
	dist[0] = 0;
	int ans = Dijkstra();

	if(ans == INF)cout << -1 << endl;
	else cout << ans << endl;
	clock_t end = clock();
	double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
	cout << "time " << time << endl;
	return 0;
}
