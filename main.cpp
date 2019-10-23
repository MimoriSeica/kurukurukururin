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
#define eq(n,m) (abs((n)-(m)) < EPS)

typedef long long ll;
typedef pair<ll, ll> P;

const double EPS = 1e-8;
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

// 扇型、中心と半径、二つの端点
// 現在中心角が180未満の前提
struct sector {
	point o;
	point a, b;
	double r;
	sector(){}
	sector(point O, point A, point B, double _r) :o(O), a(A), b(B), r(_r) {}
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

//角度足し算
double add_rad(double a,double b){
	double ret = a + b;
	if(ret > 2 * PI)ret -= 2 * PI;
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
	if (dot(b, c) - EPS < 0)     return +2;       // c--a--b on line
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
	auto a = s[0] - p;
	auto b = s[1] - p;
	return (abs(cross(a, b)) < EPS && dot(a, b) <= EPS); // triangle inequality
}
//端点の交差も考える
bool intersectSS(const segment &s, const segment &t) {
	return ccw(s[0], s[1], t[0]) * ccw(s[0], s[1], t[1]) <= 0 &&
		ccw(t[0], t[1], s[0]) * ccw(t[0], t[1], s[1]) <= 0;
}
//端点の交差hは考えない
bool strictIntersectSS(const segment &s, const segment &t) {
	return ccw(s[0], s[1], t[0]) * ccw(s[0], s[1], t[1]) == -1 &&
		ccw(t[0], t[1], s[0]) * ccw(t[0], t[1], s[1]) == -1;
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

/*多角形内包判定
half-line crossing method
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
			if (cross(a, b) < 0) in = !in;
		if (abs(cross(a, b)) < EPS && dot(a, b) <= EPS) return 1;
	}
	return in ? 2 : 0;
}

// 厳密に入っている時のみTrue
bool contain_sector(const sector &sec, point &p){
	if(abs(p - sec.o) + EPS > sec.r)return false;
	point vec = p - sec.o;
	point vecA = sec.a - sec.o;
	point vecB = sec.b - sec.o;
	if(angle(vec, vecA) < angle(vecA, vecB) && angle(vec, vecB) < angle(vecA, vecB))return true;
	return false;
}

//交点
point crosspointSS(const segment &l, const segment &m) {
	double A = cross(l[1] - l[0], m[1] - m[0]);
	double B = cross(l[1] - l[0], l[1] - m[0]);
	if (abs(A) < EPS && abs(B) < EPS) return m[0]; // same line
	if (abs(A) < EPS) return point(INF,INF); // !!!PRECONDITION NOT SATISFIED!!!
	return m[0] + B / A * (m[1] - m[0]);
}

vector<point> crosspointCL(const circle &c, const segment &l) {
	auto ret = vector<point>(2, point(INF, INF));
	auto pro_p = projection(l, c.p);
	auto dist = distanceLP(l, c.p);
	if(abs(dist - c.r) < EPS){
		ret[0] = pro_p;
		return ret;
	}
	if(c.r < dist){
		return ret;
	}
	point vec = (l[1] - l[0]) * sqrt(c.r * c.r - dist * dist) / abs(l[1] - l[0]);
	ret[0] = pro_p + vec;
	ret[1] = pro_p - vec;
	return ret;
}

vector<point> crosspointCC(const circle c1, const circle c2) {
	auto ret = vector<point>(2, point(INF, INF));
	auto dist = abs(c2.p - c1.p);
	if(eq(dist, c1.r + c2.r) || eq(dist, abs(c2.r - c1.r))){
		auto tmp = c2.p - c1.p;
		ret[0] = c1.p + tmp * (c1.r / dist);
		return ret;
	}
	if(c1.r + c2.r < dist || dist < abs(c1.r - c2.r)){
		return ret;
	}
	auto alpha = acos((c1.r * c1.r + dist * dist - c2.r * c2.r) / (2 * c1.r * dist));
	auto theta = atan2(c2.p.imag() - c1.p.imag(), c2.p.real() - c1.p.real());
	ret[0] = c1.p + point(cos(theta + alpha) * c1.r, sin(theta + alpha) * c1.r);
	ret[1] = c1.p + point(cos(theta - alpha) * c1.r, sin(theta - alpha) * c1.r);
	return ret;
}

bool isOnSector(const sector sec, const point p) {
	point vec = p - sec.o;
	point vecA = sec.a - sec.o;
	point vecB = sec.b - sec.o;
	if(eq(angle(vec, vecA) + angle(vec, vecB), angle(vecA, vecB)))return true;
	return false;
}

vector<point> crosspointSecS(const sector sec, const segment s) {
	circle c = circle(sec.o, sec.r);
	auto ret = crosspointCL(c, s);
	point inf = point(INF, INF);
	REP(i, 2){
		if(eq(ret[i], inf))continue;
		if(!isOnSector(sec, ret[i])){
			ret[i] = inf;
			continue;
		}
		if(!intersectSP(s, ret[i])){
			ret[i] = inf;
		}
	}
	return ret;
}
vector<point> crosspointSecSec(const sector sec1, const sector sec2) {
	circle c1 = circle(sec1.o, sec1.r);
	circle c2 = circle(sec2.o, sec2.r);
	auto ret = crosspointCC(c1, c2);
	point inf = point(INF, INF);
	REP(i, 2){
		if(!isOnSector(sec1, ret[i])){
			ret[i] = inf;
			continue;
		}
		if(!isOnSector(sec2, ret[i])){
			ret[i] = inf;
		}
	}
	return ret;
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

//見えるか(可視グラフ用)
bool block_off(const point &a, const point &b, const vector<point> &obj) {
  point m = (a + b) * 0.5;
  bool on = false, in = false;
  for (int j = 0; j < obj.size(); ++j) {
    point c = curr(obj,j), d = next(obj,j);
    if (imag(d) < imag(c)) swap(c, d);
    if (cross(a-c,b-c) * cross(a-d,b-d) < 0 &&    // strictly intersect.
        cross(c-a,d-a) * cross(c-b,d-b) < 0) return true;
    if (cross(a-c,b-c) == 0 && dot(a-c,b-c) < 0) return true;
    if (imag(c) <= imag(m) && imag(m) < imag(d))  // strictly contain.
      if (cross(c-m,d-m) < 0) in = !in;
    if (cross(c-m,d-m) == 0 && dot(c-m,d-m) <= EPS) on = true;
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
      Q.push_back(crosspointSS(segment(A, B), l));
  }
  return Q;
}

double L, sx, sy, gx, gy;
int point_size = 0, n, r;
vector<int> dist(10000);
vector<int> goal_list;
vector<segment> seg_list;
vector<vector<P>> v(10000);
// 各回転時点での点候補(点, id)
vector<vector<pair<point, int>>> r_point_list(20);
// 回転代表点でないものの点数
vector<int> r_def_point_size(20);
// 各回転時点での四角形
vector<vector<vector<point>>> r_square_list(20);
// 各回転時点での線分
vector<vector<segment>> r_segment_list(20);
// 各回転時点での扇
vector<vector<sector>> r_sector_list(20);

bool can_rotate(point p, int rad_num);
bool check_visible(point p_a, point p_b, int rad_num);

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
	point vec = point(cos(rad), sin(rad)) * L;

	REP(i, n) {
		segment now = seg_list[i];
		point A = now[0];
		point B = now[1];
		point now_v = A - B;
		vector<point> tmp;
		if (isParallel(now_v,vec)) {
			point rotated_vec = point(-sin(rad), cos(rad)) * EPS_GIG;
			point C, D;
			if(abs(A - (B + vec)) > abs(A - (B - vec)))C = B + vec;
			else C = B - vec;
			if(abs(B - (A + vec)) > abs(B - (A - vec)))D = A + vec;
			else D = A - vec;

			tmp.PB(C + rotated_vec);
			tmp.PB(C - rotated_vec);
			tmp.PB(D + rotated_vec);
			tmp.PB(D - rotated_vec);
		}
		else {
			auto eps_vec = now_v * EPS_GIG / abs(now_v);
			tmp.PB(A + vec + eps_vec);
			tmp.PB(A - vec + eps_vec);
			tmp.PB(B + vec - eps_vec);
			tmp.PB(B - vec - eps_vec);
		}
		tmp = convex_hull(tmp);
		r_square_list[rad_num].EB(tmp);
		REP(j, 4) {
			r_segment_list[rad_num].EB(tmp[j], tmp[(j + 1) % 4]);
		}
	}
}

void add_sector(vector<sector> &sec_v, sector sec, int rad_num, int square_num){
	if(contains(r_square_list[rad_num][square_num], sec.b) != 0){
		swap(sec.a, sec.b);
	}
	sec_v.EB(sec);
}

void make_r_sector_list(int rad_num){
	double rad = (double)rad_num * PI / r;
	double next_rad = (double)(rad_num + 1) * PI / r;
	point vec = point(cos(rad), sin(rad)) * L;
	point next_vec = point(cos(next_rad), sin(next_rad)) * L;

	REP(i, seg_list.size()){
		REP(j, 2){
			auto p = seg_list[i][j];
			auto sec = sector(p, p + vec, p + next_vec, L);
			add_sector(r_sector_list[rad_num], sec, rad_num, i);
			sec = sector(p, p - vec, p - next_vec, L);
			add_sector(r_sector_list[rad_num], sec, rad_num, i);
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
			auto p = crosspointSS(seg_a, seg_b);
			add_point_list(rad_num, p);
		}
	}
	REP(i, r_square_list[rad_num].size()){
		REP(j, 4){
			add_point_list(rad_num, r_square_list[rad_num][i][j]);
		}
	}
	r_def_point_size[rad_num] = r_point_list[rad_num].size();
}

void add_rotate_point(point p_a, int id_a, int rad_num){
	int next_rad = (rad_num + 1) % r;
	REP(i, r_def_point_size[rad_num]){
		auto p_b = r_point_list[rad_num][i].FI;
		auto id_b = r_point_list[rad_num][i].SE;
		if(check_visible(p_a, p_b, rad_num)){
			v[id_b].EB(id_a, 0);
			v[id_a].EB(id_b, 0);
		}
	}
	REP(i, r_def_point_size[next_rad]){
		auto p_b = r_point_list[next_rad][i].FI;
		auto id_b = r_point_list[next_rad][i].SE;
		if(check_visible(p_a, p_b, next_rad)){
			v[id_a].EB(id_b, 1);
		}
	}
}

void make_rotate_point(int rad_num) {
	int next_rad = (rad_num + 1) % r;
	int pre_rad = (rad_num + r - 1) % r;
	REP(i, r_def_point_size[rad_num]){
		auto p_a = r_point_list[rad_num][i].FI;
		auto id_a = r_point_list[rad_num][i].SE;
		if(can_rotate(p_a, rad_num)){
			REP(j, r_def_point_size[next_rad]){
				auto p_b = r_point_list[next_rad][j].FI;
				auto id_b = r_point_list[next_rad][j].SE;
				if(check_visible(p_a, p_b, next_rad)){
					v[id_a].EB(id_b, 1);
				}
			}
		}
		if(can_rotate(p_a, pre_rad)){
			REP(j, r_def_point_size[pre_rad]){
				auto p_b = r_point_list[pre_rad][j].FI;
				auto id_b = r_point_list[pre_rad][j].SE;
				if(check_visible(p_a, p_b, pre_rad)){
					v[id_b].EB(id_a, 1);
				}
			}
		}
	}
	REP(i, r_segment_list[rad_num].size()) {
		auto seg_a = r_segment_list[rad_num][i];
		REP(j, r_segment_list[next_rad].size()){
			auto seg_b = r_segment_list[next_rad][j];
			if(!intersectSS(seg_a, seg_b))continue;
			auto p_a = crosspointSS(seg_a, seg_b);
			if(!can_rotate(p_a, rad_num))continue;
			auto id_a = point_size;
			add_point_list(rad_num, p_a);
			add_rotate_point(p_a, id_a, rad_num);
		}
	}
	REP(i, r_sector_list[rad_num].size()){
		auto sc1 = r_sector_list[rad_num][i];
		point inf = point(INF, INF);
		REP(j, i){
			auto sc2 = r_sector_list[rad_num][j];
			auto tmp = crosspointSecSec(sc1, sc2);
			REP(ii, 2){
				if(eq(inf, tmp[ii]))continue;
				if(!can_rotate(tmp[ii], rad_num))continue;
				auto p_a = tmp[ii];
				auto id_a = point_size;
				add_point_list(rad_num, p_a);
				add_rotate_point(p_a, id_a, rad_num);
			}
		}
		REP(j, r_segment_list[rad_num].size()){
			auto seg = r_segment_list[rad_num][j];
			auto tmp = crosspointSecS(sc1, seg);
			REP(ii, 2){
				if(eq(inf, tmp[ii]))continue;
				if(!can_rotate(tmp[ii], rad_num))continue;
				auto p_a = tmp[ii];
				auto id_a = point_size;
				add_point_list(rad_num, p_a);
				add_rotate_point(p_a, id_a, rad_num);
			}
		}
		REP(j, r_segment_list[next_rad].size()){
			auto seg = r_segment_list[next_rad][j];
			auto tmp = crosspointSecS(sc1, seg);
			REP(ii, 2){
				if(eq(inf, tmp[ii]))continue;
				if(!can_rotate(tmp[ii], rad_num))continue;
				auto p_a = tmp[ii];
				auto id_a = point_size;
				add_point_list(rad_num, p_a);
				add_rotate_point(p_a, id_a, rad_num);
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
			segment tmp_seg = segment(sq[j], sq[(j+1)%4]);
			if(strictIntersectSS(seg, tmp_seg))return false;
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

bool can_rotate(point p, int rad_num) {
	int next_rad = (rad_num + 1) % r;
	REP(i, r_square_list[rad_num].size()){
		if(contains(r_square_list[rad_num][i], p) == 2)return false;
	}
	REP(i, r_square_list[next_rad].size()){
		if(contains(r_square_list[next_rad][i], p) == 2)return false;
	}
	REP(i, r_sector_list[rad_num].size()){
		auto sec = r_sector_list[rad_num][i];
		if(contain_sector(sec, p))return false;
	}
	return true;
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
		m = r_segment_list[(rad_num + 1) % r].size();
		cout << m << endl;
		REP(i, m){
			segment seg = r_segment_list[(rad_num + 1) % r][i];
			cout << seg[0].real() << " " << seg[0].imag() << " " << seg[1].real() << " " << seg[1].imag() << endl;
		}
		m = r_sector_list[rad_num].size();
		cout << 2 * m << endl;
		REP(i, m){
			auto sec = r_sector_list[rad_num][i];
			cout << sec.o.real() << " " << sec.o.imag() << " " << sec.a.real() << " " << sec.a.imag() << endl;
			cout << sec.o.real() << " " << sec.o.imag() << " " << sec.b.real() << " " << sec.b.imag() << endl;
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
	exit(0);
}

int main(){
	cin.tie(0);cout.tie(0);ios::sync_with_stdio(false);
	cin >> L >> r;L += EPS_GIG;
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
	REP(i, r){
		make_visible_graph(i);
	}
	REP(i, r){
		make_rotate_point(i);
	}

	//output_visible_graph();

	//cout << "point_size " << point_size << endl;

	REP(i, point_size)dist[i] = INF;
	dist[0] = 0;
	int ans = Dijkstra();

	if(ans == INF)cout << -1 << endl;
	else cout << ans << endl;
	clock_t end = clock();
	double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
	//cout << "time " << time << endl;
	return 0;
}
