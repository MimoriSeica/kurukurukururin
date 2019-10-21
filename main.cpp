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
#define X real()
#define Y imag()

typedef long long ll;
typedef pair<ll, ll> P;

const double EPS = 1e-10;
const double EPS_GIG = 1e-5;
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
	if (cross(b, c) > 0)   return +1;       // counter clockwise
	if (cross(b, c) < 0)   return -1;       // clockwise
	if (dot(b, c) < 0)     return +2;       // c--a--b on line
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
	return abs(s[0] - p) + abs(s[1] - p) - abs(s[1] - s[0]) < EPS; // triangle inequality
}
//端点の交差も考える
bool intersectSS(const segment &s, const segment &t) {
	if(intersectSP(s,t[0]) || intersectSP(s,t[1]) || intersectSP(t,s[0]) || intersectSP(t,s[1]))return true;
	return ccw(s[0], s[1], t[0])*ccw(s[0], s[1], t[1]) <= 0 &&
		ccw(t[0], t[1], s[0])*ccw(t[0], t[1], s[1]) <= 0;
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
		if (imag(a) <= 0 && 0 < imag(b))
			if (cross(a, b) < 0) in = !in;
		if (cross(a, b) == 0 && dot(a, b) <= 0) return 1;
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

//中心と二つの端点
struct sector {
	point o;
	point a, b;
	sector(){}
	sector(point O, point A, point B) :o(O), a(A), b(B) {}
};

//回転不可領域
struct rotate_poly{
	rotate_poly(){}
	sector sec[4];
	vector<point> p;
};

int node_index;
vector<int> g_index;
vector < vector<pair<int, int>>> master_v(120 * 600);
vector<vector<segment>> r_seg(12);
vector<vector<segment>> rotate_seg(12);
vector<vector<int>> rn_index(12);
vector<vector<point>> rn_pos(12);
vector<vector<rotate_poly>> r_poly(12);
vector<vector<vector<point>>> r_squ(12);
vector<segment> in_seg;
point st,go;
vector<int> master_d(12*600);
double L;
int r, sx, sy, gx, gy, n;
set<pair<pair<int,int>,int>> edgeset;

bool contain_sector(sector &sec,point &p){
	point o = sec.o;
	point a = sec.a;
	point b = sec.b;
	if(distancePP(p,o) > L )return false;
	point vec = p - o;
	point vecA = a - o;
	point vecB = b - o;
	if(angle(vec,vecA) <= angle(vecA,vecB) && angle(vec,vecB) <= angle(vecA,vecB))return true;
	return false;
}

void check_visible(int index,point &p,int p_ind,int dist){
	REP(i,rn_pos[index].size()){
		if(rn_index[index][i] == p_ind)continue;
		int a_ind = rn_index[index][i];
		point pa = rn_pos[index][i];
		segment n_seg = segment(pa,p);
		bool flag = true;
		REP(k,r_squ[index].size()){
			if(block_off(p,pa,r_squ[index][k])){
				flag = false;
				break;
			}
		}
		if(flag){
			if(edgeset.count(MP(MP(min(a_ind,p_ind),max(a_ind,p_ind)),dist))){
				continue;
			}
			edgeset.insert(MP(MP(min(a_ind,p_ind),max(a_ind,p_ind)),dist));
			if(dist != 1)master_v[a_ind].PB(MP(p_ind,-dist));
			if(dist != -1)master_v[p_ind].PB(MP(a_ind,dist));
		}
	}
}

rotate_poly make_rotate_poly(point vec,point next_vec,segment seg,int num){
	point now_vec = seg[1] - seg[0];
	point r_vec;
	double radA = angle(now_vec,vec);
	double radB = angle(now_vec,next_vec);
	if((radA <= PI / 2 && radB >= PI / 2)||(radA >= PI/2 && radB <= PI/2)){
		double tmp_rad = add_rad(atan2(now_vec.imag(),now_vec.real()),PI/2);
		r_vec = rotate(point(L,0),tmp_rad);
	}
	else{
		if(abs(radA - PI / 2) < abs(radB - PI / 2)){
			r_vec = vec;
		}
		else{
			r_vec = next_vec;
		}
	}
	rotate_poly ret = rotate_poly();

	vector<point> tmp;
	tmp.PB(seg[0] + r_vec);
	tmp.PB(seg[0] - r_vec);
	tmp.PB(seg[1] - r_vec);
	tmp.PB(seg[1] + r_vec);
	tmp = convex_hull(tmp);
	ret.p = tmp;

	REP(i,4){
		rotate_seg[num].PB(segment(tmp[i],tmp[(i+1)%4]));
	}

	ret.sec[0] = sector(seg[0],seg[0] + vec,seg[0] + next_vec);
	ret.sec[1] = sector(seg[0],seg[0] - vec,seg[0] - next_vec);
	ret.sec[2] = sector(seg[1],seg[1] + vec,seg[1] + next_vec);
	ret.sec[3] = sector(seg[1],seg[1] - vec,seg[1] - next_vec);

	return ret;
}
void make_visible_graph(int num) {

	//スタートとゴール
	rn_pos[num].PB(st);
	rn_pos[num].PB(go);
	rn_index[num].PB(node_index);
	node_index++;
	rn_index[num].PB(node_index);
	g_index.PB(node_index);
	node_index++;

	//不可領域の作成
	double rad = num * PI / (double)r;
	double next_rad = ((num + 1) % r) * PI / (double)r;
	point vec = rotate(point(L, 0), rad);
	point next_vec = rotate(point(L, 0), next_rad);


	REP(i, n) {
		segment now = in_seg[i];
		point A = now[0];
		point B = now[1];
		point now_v = B - A;
		vector<point> tmp;
		if (isParallel(now_v,vec)) {
			point rotated_vec = rotate(point(EPS_GIG,0),add_rad(rad,PI/2));
			point C = A + vec;
			point D = B - vec;
			A -= vec;
			B += vec;
			if(abs(A - B) < abs(C - D)){A = C;B = D;}
			tmp.PB(A + rotated_vec);
			tmp.PB(A - rotated_vec);
			tmp.PB(B + rotated_vec);
			tmp.PB(B - rotated_vec);
		}
		else {
			tmp.PB(A + vec);
			tmp.PB(A - vec);
			tmp.PB(B + vec);
			tmp.PB(B - vec);
		}
		tmp = convex_hull(tmp);
		r_squ[num].PB(tmp);
		REP(j, 4) {
			r_seg[num].PB(segment(tmp[j], tmp[(j + 1) % 4]));
		}
	}

	REP(i, r_seg[num].size()) {
		for(int j = i+1;j <  r_seg[num].size();j++) {
			segment Aseg = r_seg[num][i];
			segment Bseg = r_seg[num][j];
			point tmp = crosspoint(Aseg,Bseg);
			if(tmp == point(INF,INF))continue;
			rn_pos[num].PB(tmp);
			rn_index[num].PB(node_index);
			node_index++;
		}
	}

	//可視グラフ全探索
	REP(i,rn_pos[num].size()){
		check_visible(num,rn_pos[num][i],rn_index[num][i],0);
	}
}

void Dijkstra() {
	REP(i,master_d.size())master_d[i] = INF;
	master_d[0] = 0;
	priority_queue<P,vector<P>,greater<P>> pq;
	pq.push(MP(0,0));
	while(!pq.empty()){
		P now = pq.top();pq.pop();
		int now_d = now.FI;
		int node = now.SE;
		if(now_d > master_d[node])continue;
		REP(i,master_v[node].size()){
			int aite = master_v[node][i].FI;
			int dist = master_v[node][i].SE;
			if(master_d[aite] > master_d[node] + dist){
				master_d[aite] = master_d[node] + dist;
				pq.push(MP(master_d[aite],aite));
			}
		}
	}
}

bool contain_rotate_poly(rotate_poly poly,point p){
	REP(i,4){
		if(contain_sector(poly.sec[i],p))return true;
	}
	if(contains(poly.p,p) == 2)return true;

	return false;
}

bool check_totate_poly(int num,point p){
	REP(i,r_poly[num].size()){
		if(contain_rotate_poly(r_poly[num][i],p))return true;
	}
	return false;
}

void make_rotate_graph(int num){

	//不可領域の作成
	double rad = num * PI / (double)r;
	double next_rad = (num + 1) * PI / (double)r;
	point vec = rotate(point(L, 0), rad);
	point next_vec = rotate(point(L, 0), next_rad);

	//回転不可領域作成
	REP(i,in_seg.size()){
		r_poly[num].PB(make_rotate_poly(vec,next_vec,in_seg[i],num));
	}

	REP(i, rotate_seg[num].size()) {
		for(int j = i+1;j <  rotate_seg[num].size();j++) {
			segment Aseg = rotate_seg[num][i];
			segment Bseg = rotate_seg[num][j];
			point tmp = crosspoint(Aseg,Bseg);
			if(tmp == point(INF,INF))continue;
			if(check_totate_poly(num,tmp))continue;
			check_visible(num,tmp,node_index,-1);
			check_visible((num+1)%r,tmp,node_index,1);
			node_index++;
		}
	}
}

int main() {
	//入力
	cin >> L >> r >> sx >> sy >> gx >> gy >> n;
	st = point(sx,sy);
	go = point(gx,gy);
	REP(i, n) {
		int a, b, c, d; cin >> a >> b >> c >> d;
		point A = point(a, b);
		point B = point(c, d);
		in_seg.PB(segment(A, B));
	}

	//各角度について、可視グラフを作成
	REP(i, r){
		make_visible_graph(i);
	}

	REP(i,r){
		make_rotate_graph(i);
	}

	//ダイクストラ
	Dijkstra();
	int ans = INF;
	REP(i,g_index.size()){
		ans = min(ans,master_d[g_index[i]]);
	}

	if(ans == INF){
		cout << -1 << endl;
	}
	else{
		cout << ans / 2 << endl;
	}

	return 0;
}
