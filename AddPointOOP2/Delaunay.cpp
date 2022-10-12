#include "Delaunay.h"
#include"Point.h"
#include"Edge.h"
#include"Triangle.h"
#include"Point.h"
#include"Block.h"
#include"Blocks.h"

#include<vector>
#include<algorithm>
#include<iostream>
#include<map>
#include<set>

using namespace std;

typedef pair<Edge, Triangle> PLT;


bool Delaunay::almost_equal(const Point& p1, const Point& p2) const		//点
{

	return std::abs(p1.x - p2.x) < 0.0001f &&
		std::abs(p1.y - p2.y) < 0.00001f &&
		std::abs(p1.z - p2.z) < 0.00001f;
}

bool Delaunay::almost_equal(const Edge& e1, const Edge& e2)const		//线
{

	return	(almost_equal(_points[e1.p1], _points[e2.p1]) && almost_equal(_points[e1.p2], _points[e2.p2])) ||
		(almost_equal(_points[e1.p1], _points[e2.p2]) && almost_equal(_points[e1.p2], _points[e2.p1]));
}

bool Delaunay::almost_equal(const Triangle& t1, const Triangle& t2)const		//面
{

	return	(almost_equal(_points[t1.a], _points[t2.a]) || almost_equal(_points[t1.a], _points[t2.b]) || almost_equal(_points[t1.a], _points[t2.c])) &&
		(almost_equal(_points[t1.b], _points[t2.a]) || almost_equal(_points[t1.b], _points[t2.b]) || almost_equal(_points[t1.b], _points[t2.c])) &&
		(almost_equal(_points[t1.c], _points[t2.a]) || almost_equal(_points[t1.c], _points[t2.b]) || almost_equal(_points[t1.c], _points[t2.c]));
}

bool Delaunay::ContainVertex(const Triangle& t, const Point& p) {
	//判断三角形中是否包含这个点
	return almost_equal(_points[t.a], p) || almost_equal(_points[t.b], p) || almost_equal(_points[t.c], p);
}

Point Delaunay::VectorCreated2D(Point a, Point b) {
	Point t;
	t.x = a.x - b.x;
	t.y = a.y - b.y;
	t.z = 0;
	return t;
}

Point Delaunay::VectorCreated3D(Point a, Point b) {
	Point t;
	t.x = a.x - b.x;
	t.y = a.y - b.y;
	t.z = a.z - b.z;
	return t;
}

Point Delaunay::VectorCross2D(Point a, Point b) {
	Point tt;
	a.z = 0, b.z = 0;
	tt.x = a.y * b.z - b.y * a.z;
	tt.y = b.x * a.z - a.x * b.z;
	tt.z = a.x * b.y - b.x * a.y;
	return tt;
}

Point Delaunay::VectorCross3D(Point a, Point b) {
	Point tt;
	tt.x = a.y * b.z - b.y * a.z;
	tt.y = b.x * a.z - a.x * b.z;
	tt.z = a.x * b.y - b.x * a.y;
	return tt;
}

bool Delaunay::VectorSameDirection2D(Point a, Point b) {
	if ((a.z >= 0 && b.z >= 0)) {
		return true;
	}
	else if ((a.z < 0 && b.z < 0)) {
		return true;
	}
	else
		return false;
}

bool Delaunay::VectorSameDirection3D(Point a, Point b) {
	double flag1 = a.x * b.y;
	double flag2 = a.y * b.x;
	double flag3 = a.x * b.z;
	double flag4 = b.x * a.z;
	if (flag1 == flag2 && flag3 == flag4) {
		if (abs(a.x) != 0 && abs(b.x) != 0) {
			if (a.x * 1.0 / b.x >= 0.0000000000001)
				return true;
			else
				return false;
		}
		if (abs(a.y) != 0 && abs(b.y) != 0) {
			if (a.y * 1.0 / b.y >= 0.0000000000001)
				return true;
			else
				return false;
		}
		if (abs(a.z) != 0 && abs(b.z) != 0) {
			if (a.z * 1.0 / b.z >= 0.0000000000001)
				return true;
			else
				return false;
		}
		return true;
	}
}

bool Delaunay::CircumCircleContainsVector(Triangle t, Point p) {
	Point ab = VectorCreated2D(_points[t.a], _points[t.b]);
	Point ap = VectorCreated2D(_points[t.a], p);
	Point ac = VectorCreated2D(_points[t.a], _points[t.c]);
	Point bc = VectorCreated2D(_points[t.b], _points[t.c]);
	Point bp = VectorCreated2D(_points[t.b], p);
	Point ba = VectorCreated2D(_points[t.b], _points[t.a]);
	bool flag1 = false;
	bool flag2 = false;
	bool flag3 = false;
	if (VectorSameDirection2D(VectorCross2D(ab, ap), VectorCross2D(ab, ac)) &&
		VectorSameDirection2D(VectorCross2D(ac, ap), VectorCross2D(ac, ab)) &&
		VectorSameDirection2D(VectorCross2D(bc, bp), VectorCross2D(bc, ba))) {
		return true;
	}

}

bool Delaunay::CircumCircleContains(const Triangle& t, const Point& p) {
	//判断一个点是否在三角形内

	const double aa = _points[t.a].Norm2();		//点p1的平方和
	const double bb = _points[t.b].Norm2();
	const double cc = _points[t.c].Norm2();

	const double ax = _points[t.a].x;		//a点的x坐标
	const double ay = _points[t.a].y;
	//const double az = _points[t.a].z;
	const double bx = _points[t.b].x;
	const double by = _points[t.b].y;
	//const double bz = _points[t.b].z;
	const double cx = _points[t.c].x;
	const double cy = _points[t.c].y;
	//const double cz = _points[t.c].z;

	//求取圆心坐标
	const double circum_x = (aa * (cy - by) + bb * (ay - cy) + cc * (by - ay)) / (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
	const double circum_y = (aa * (cx - bx) + bb * (ax - cx) + cc * (bx - ax)) / (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));
	//const double circum_z = (aa * (cz - bz) + bb * (az - cz) + cc * (bz - az)) / (ay * (cz - bz) + by * (az - cz) + cy * (bz - az));

	const Point circumCenter(circum_x / 2, circum_y / 2, 0);		//圆心
	const double circumRadius = _points[t.a].Distance2(circumCenter);		//a点到圆心的距离
	const double dist = p.Distance2(circumCenter);	//p点到圆心的距离

	//cout << "dist:  " << dist << "   R:  " << circumRadius << endl;

	return dist <= circumRadius;
}

double Delaunay::JudgePointLine(Point p1, Point p2, Point p) {
	double temp = (p2.x - p1.x) * (p.y - p1.y) - (p.x - p1.x) * (p2.y - p1.y);		//使用向量叉乘来判别方向

	//cout << "JudgePointLine: " << temp << endl;
	return temp;
}

//三阶行列式
double Delaunay::Det3(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3) {
	return a1 * (b2 * c3 - c2 * b3) + a2 * (b3 * c1 - b1 * c3) + a3 * (b1 * c2 - c1 * b2);
}

//P是否在三角形外接圆内的判别准则,true表示圆外
bool Delaunay::JudgeCricumcircle(Point P, Triangle Tri) {
	double ax = _points[Tri.a].x;
	double ay = _points[Tri.a].y;
	double aa = pow(ax, 2) + pow(ay, 2);
	double bx = _points[Tri.b].x;
	double by = _points[Tri.b].y;
	double bb = pow(bx, 2) + pow(by, 2);
	double cx = _points[Tri.c].x;
	double cy = _points[Tri.c].y;
	double cc = pow(cx, 2) + pow(cy, 2);
	double px = P.x;
	double py = P.y;
	double pp = pow(px, 2) + pow(py, 2);
	double temp = ax * Det3(by, bb, 1, cy, cc, 1, px, pp, 1) - ay * Det3(bx, bb, 1, cx, cc, 1, px, pp, 1)
		+ aa * Det3(bx, by, 1, cx, cy, 1, px, py, 1) - Det3(bx, by, bb, cx, cy, cc, px, py, pp);

	cout << ax << " " << ay << " " << aa << " " << bx << " " << by << " "
		<< bb << " " << cx << " " << cy << " " << cc << " " << px << " " << py << " " << pp << " " << temp << endl;
	if (temp < 0)
		return true;
	else
		return false;

}

bool Delaunay::PointInTri(Point p, Triangle tri) {
	Point p1 = _points[tri.a];
	Point p2 = _points[tri.b];
	Point p3 = _points[tri.c];

	int flag = 0;

	if (JudgePointLine(p1, p2, p) <= 0) {	//位于这条边的右侧
		flag++;
	}
	if (JudgePointLine(p2, p3, p) <= 0) {
		flag++;
	}
	if (JudgePointLine(p3, p1, p) <= 0) {
		flag++;
	}
	if (flag == 3)
		return true;
	else
		return false;
}

void Delaunay::Init(const vector<Point>& points) {
	//将传进来的点集保存至_points中
	_points = points;

	//构建一个矩形包围盒
	double minX = _points[0].x;
	double minY = _points[0].y;
	double minZ = _points[0].z;
	double maxX = minX;
	double maxY = minY;
	double maxZ = minZ;

	for (int i = 0; i < _points.size(); i++)		//挑选出x,y中的最大最小值
	{

		if (_points[i].x < minX) minX = _points[i].x;
		if (_points[i].y < minY) minY = _points[i].y;
		if (_points[i].x > maxX) maxX = _points[i].x;
		if (_points[i].y > maxY) maxY = _points[i].y;
		//if (_points[i].z > maxZ) maxZ = _points[i].z;
		//if (_points[i].z > maxZ) maxZ = _points[i].z;
	}

	const double dx = maxX - minX;
	const double dy = maxY - minY;
	const double dz = maxZ - minZ;
	const double deltaMax = std::max(dx, dy);
	const double midx = (minX + maxX) / 2;
	const double midy = (minY + maxY) / 2;
	const double midz = (minZ + maxZ) / 2;

	double delta = 2;		//设一个小增量，减小误差

	Point Vertex1(minX - delta, minY - delta, 0);		//矩形包围盒的顶点，左下至右下逆时针顺序
	Point Vertex2(minX - delta, maxY + delta, 0);
	Point Vertex3(maxX + delta, maxY + delta, 0);
	Point Vertex4(maxX + delta, minY - delta, 0);

	Vertex1.Num = _points.size();
	Vertex2.Num = _points.size() + 1;
	Vertex3.Num = _points.size() + 2;
	Vertex4.Num = _points.size() + 3;


	_points.push_back(Vertex1);		
	_points.push_back(Vertex2);
	_points.push_back(Vertex3);
	_points.push_back(Vertex4);

	_vertexs.push_back(Vertex1);	//把这四个顶点加入vertexs
	_vertexs.push_back(Vertex2);
	_vertexs.push_back(Vertex3);
	_vertexs.push_back(Vertex4);

	Triangle tri1 = Triangle(_points.size() - 4, _points.size() - 3, _points.size() - 2);
	Triangle tri2 = Triangle(_points.size() - 2, _points.size() - 1, _points.size() - 4);

	_triangles.push_back(tri1);		//把矩形沿对角线分出来的两个三角形放进去
	_triangles.push_back(tri2);		//注意这个点是有顺序的

	//cout << "tri1: " << _triangles[0].a << "  " << _triangles[0].b << "  " << _triangles[0].c << endl;
	//cout << "tri2: " << _triangles[1].a << "  " << _triangles[1].b << "  " << _triangles[1].c << endl;

	//Edge e1(Vertex1.Num, Vertex2.Num);
	//Edge e2(Vertex2.Num, Vertex3.Num);
	//Edge e3(Vertex3.Num, Vertex1.Num);
	//Edge e4(Vertex3.Num, Vertex4.Num);
	//Edge e5(Vertex4.Num, Vertex1.Num);
	//Edge e6(Vertex1.Num, Vertex3.Num);

	LT.insert(PLT(tri1.edges[2], tri2));
	LT.insert(PLT(tri2.edges[2], tri1));

}

Triangle Delaunay::PointLoacation(Point p, Triangle tri, int x) {


	while (PointInTri(p, tri) == false) {
		Point p1 = _points[tri.a];
		Point p2 = _points[tri.b];
		Point p3 = _points[tri.c];

		int flag = 0;

		//cout << "tri:  " << tri.a << "  " << tri.b << "  " << tri.c << endl;

		if (JudgePointLine(p1, p2, p) > 0) {	//位于这条边的左侧
			tri = LT.at(tri.edges[0]);
			continue;
		}
		if (JudgePointLine(p2, p3, p) > 0) {
			tri = LT.at(tri.edges[1]);
			continue;
		}

		if (JudgePointLine(p3, p1, p) > 0) {
			tri = LT.at(tri.edges[2]);
			continue;
		}

	}

	return tri;
}

//void Delaunay::PointLoacation(Point p,  Triangle& tri, int x) {
//
//	cout << "ccccc" << x << endl;
//
//	cout << "tri:  " << tri.a << "  " << tri.b << "  " << tri.c << endl;
//
//	Point p1 = _points[tri.a];
//	Point p2 = _points[tri.b];
//	Point p3 = _points[tri.c];
//
//	int flag = 0;
//
//	if (JudgePointLine(p1, p2, p) <= 0) {	//位于这条边的右侧
//		flag++;
//	}
//	else {									//位于这条边的左侧的话，递归查找这条边外部的邻接三角形
//
//		PointLoacation(p, LT.at(tri.edges[0]), x + 1);
//	}
//
//	if (JudgePointLine(p2, p3, p) <= 0) {
//		flag++;
//	}
//	else {									//位于这条边的左侧的话，递归查找这条边外部的邻接三角形
//		PointLoacation(p, LT.at(tri.edges[1]), x + 1);
//	}
//
//	if (JudgePointLine(p3, p1, p) <= 0) {
//		flag++;
//	}
//	else {									//位于这条边的左侧的话，递归查找这条边外部的邻接三角形
//		PointLoacation(p, LT.at(tri.edges[2]), x + 1);
//	}
//
//	if (flag == 3)
//		return;
//
//}

//void Delaunay::JudgeLine(Edge ln, Point P, vector<Edge>& Poly) {
//
//	if (LT.find(ln) == LT.end())
//		return;
//	else {
//		Triangle tri = LT.at(ln);////////////////////////$
//
//		if (JudgeCricumcircle(P, tri) == true) {    //在圆内或圆上
//
//			//找其他两边
//			for (int i = 0; i <= 2; i++) {
//				Point P1, P2;
//				if (i == 0) {
//					P1 = _points[tri.a];
//					P2 = _points[tri.b];
//				}
//				else if (i == 1) {
//					P1 = _points[tri.b];
//					P2 = _points[tri.c];
//				}
//				else {
//					P1 = _points[tri.c];
//					P2 = _points[tri.a];
//				}
//
//				double temp = (P2.x - P1.x) * (P.y - P1.y) - (P.x - P1.x) * (P2.y - P1.y);
//				if (temp > 0) {//如果插入点位于这条边左侧的话
//
//					//删除这条边,添加其他两条边
//					int la = i + 1, lb = i + 2;
//					if (la >= 3)
//						la -= 3;
//					if (lb >= 3)
//						lb -= 3;
//					Poly.LineN.insert(Poly.LineN.begin() + i, tri.Line3[la]);
//					Poly.LineN.insert(Poly.LineN.begin() + i + 1, tri.Line3[lb]);
//					Poly.LineN.erase(Poly.LineN.begin() + i);
//
//					//接着搜索其他两条边
//					JudgeLine(tri.Line3[la], P, Poly);
//					JudgeLine(tri.Line3[lb], P, Poly);
//				}
//
//
//			}
//		}
//
//		else {  //在外接圆外
//			return;
//		}
//	}
//
//
//}

void Delaunay::DeleteTriangleLT(Triangle tri) {

	vector<Edge>e;

	Edge e1(tri.a, tri.b);
	Edge e2(tri.b, tri.c);
	Edge e3(tri.c, tri.a);

	Edge e4(tri.b, tri.a);
	Edge e5(tri.c, tri.b);
	Edge e6(tri.a, tri.c);

	e.push_back(e1), e.push_back(e2), e.push_back(e3);
	e.push_back(e4), e.push_back(e5), e.push_back(e6);

	for (int i = 0; i < e.size(); i++) {
		if (LT.find(e[i]) != LT.end()) {
			LT.erase(e[i]);
		}
	}

}

void Delaunay::CreatePoly(Point p, vector<Edge>& Poly) {

	for (int i = 0; i < Poly.size(); i++) {
		if (LT.find(Poly[i]) != LT.end()) {			//如果这条边有邻接三角形，才考虑是否在其外接圆内

			Triangle tri = LT.at(Poly[i]);		//tri表示这条边外部的邻接三角形

			if (CircumCircleContains(tri, p) == false) {	//p位于圆外,需要覆盖指向这条边内的那对键值

				if (LT.find(Poly[i]) != LT.end()) {		//如果这条边外部有三角形且未加入map的话，map中插入这对键值

					//LT.insert(PLT(Poly[i], Triangle(p.Num, Poly[i].p2, Poly[i].p1)));
					LT[Edge(Poly[i].p2, Poly[i].p1)] = Triangle(p.Num, Poly[i].p1, Poly[i].p2);

				}

				continue;
			}
			else {												//p位于圆内

				Triangle tt(tri.a, tri.b, tri.c);

				//需要删除这个邻接三角形，并在map中去除这对键值(注意键值对的删除，这个三角形的每条边的两对键值都需要删，因为边去掉后，两边的邻接三角形都应该删除
				vector<Triangle>::iterator iter = find(_triangles.begin(), _triangles.end(), tt);
				_triangles.erase(iter);

				LT.erase(Poly[i]);
				LT.erase(Edge(Poly[i].p2, Poly[i].p1));

				//先插入新边，再删除原边，注意顺序

				int order = 0;
				//找新边,并添加


				int bgid = Poly[i].p1;		//存储第一个顶点的位置
				for (int k = 0; k < 3; k++) {
					if (tri.a == bgid)
						bgid = k;
					if (tri.b == bgid)
						bgid = k;
					if (tri.c == bgid)
						bgid = k;
				}

				//添加新边
				int j = bgid;
				int num = 3;
				for (int j = 0; j < 3; j++) {
					Point P1, P2;
					int idP1, idP2;
					if (j == 0) {
						P1 = _points[tri.a];
						P2 = _points[tri.b];
					}
					else if (j == 1) {
						P1 = _points[tri.b];
						P2 = _points[tri.c];

					}
					else {
						P1 = _points[tri.c];
						P2 = _points[tri.a];
					}

					idP1 = P1.Num;
					idP2 = P2.Num;

					if (JudgePointLine(P1, P2, p) < 0) {
						order++;
						Edge e(idP1, idP2);
						Poly.insert(Poly.begin() + i + order, e);		//这里要添加两次///////记得检查一下插入位置是否正确
						if (order == 2)
							break;
					}
				}

				//找旧边，并删除
				num = 3;
				j = bgid;
				for (int j = 0; j < 3; j++) {
					Point P1, P2;
					int idP1, idP2;
					if (j == 0) {
						P1 = _points[tri.a];
						P2 = _points[tri.b];
					}
					else if (j == 1) {
						P1 = _points[tri.b];
						P2 = _points[tri.c];

					}
					else {
						P1 = _points[tri.c];
						P2 = _points[tri.a];
					}

					idP1 = P1.Num;
					idP2 = P2.Num;

					if (JudgePointLine(P1, P2, p) > 0) {
						Edge e(idP1, idP2);
						Poly.erase(Poly.begin() + i);
						break;
					}
				}

				//使下轮循环求的是插入的新边						
				i--;
			}
		}
		else {					//无邻接三角形，不用管
			continue;
		}
	}

	vector<Edge> Polyf;
	Polyf.push_back(Poly[0]);

	//给Poly排个序，使Poly首尾相接
	for (int i = 1; i < Poly.size(); i++) {
		int t = Polyf[i-1].p2;
		int j,id1,id2;
		for (j = 0; j < Poly.size(); j++) {
			if (Poly[j].p1 == t)
			{
				t = j;
				break;
			}
				
		}
		Polyf.push_back(Poly[t]);
	}
	for (int i = 0; i < Polyf.size(); i++) {
		Poly[i].p1 = Polyf[i].p1;
		Poly[i].p2 = Polyf[i].p2;
	}

	cout << "该点的影响域Poly创建完成！" << endl;
	for (int i = 0; i < Poly.size(); i++) {
		cout << "    " << Poly[i].p1 << "     " << Poly[i].p2 << endl;
	}

}

void Delaunay::LocalConstruction(Point p) {

	//点定位时初始三角形的设置											//////////点定位时初始三角形的设置可优化
	Triangle Tri = _triangles.front();

	//以包含该点的三角形为初始多边形,把该三角形内的顶点和相邻三角形按顺序放到Poly中
	Triangle InitTri = PointLoacation(p, Tri, 1);

	cout << "该点所在的三角形已找到： " << InitTri.a << "  " << InitTri.b << "  " << InitTri.c << endl;

	//需要删除这个作为初始多边形的三角形								
	vector<Triangle>::iterator iter = find(_triangles.begin(), _triangles.end(), Triangle(InitTri.a, InitTri.b, InitTri.c));
	_triangles.erase(iter);

	vector<Edge> Poly;
	for (int i = 0; i <= 2; i++) {
		// Poly.PolyVertexs.insert(Poly.PolyVertexs.end(),InitTri.Vertex3[i]);
		// Poly.PolyTriangles.insert(Poly.PolyTriangles.end(),InitTri.Triangle3[i]);
		Poly.insert(Poly.end(), InitTri.edges[i]);

	}

	//确定p的影响域
	CreatePoly(p, Poly);

	//连接p和Poly中的所有顶点，把形成的三角形放入三角网
	int tlen = _triangles.size();
	for (int i = 0; i < Poly.size(); i++) {
		int idp1 = Poly[i].p1;
		int idp2 = Poly[i].p2;
		Triangle tri(p.Num, idp1, idp2);
		_triangles.insert(_triangles.begin(), tri);
	}

	//确认新加入的三角形之间的邻接关系											
	for (int i = 0; i < _triangles.size() - tlen; i++) {
		if (i < _triangles.size() - tlen - 1) {
			Edge e(_triangles[i].a, _triangles[i].b);
			LT.insert(PLT(e, _triangles[i + 1]));
		}
		if (i >= 1) {
			Edge e1(_triangles[i].c, _triangles[i].a);
			LT.insert(PLT(e1, _triangles[i - 1]));
		}
	}
	int len = _triangles.size();
	Edge e(_triangles[len - tlen - 1].a, _triangles[len - tlen - 1].b);
	LT.insert(PLT(e, _triangles[0]));

	Edge e1(_triangles[0].c, _triangles[0].a);
	LT.insert(PLT(e1, _triangles[len - tlen - 1]));

}

void Delaunay::Triangulation(const vector<Point>& points) {			

	int length = points.size();
	Init(points);
	for (int i = 0; i < points.size(); i++) {
		Point p = _points[i];
		LocalConstruction(p);
		cout << "第" << i + 1 << "个点  " << _points[i].x << " " << _points[i].y << " " << "的插入结束----------------------------------：" << endl;
		ShowTriangles();
		CreateEdges();
		Vertify(i);
	}

	//删除和矩形包围盒相关的三角形
	// 
	//最后删除所有包含超三角三个顶点的三角形,类中的每一个函数都有一个this指针
	Point p1 = _vertexs[0], p2 = _vertexs[1], p3 = _vertexs[2], p4 = _vertexs[3];
	_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [p1, p2, p3,p4, this](Triangle& t) {

	return ContainVertex(t, p1) || ContainVertex(t, p2) || ContainVertex(t, p3) || ContainVertex(t,p4);
	}), end(_triangles));
}

void Delaunay::CreateEdges() {
	set<Edge>edges;
	for (int i = 0; i < _triangles.size(); i++) {
		Point P1, P2;
		int idP1, idP2;
		Triangle tri = _triangles[i];
		for (int j = 0; j <= 2; j++) {
			if (j == 0) {
				P1 = _points[tri.a];
				P2 = _points[tri.b];
			}
			else if (j == 1) {
				P1 = _points[tri.b];
				P2 = _points[tri.c];

			}
			else {
				P1 = _points[tri.c];
				P2 = _points[tri.a];
			}

			idP1 = P1.Num;
			idP2 = P2.Num;
			edges.insert(Edge(idP1, idP2));
		}
	}

	int edgesNum = (edges.size() - 4) / 2;

	cout << "edges:      " << edges.size() <<"  " << edgesNum <<endl;

	if (edgesNum != LT.size() / 2) {
		cout << "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ" << endl;
	}

}


void Delaunay::Vertify(int num) {

	for (int i = 0; i < _triangles.size(); i++) {
		int p1 = _triangles[i].a, p2 = _triangles[i].b, p3 = _triangles[i].c;
		for (int j = 0; j <= num; j++) {
			int p = _points[j].Num;
			
			if (p != p1 && p != p2 && p != p3) {
				if (CircumCircleContainsVector(_triangles[i], _points[j]) == true) {
					cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" << endl;
					cout << "num: " << num << "  Point: " << _points[j].Num << "  Tri: " << i << endl;
				}
			}
			
		}
	}

}

const vector<Triangle>& Delaunay::triangulate(const vector<Point>& points) {
	//对给出的点集进行三角化
	//将传进来的点集保存至_points中
	_points = points;

	//构建一个超三角
	double minX = _points[0].x;
	double minY = _points[0].y;
	double minZ = _points[0].z;
	double maxX = minX;
	double maxY = minY;
	double maxZ = minZ;

	for (int i = 0; i < _points.size(); i++)		//挑选出x,y中的最大最小值
	{

		if (_points[i].x < minX) minX = _points[i].x;
		if (_points[i].y < minY) minY = _points[i].y;
		if (_points[i].x > maxX) maxX = _points[i].x;
		if (_points[i].y > maxY) maxY = _points[i].y;
		//if (_points[i].z > maxZ) maxZ = _points[i].z;
		//if (_points[i].z > maxZ) maxZ = _points[i].z;
	}

	const double dx = maxX - minX;
	const double dy = maxY - minY;
	const double dz = maxZ - minZ;
	const double deltaMax = std::max(dx, dy);
	const double midx = (minX + maxX) / 2;
	const double midy = (minY + maxY) / 2;
	const double midz = (minZ + maxZ) / 2;

	double delta = 5;		//设一个小增量，减小误差

	const Point Vertex1(minX - delta, minY - delta, 0);		//矩形包围盒的顶点，左下至右下逆时针顺序
	const Point Vertex2(minX - delta, maxY + delta, 0);
	const Point Vertex3(maxX + delta, maxY + delta, 0);
	const Point Vertex4(maxX + delta, minY - delta, 0);

	_points.push_back(Vertex1);		//将超三角形的顶点添加到_points中，但是有感觉这种方法笨笨的，不过一时之间又想不到其他的方法
	_points.push_back(Vertex2);
	_points.push_back(Vertex3);
	_points.push_back(Vertex4);

	_triangles.push_back(Triangle(_points.size() - 4, _points.size() - 3, _points.size() - 2));		//把矩形沿对角线分出来的两个三角形放进去
	_triangles.push_back(Triangle(_points.size() - 1, _points.size() - 4, _points.size() - 2));		//注意这个点是有顺序的

	//遍历点集中的所有点
	for (int i = 0; i < _points.size() - 4; i++) {

	}

	//遍历点集中所有的点
	for (size_t i = 0; i < _points.size() - 4; i++)
	{
		std::vector<Edge> polygon;		//存储包含插入点的三角形的边

		//遍历三角形列表，判断三角形是否包含插入点
		for (auto& t : _triangles)
		{
			if (CircumCircleContainsVector(t, _points[i]))
			{
				t.isBad = true;

				polygon.push_back(Edge{
					t.a, t.b });
				polygon.push_back(Edge{
					t.b, t.c });
				polygon.push_back(Edge{
					t.c, t.a });
			}
		}

		//如果三角形中含有插入点，则将其从三角形列表中删除
		_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [](Triangle& t) {

			return t.isBad;
			}), end(_triangles));
		//遍历边表，检查是否有相同的边
		for (auto e1 = begin(polygon); e1 != end(polygon); ++e1)
		{

			for (auto e2 = e1 + 1; e2 != end(polygon); ++e2)
			{

				if (almost_equal(*e1, *e2))
				{

					e1->isBad = true;
					e2->isBad = true;
				}
			}
		}

		//移除相同的边
		polygon.erase(std::remove_if(begin(polygon), end(polygon), [](Edge& e) {

			return e.isBad;
			}), end(polygon));

		//添加新的三角形
		for (const auto e : polygon)
			_triangles.push_back(Triangle(e.p1, e.p2, i));
	}

	//最后删除所有包含超三角三个顶点的三角形,类中的每一个函数都有一个this指针
	//_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [p1, p2, p3, this](Triangle& t) {

	//	return ContainVertex(t, p1) || ContainVertex(t, p2) || ContainVertex(t, p3);
	//	}), end(_triangles));


	//添加三角形的边至_edges边表中
	for (const auto t : _triangles)
	{

		//_edges.push_back(Edge{
		//	t.a, t.b });		//使用初始化列表来实现
		//_edges.push_back(Edge{
		//	t.b, t.c });
		//_edges.push_back(Edge{
		//	t.c, t.a });
	}

	//将超三角形中的点弹出
	_points.pop_back();
	_points.pop_back();
	_points.pop_back();

	return _triangles;
}

void Delaunay::Test(Point p) {
	Triangle tri(0, 1, 6);
	//JudgeCricumcircle(p, tri);
	cout << p.x << "   " << p.y << endl;
	if (CircumCircleContains(tri, p) == true)
		cout << "true" << endl;
	else
		cout << "false" << endl;
}

void Delaunay::ShowTriangles() {
	cout << "ShwoTriangles:" << endl;
	for (int i = 0; i < _triangles.size(); i++) {
		
		cout << "第" << i << "个三角形： " << _triangles[i].a << "  " << _triangles[i].b << "  " << _triangles[i].c << endl;
	}
}

const vector<Triangle >& Delaunay::GetTriangles() const {

	return _triangles;
}
const set<Edge >& Delaunay::GetEdges() const {

	return _edges;
}
const vector<Point >& Delaunay::GetVertices() const {

	return _points;
}



