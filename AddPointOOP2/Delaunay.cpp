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

using namespace std;

typedef pair<Edge, Triangle*> PLT;


bool Delaunay::almost_equal(const Point& p1, const Point& p2) const		//��
{

	return std::abs(p1.x - p2.x) < 0.0001f &&
		std::abs(p1.y - p2.y) < 0.00001f &&
		std::abs(p1.z - p2.z) < 0.00001f;
}

bool Delaunay::almost_equal(const Edge& e1, const Edge& e2)const		//��
{

	return	(almost_equal(_points[e1.p1], _points[e2.p1]) && almost_equal(_points[e1.p2], _points[e2.p2])) ||
		(almost_equal(_points[e1.p1], _points[e2.p2]) && almost_equal(_points[e1.p2], _points[e2.p1]));
}

bool Delaunay::almost_equal(const Triangle& t1, const Triangle& t2)const		//��
{

	return	(almost_equal(_points[t1.a], _points[t2.a]) || almost_equal(_points[t1.a], _points[t2.b]) || almost_equal(_points[t1.a], _points[t2.c])) &&
		(almost_equal(_points[t1.b], _points[t2.a]) || almost_equal(_points[t1.b], _points[t2.b]) || almost_equal(_points[t1.b], _points[t2.c])) &&
		(almost_equal(_points[t1.c], _points[t2.a]) || almost_equal(_points[t1.c], _points[t2.b]) || almost_equal(_points[t1.c], _points[t2.c]));
}

bool Delaunay::ContainVertex(const Triangle& t, const Point& p) {
	//�ж����������Ƿ���������
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
	//�ж�һ�����Ƿ�����������

	const double aa = _points[t.a].Norm2();		//��p1��ƽ����
	const double bb = _points[t.b].Norm2();
	const double cc = _points[t.c].Norm2();

	const double ax = _points[t.a].x;		//a���x����
	const double ay = _points[t.a].y;
	const double az = _points[t.a].z;
	const double bx = _points[t.b].x;
	const double by = _points[t.b].y;
	const double bz = _points[t.b].z;
	const double cx = _points[t.c].x;
	const double cy = _points[t.c].y;
	const double cz = _points[t.c].z;

	//��ȡԲ������
	const double circum_x = (aa * (cy - by) + bb * (ay - cy) + cc * (by - ay)) / (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
	const double circum_y = (aa * (cx - bx) + bb * (ax - cx) + cc * (bx - ax)) / (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));
	const double circum_z = (aa * (cz - bz) + bb * (az - cz) + cc * (bz - az)) / (ay * (cz - bz) + by * (az - cz) + cy * (bz - az));

	const Point circumCenter(circum_x / 2, circum_y / 2, circum_z / 2);		//Բ��
	const double circumRadius = _points[t.a].Distance2(circumCenter);		//a�㵽Բ�ĵľ���
	const double dist = p.Distance2(circumCenter);	//p�㵽Բ�ĵľ���

	return dist <= circumRadius;
}

double Delaunay::JudgePointLine(Point p1, Point p2, Point p) {
	double temp = (p2.x - p1.x) * (p.y - p1.y) - (p.x - p1.x) * (p2.y - p1.y);		//ʹ������������б���
	return temp;
}


void Delaunay::Init(const vector<Point>& points) {
	//���������ĵ㼯������_points��
	_points = points;

	//����һ�����ΰ�Χ��
	double minX = _points[0].x;
	double minY = _points[0].y;
	double minZ = _points[0].z;
	double maxX = minX;
	double maxY = minY;
	double maxZ = minZ;

	for (int i = 0; i < _points.size(); i++)		//��ѡ��x,y�е������Сֵ
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

	double delta = 5;		//��һ��С��������С���

	Point Vertex1(minX - delta, minY - delta, 0);		//���ΰ�Χ�еĶ��㣬������������ʱ��˳��
	Point Vertex2(minX - delta, maxY + delta, 0);
	Point Vertex3(maxX + delta, maxY + delta, 0);
	Point Vertex4(maxX + delta, minY - delta, 0);

	Vertex1.Num = _points.size();
	Vertex2.Num = _points.size() + 1;
	Vertex3.Num = _points.size() + 2;
	Vertex4.Num = _points.size() + 3;


	_points.push_back(Vertex1);		//���������εĶ�����ӵ�_points�У������ио����ַ��������ģ�����һʱ֮�����벻�������ķ���
	_points.push_back(Vertex2);
	_points.push_back(Vertex3);
	_points.push_back(Vertex4);

	Triangle tri1 = Triangle(_points.size() - 4, _points.size() - 3, _points.size() - 2);
	Triangle tri2 = Triangle(_points.size() - 1, _points.size() - 4, _points.size() - 2);

	_triangles.push_back(Triangle(_points.size() - 4, _points.size() - 3, _points.size() - 2));		//�Ѿ����ضԽ��߷ֳ��������������ηŽ�ȥ
	_triangles.push_back(Triangle(_points.size() - 1, _points.size() - 4, _points.size() - 2));		//ע�����������˳���

	//Edge e1(Vertex1.Num, Vertex2.Num);
	//Edge e2(Vertex2.Num, Vertex3.Num);
	//Edge e3(Vertex3.Num, Vertex1.Num);
	//Edge e4(Vertex3.Num, Vertex4.Num);
	//Edge e5(Vertex4.Num, Vertex1.Num);
	//Edge e6(Vertex1.Num, Vertex3.Num);

	LT.insert(PLT(tri1.edges[2], &tri2));
	LT.insert(PLT(tri2.edges[2], &tri1));

	Triangle* ptri = LT.at(tri1.edges[2]);
	cout<<" ptri:  " << ptri->a <<" "<<ptri->b <<" "<<ptri->c<< endl;

	Triangle* pttri = LT.at(tri2.edges[2]);
	cout << " pptri:  " << pttri->a << " " << pttri->b << " " << pttri->c << endl;

}

Triangle&  Delaunay::pointloacation(Point p, Triangle& tri, int x) {

	cout << "ccccc" << x << endl;

	cout << "tri:  " << tri.a << "  " << tri.b << "  " << tri.c << endl;

	Point p1 = _points[tri.a];
	Point p2 = _points[tri.b];
	Point p3 = _points[tri.c];

	int flag = 0;

	if (JudgePointLine(p1, p2, p) < 0) {	//λ�������ߵ��Ҳ�
		flag++;
	}
	else {									//λ�������ߵ����Ļ����ݹ�����������ⲿ���ڽ�������
		Triangle* pptr = LT.at(tri.edges[0]);
		cout << "pptr:  " << pptr->a << " " << pptr->b << " " << pptr->c;
		pointloacation(p, *(LT.at(tri.edges[0])), x + 1);
	}

	if (JudgePointLine(p2, p3, p) < 0) {
		flag++;
	}
	else {									//λ�������ߵ����Ļ����ݹ�����������ⲿ���ڽ�������
		pointloacation(p, *LT.at(tri.edges[1]), x + 1);
	}

	if (JudgePointLine(p3, p1, p) < 0) {
		flag++;
	}
	else {									//λ�������ߵ����Ļ����ݹ�����������ⲿ���ڽ�������
		Triangle* pptr = LT.at(tri.edges[2]);
		cout << "pptr:  " << pptr->a << " " << pptr->b << " " << pptr->c;
		pointloacation(p, *LT.at(tri.edges[2]), x + 1);
	}

	if (flag == 3)
		return tri;
}

const vector<Triangle>& Delaunay::triangulate(const vector<Point>& points) {
	//�Ը����ĵ㼯�������ǻ�
	//���������ĵ㼯������_points��
	_points = points;

	//����һ��������
	double minX = _points[0].x;
	double minY = _points[0].y;
	double minZ = _points[0].z;
	double maxX = minX;
	double maxY = minY;
	double maxZ = minZ;

	for (int i = 0; i < _points.size(); i++)		//��ѡ��x,y�е������Сֵ
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

	double delta = 5;		//��һ��С��������С���

	const Point Vertex1(minX - delta, minY - delta, 0);		//���ΰ�Χ�еĶ��㣬������������ʱ��˳��
	const Point Vertex2(minX - delta, maxY + delta, 0);
	const Point Vertex3(maxX + delta, maxY + delta, 0);
	const Point Vertex4(maxX + delta, minY - delta, 0);

	_points.push_back(Vertex1);		//���������εĶ�����ӵ�_points�У������ио����ַ��������ģ�����һʱ֮�����벻�������ķ���
	_points.push_back(Vertex2);
	_points.push_back(Vertex3);
	_points.push_back(Vertex4);

	_triangles.push_back(Triangle(_points.size() - 4, _points.size() - 3, _points.size() - 2));		//�Ѿ����ضԽ��߷ֳ��������������ηŽ�ȥ
	_triangles.push_back(Triangle(_points.size() - 1, _points.size() - 4, _points.size() - 2));		//ע�����������˳���

	//�����㼯�е����е�
	for (int i = 0; i < _points.size() - 4; i++) {

	}

	//�����㼯�����еĵ�
	for (size_t i = 0; i < _points.size() - 4; i++)
	{
		std::vector<Edge> polygon;		//�洢���������������εı�

		//�����������б��ж��������Ƿ���������
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

		//����������к��в���㣬������������б���ɾ��
		_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [](Triangle& t) {

			return t.isBad;
			}), end(_triangles));
		//�����߱�����Ƿ�����ͬ�ı�
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

		//�Ƴ���ͬ�ı�
		polygon.erase(std::remove_if(begin(polygon), end(polygon), [](Edge& e) {

			return e.isBad;
			}), end(polygon));

		//����µ�������
		for (const auto e : polygon)
			_triangles.push_back(Triangle(e.p1, e.p2, i));
	}

	//���ɾ�����а������������������������,���е�ÿһ����������һ��thisָ��
	//_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [p1, p2, p3, this](Triangle& t) {

	//	return ContainVertex(t, p1) || ContainVertex(t, p2) || ContainVertex(t, p3);
	//	}), end(_triangles));


	//��������εı���_edges�߱���
	for (const auto t : _triangles)
	{

		_edges.push_back(Edge{
			t.a, t.b });		//ʹ�ó�ʼ���б���ʵ��
		_edges.push_back(Edge{
			t.b, t.c });
		_edges.push_back(Edge{
			t.c, t.a });
	}

	//�����������еĵ㵯��
	_points.pop_back();
	_points.pop_back();
	_points.pop_back();

	return _triangles;
}

const vector<Triangle >& Delaunay::GetTriangles() const {

	return _triangles;
}
const vector<Edge >& Delaunay::GetEdges() const {

	return _edges;
}
const vector<Point >& Delaunay::GetVertices() const {

	return _points;
}



