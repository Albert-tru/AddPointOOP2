#include<iostream>
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

int main() {
	vector<Point>points;
	Point p1(10, 100, 0);
	Point p2(100, 10, 0);
	Point p3(10, 10, 0);
	Point p4(100, 100, 0);
	Point p(40, 20, 0);
	p1.Num = 1;
	p2.Num = 2;
	p3.Num = 3;
	p4.Num = 4;
	points.push_back(p1);
	points.push_back(p2);
	points.push_back(p3);
	points.push_back(p4);

	Delaunay delaunay;
	delaunay.Init(points);
	vector<Triangle> tris = delaunay.GetTriangles();

	cout << "tris:  " << tris[0].a << "  " << tris[0].b << "  " << tris[0].c << endl;

	Triangle tri = delaunay.pointloacation(p, tris[0], 1);



	cout << "tri.size() = " << tris.size() << endl;

	return 0;
}