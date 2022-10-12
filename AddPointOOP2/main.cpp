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
#include <fstream>
#include <sstream>

using namespace std;

int main() {
	vector<Point>points;
	//Point p1(10, 50, 0);
	//Point p2(70, 8, 0);
	//Point p3(43, 52, 0);
	//Point p4(23, 66, 0);
	//Point p5(28, 5, 0);
	//p1.Num = 0;
	//p2.Num = 1;
	//p3.Num = 2;
	//p4.Num = 3;
	//points.push_back(p1);
	//points.push_back(p2);
	//points.push_back(p3);
	//points.push_back(p4);


	double arr1[24][3];
	ifstream infile;
	infile.open("归一化xy.txt");
	if (!infile) cout << "error" << endl;
	double t1, t2, t3, t4, t5, t6, t7;
	int i = 0;
	int j = 0;
	while (infile >> t1 >> t2) {
		cout << t1 << "         " << t2 << endl;
		arr1[i][0] = t1;
		arr1[i][1] = t2;
		arr1[i][2] = 0;
		j++;
		cout << arr1[i][0] << " " << arr1[i][1] << " " << arr1[i][2] << " " << j << endl;
		i++;

	}
	for (int i = 0; i < j; i++)
	{
		Point tp; tp.Num = i;
		tp.x = arr1[i][0], tp.y = arr1[i][1], tp.z = 0;
		points.insert(points.end(), tp);
		//ps1->InsertNextPoint(arr1[i][0], arr1[i][1], arr1[i][2]);
	}

	Delaunay delaunay;
	//delaunay.Init(points);
	//vector<Triangle> tris = delaunay.GetTriangles();

	//cout << "tris:  " << tris[0].a << "  " << tris[0].b << "  " << tris[0].c << endl;

	//Triangle tri = delaunay.PointLoacation(p, tris[0], 1);

	//cout << "main--tri: " << tri.a << "  " << tri.b << "  " << tri.c << endl;

	delaunay.Triangulation(points);

	//delaunay.Test(p3);
	//delaunay.Test(p2);
	//delaunay.Test(p3);
	//delaunay.Test(p4);
	//delaunay.Test(p5);


	vector<Triangle> tris = delaunay.GetTriangles();

	for (int i = 0; i < tris.size(); i++) {

		cout << "第" << i << "个三角形： " << tris[i].a << "  " << tris[i].b << "  " << tris[i].c << endl;
	}

	return 0;
}