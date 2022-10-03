#pragma once
#ifndef DELAUNARY_H
#define DELAUNARY_H

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

class Delaunay
{

	vector<Triangle> _triangles;		//三角形面片集合
	vector<Edge> _edges;		//边集合
	vector<Point> _points;	//点集合

	//判断点、线、面相同的重载函数
	bool almost_equal(const Point& p1, const Point& p2) const;		//点

	bool almost_equal(const Edge& e1, const Edge& e2)const;		//线

	bool almost_equal(const Triangle& t1, const Triangle& t2)const;		//面

	bool ContainVertex(const Triangle& t, const Point& p);

	Point VectorCreated2D(Point a, Point b);

	Point VectorCreated3D(Point a, Point b);

	Point VectorCross2D(Point a, Point b);

	Point VectorCross3D(Point a, Point b);

	bool VectorSameDirection2D(Point a, Point b);

	bool VectorSameDirection3D(Point a, Point b);

	bool CircumCircleContainsVector(Triangle t, Point p);

	bool CircumCircleContains(const Triangle& t, const Point& p);

	double JudgePointLine(Point p1, Point p2, Point p);			//判断p在p1，p2向量的左侧or右侧，>0左侧，<0右侧

public:
	map<Edge, Triangle*> LT;
	Delaunay() = default;
	Delaunay(const Delaunay&) = delete;		//不使用默认的拷贝构造函数
	Delaunay(Delaunay&&) = delete;


	//初始化，接收传进来的points，建立矩形包围盒
	//相关数据：Vertexs1-4四个顶点
	void Init(const vector<Point>& points);

	//点定位，从tri递归搜索p所在的三角形，并返回该三角形的引用？
	Triangle& pointloacation(Point p, Triangle& tri, int x);

	const vector<Triangle>& triangulate(const vector<Point>&points);

	Delaunay& operator=(const Delaunay&) = delete;
	Delaunay& operator=(Delaunay&&) = delete;

	const vector<Triangle >& GetTriangles() const;
	const vector<Edge >& GetEdges() const;
	const vector<Point >& GetVertices() const;
};


#endif // !DELAUNARY_H