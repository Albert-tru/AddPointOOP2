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

	vector<Triangle> _triangles;		//��������Ƭ����
	vector<Edge> _edges;		//�߼���
	vector<Point> _points;	//�㼯��

	//�жϵ㡢�ߡ�����ͬ�����غ���
	bool almost_equal(const Point& p1, const Point& p2) const;		//��

	bool almost_equal(const Edge& e1, const Edge& e2)const;		//��

	bool almost_equal(const Triangle& t1, const Triangle& t2)const;		//��

	bool ContainVertex(const Triangle& t, const Point& p);

	Point VectorCreated2D(Point a, Point b);

	Point VectorCreated3D(Point a, Point b);

	Point VectorCross2D(Point a, Point b);

	Point VectorCross3D(Point a, Point b);

	bool VectorSameDirection2D(Point a, Point b);

	bool VectorSameDirection3D(Point a, Point b);

	bool CircumCircleContainsVector(Triangle t, Point p);

	bool CircumCircleContains(const Triangle& t, const Point& p);

	double JudgePointLine(Point p1, Point p2, Point p);			//�ж�p��p1��p2���������or�Ҳ࣬>0��࣬<0�Ҳ�

public:
	map<Edge, Triangle*> LT;
	Delaunay() = default;
	Delaunay(const Delaunay&) = delete;		//��ʹ��Ĭ�ϵĿ������캯��
	Delaunay(Delaunay&&) = delete;


	//��ʼ�������մ�������points���������ΰ�Χ��
	//������ݣ�Vertexs1-4�ĸ�����
	void Init(const vector<Point>& points);

	//�㶨λ����tri�ݹ�����p���ڵ������Σ������ظ������ε����ã�
	Triangle& pointloacation(Point p, Triangle& tri, int x);

	const vector<Triangle>& triangulate(const vector<Point>&points);

	Delaunay& operator=(const Delaunay&) = delete;
	Delaunay& operator=(Delaunay&&) = delete;

	const vector<Triangle >& GetTriangles() const;
	const vector<Edge >& GetEdges() const;
	const vector<Point >& GetVertices() const;
};


#endif // !DELAUNARY_H