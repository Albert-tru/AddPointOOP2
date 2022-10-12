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
#include<set>

using namespace std;

class Delaunay
{
	vector<Point> _vertexs;			//���ΰ�Χ���ĸ�����
	vector<Triangle> _triangles;		//��������Ƭ����
	set<Edge> _edges;		//�߼���
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

	bool CircumCircleContains(const Triangle& t, const Point& p);	//true��ʾԲ��

	double JudgePointLine(Point p1, Point p2, Point p);			//�ж�p��p1��p2���������or�Ҳ࣬>0��࣬<0�Ҳ�

	//��������ʽ
	double Det3(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3);

	//P�Ƿ������������Բ�ڵ��б�׼��,true��ʾԲ��
	bool JudgeCricumcircle(Point P, Triangle Tri);

	bool PointInTri(Point p, Triangle tri);


public:
	map<Edge, Triangle> LT;
	Delaunay() = default;
	Delaunay(const Delaunay&) = delete;		//��ʹ��Ĭ�ϵĿ������캯��
	Delaunay(Delaunay&&) = delete;


	//һ����ʼ�������մ�������points���������ΰ�Χ��
	//	������ݣ�Vertexs1-4�ĸ�����
	void Init(const vector<Point>& points);

	//�������������--�ֿ��㷨
																				///////�ֿ��㷨��ûʵ��
	void PointsResorted(){}

	//�����㶨λ����tri�ݹ�����p���ڵ������Σ������ظ������ε����ã�
	Triangle PointLoacation(Point p, Triangle tri, int x);

	//�ġ�ȷ����P��Ӱ����Poly�������оֲ��ع�

	//4.1 ȷ��P��Ӱ����

	void DeleteTriangleLT(Triangle tri);

	void CreatePoly(Point p, vector<Edge>& Poly);
	
	//4.2 ��ʽ�ֲ��ع�
	void LocalConstruction(Point P);

	//�塢��ʽ���������ʷ�
	void Triangulation(const vector<Point>& points);

	const vector<Triangle>& triangulate(const vector<Point>&points);

	void CreateEdges();
	void Test(Point p);
	void Vertify(int);		//У�飬�ж����������Ƿ���������
	void ShowTriangles();

	Delaunay& operator=(const Delaunay&) = delete;
	Delaunay& operator=(Delaunay&&) = delete;

	const vector<Triangle >& GetTriangles() const;
	const set<Edge >& GetEdges() const;
	const vector<Point >& GetVertices() const;
};


#endif // !DELAUNARY_H