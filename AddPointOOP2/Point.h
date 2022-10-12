#pragma once
#ifndef POINT
#define POINT

#include <math.h>


class Point
{

public:
	int Num;
	double x;
	double y;
	double z;

	Point() = default;		//使用默认的构造函数
	Point(const Point& point) = default;		//使用默认的拷贝构造函数
	Point(const double point_x, const double point_y, const double point_z) //一般构造函数
		:x(point_x), y(point_y), z(point_z)
		{
		}

	double Distance2(const Point& point) const {//
		//两点之间距离的平方
		const double dx = this->x - point.x;
		const double dy = this->y - point.y;
		//const double dz = this->z - point.z;
		//return dx * dx + dy * dy + dz * dz;
		return dx * dx + dy * dy;
	}

	double Distance(const Point& point)const {

		return sqrt(Distance2(point));
	}

	double Norm2() const {
		//点到原点的距离平方
		return x * x + y * y + z * z;
	}

	Point& operator= (const Point&) = default;	//使用默认的拷贝构造函数
	Point& operator= (Point&&) = default;
	bool operator == (const Point& point) const {
		//重载==符号，用于判断两个点是否相同
		return (this->x == point.x) && (this->y == point.y) && (this->z == point.z);
	}
};



#endif // !POINT


