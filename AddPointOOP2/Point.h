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

	Point() = default;		//ʹ��Ĭ�ϵĹ��캯��
	Point(const Point& point) = default;		//ʹ��Ĭ�ϵĿ������캯��
	Point(const double point_x, const double point_y, const double point_z) //һ�㹹�캯��
		:x(point_x), y(point_y), z(point_z)
		{
		}

	double Distance2(const Point& point) const {//
		//����֮������ƽ��
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
		//�㵽ԭ��ľ���ƽ��
		return x * x + y * y + z * z;
	}

	Point& operator= (const Point&) = default;	//ʹ��Ĭ�ϵĿ������캯��
	Point& operator= (Point&&) = default;
	bool operator == (const Point& point) const {
		//����==���ţ������ж��������Ƿ���ͬ
		return (this->x == point.x) && (this->y == point.y) && (this->z == point.z);
	}
};



#endif // !POINT


