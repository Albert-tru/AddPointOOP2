#pragma once


#include"Point.h"
#include<iostream>
#include<vector>

using namespace std;

class Block
{
public:
	vector<Point> _blockpoints;			//该block内的点集
	int row;		//该block所在的行数
	int col;		//该block所在的列数

};



