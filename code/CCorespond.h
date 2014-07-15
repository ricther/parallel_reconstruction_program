#pragma once
//author by xu liu 07/10/13
//the class of the correspondlayer
#include <vector>
#include <map>
#include "CPoint.h"
typedef struct
{
  CPoint point1;//lower layer's point
  CPoint point2;//higher layer's point
}Parapoint;
class CTriangle;
class CCorrespond
{
public:
 float  CorrespondLayer1;
 float  CorrespondLayer2;
 std:: map<int,Parapoint*> map_cor;//the parameter one is the seqID
 std:: vector<CTriangle*> vec_Triangle;
 void produce_triangle();
};
