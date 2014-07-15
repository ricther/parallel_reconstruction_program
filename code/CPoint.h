#pragma once
//athour xu liu 07/10/13
//the class of point
class CPoint
{
 public:
  CPoint(void);
  CPoint(float inix,float iniy,float iniz);
  CPoint(const CPoint& temp);
  float x;
  float y;
  float z;
  int index;
  void operator=(CPoint temp);
  void operator=(CPoint *temp);
  //  void operator=(CPoint temp);
  float operator*(CPoint& right);
  CPoint operator-(CPoint & right);
  CPoint operator+(CPoint & right);
  CPoint operator/(float r);
  static int index_counter;
  //carefully use this function, casue this will affect the index!!!!
  int get_index();
  int used_times;
};

CPoint operator*(float& left,CPoint& right);

