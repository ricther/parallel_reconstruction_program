#include "CPoint.h"

int CPoint::index_counter=0;

int CPoint::get_index()
{
  return index_counter++;
}
CPoint::CPoint()
{
  x=0;y=0;z=0;index=0;used_times=0;
}
CPoint::CPoint(float inix,float iniy,float iniz)
{
  x=inix;
  y=iniy;
  z=iniz;
  index=index;
  used_times=0;
}

CPoint::CPoint(const CPoint& temp)
{
  x=temp.x;
  y=temp.y;
  z=temp.z;
  index=temp.index;
  used_times=temp.used_times;
}
void CPoint::operator=(CPoint temp)
{
  this->x=temp.x;
  this->y=temp.y;
  this->z=temp.z;
  this->index=temp.index;
  this-> used_times=temp.used_times;
}

void CPoint:: operator=(CPoint *temp)
{
  this->x=temp->x;
  this->y=temp->y;
  this->z=temp->z;
  this->index=temp->index;
  this->used_times=temp->used_times;
}

//void CPoint:: operator=(CPoint temp)
//{
////   x=temp.x;
//   y=temp.y;
//   z=temp.z;
// }
float CPoint::operator*(CPoint &right)
{
  float temp;
  
  temp+= this->x*right.x;
  temp+=this->y*right.y;
  temp+=this->z*right.z;
  return temp;
}

CPoint CPoint::operator-(CPoint& right)
{
  CPoint temp;
  temp.x=this->x-right.x;
  temp.y=this->y-right.y;
  temp.z= this->z-right.z;
  temp.index=this->index;
  return temp;
}
CPoint CPoint::operator+(CPoint& right)
{
  CPoint temp;
  temp.x=this->x+right.x;
  temp.y=this->y+right.y;
  temp.z=this->z+right.z;
  temp.index=this->index;
  return temp;
}

CPoint CPoint::operator/(float r)
{
  CPoint temp;
  temp.x=this->x/r;
  temp.y=this->y/r;
  temp.z=this->z/r;
  temp.index=this->index;
  return temp;
}
CPoint operator*(float &left,CPoint& right)
{
  CPoint temp;
  temp.x=left*right.x;
  temp.y=left*right.y;
  temp.z=left*right.z;
  temp.index=right.index;
  return temp;
}


