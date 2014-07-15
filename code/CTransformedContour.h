#pragma once
#include "CContourDisplay.h"
#include <map>
#include "CLayer.h"
#include "CContour.h"

class CTransformedContour:public CContourDisplay
{
public:
CTransformedContour(std::map<float,CLayer*> & temp_map):CContourDisplay(temp_map){};
  virtual void set_iterator(CContour* temp)
  {
    point_itr= temp->vec_new_points.begin();
    point_etr= temp->vec_new_points.end();
    point_size= temp->vec_new_points.size();
    if (point_size==0)
    {
      point_itr= temp->vec_Points_Origin.begin();
      point_etr= temp->vec_Points_Origin.end();
      point_size= temp->vec_Points_Origin.size();
    }
  }
}
;
