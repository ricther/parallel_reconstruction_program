#pragma once
#include "CContourDisplay.h"
#include <map>
#include "CLayer.h"
#include "CContour.h"
class CNormalContour:public CContourDisplay
{
public:
CNormalContour(std::map<float,CLayer*> & temp_map):CContourDisplay(temp_map){};
  virtual void set_iterator(CContour* temp)
  {
    /* point_itr= temp->vec_Points_Origin.begin(); */
    /* point_etr= temp->vec_Points_Origin.end(); */
    /* point_size= temp->vec_Points_Origin.size(); */
    point_itr= temp->vec_points.begin();
    point_etr= temp->vec_points.end();
    point_size= temp->vec_points.size();

  }
}
;
