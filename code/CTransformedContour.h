#pragma once
#include "CContourDisplay.h"
#include <map>
#include "CLayer.h"
#include "CContour.h"
#include "CVirtualContour.h"
class CTransformedContour:public CContourDisplay
{
public:
CTransformedContour(std::map<float,CLayer*> & temp_map):CContourDisplay(temp_map){};
  virtual void set_iterator(CContour* temp,int deformID)
  {
    bool has_new_points=false;
    int size = temp->vec_virtual_contour.size();
    int i=0;
    for ( i = 0; i < size; ++i)
    {
      if (temp->vec_virtual_contour[deformID]->vec_new_points.size()==0)
      {
        deformID++;
      }
      else
      {
        has_new_points=true;
        point_itr= temp->vec_virtual_contour[deformID]->vec_new_points.begin();
        point_etr=  temp->vec_virtual_contour[deformID]->vec_new_points.end();
        point_size= temp->vec_virtual_contour[deformID]->vec_new_points.size();
        if (deformID+1<temp->vec_virtual_contour.size())
        {
          deformID++;
        }
      }
    }
    if (i>=size&&has_new_points==false)
    {
        point_itr= temp->vec_Points_Origin.begin();
        point_etr= temp->vec_Points_Origin.end();
        point_size= temp->vec_Points_Origin.size();
      
    }
  }
}
;
