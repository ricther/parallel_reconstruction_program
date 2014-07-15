#pragma once
#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkLineSource.h"
#include <vector>
#include "CCorespond.h"
#include "vtkProperty.h"
#include "vtkCellArray.h"

extern const int interval_of_line_layer;
class CLineSetDisplay
{
public:
  CLineSetDisplay();
  
  float LayerID;
  std::vector< vtkSmartPointer<vtkActor> > vec_actor;
  std::vector< vtkSmartPointer<vtkLineSource> > vec_line_source;

  vtkSmartPointer<vtkActor> m_actor;
  vtkSmartPointer<vtkPolyData> m_polydata;
  vtkSmartPointer<vtkPolyDataMapper> m_polydata_mapper;
  vtkSmartPointer<vtkPoints>m_points;
  vtkSmartPointer<vtkCellArray>m_lines;

  
  
  void set_actor_visibility(bool b_swith,int& lineID)
  {
    std::vector<vtkSmartPointer<vtkActor> >::iterator itr,etr;
    itr = vec_actor.begin();etr=vec_actor.end();
    int size=vec_actor.size();
    if (lineID<0)
    {
      lineID=0;
    }
    else if (lineID>size-1)
    {
      lineID=size-1;
    }
    int count=0;
    for (; itr!=etr; ++itr)
    {
      if (count==lineID)
      {
        (*itr)->GetProperty()->SetColor(255,0,0);
        std::cout<<"current line ID:"<<lineID<<std::endl;
      }
      else
      {
        (*itr)->GetProperty()->SetColor(255,255,255);
      }
      (*itr)->SetVisibility(b_swith);
      (*itr)->Modified();
      count++;
    }
  }
  void initial_actors(CCorrespond*);
  void initial_actor(CCorrespond*);
}
;







