#pragma once
#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkTriangle.h"
#include <vector>
#include "CCorespond.h"
#include "vtkProperty.h"
#include "CBaseDisplay.h"
class CTriangleSetDisplay:public CBaseDisplay
{
public:
  CTriangleSetDisplay();
  vtkSmartPointer<vtkActor> m_mesh_actor;
  float LayerID;
  void set_actor_visibility(bool b_swith,int& levelID)
  {
    /* std::vector<vtkSmartPointer<vtkActor> >::iterator itr,etr; */
    /* itr = vec_actor.begin();etr=vec_actor.end(); */
    /* int size=vec_actor.size(); */
    /* if (levelID<0) */
    /* { */
    /*   levelID=0; */
    /* } */
    /* else if (levelID>size-1) */
    /* { */
    /*   levelID=size-1; */
    /* } */
    /* int count=0; */
    /* for (; itr!=etr; ++itr) */
    /* { */
    /*   if (count==levelID) */
    /*   { */
    /*     (*itr)->GetProperty()->SetColor(255,0,0); */
    /*     std::cout<<"current line ID:"<<levelID<<std::endl; */
    /*   } */
    /*   else */
    /*   { */
    /*     (*itr)->GetProperty()->SetColor(255,255,255); */
    /*   } */
    /*   (*itr)->SetVisibility(b_swith); */
    /*   (*itr)->Modified(); */
    /*   count++; */
    /* } */
  }
  void initial_actors(CCorrespond*);
}
;

