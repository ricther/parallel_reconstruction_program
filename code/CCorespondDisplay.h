#pragma once
#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkLineSource.h"
#include "CCorespond.h"
#include <map>
#include "CLineSetDisplay.h"

class CShape;
class CCorespondDisplay
{
public:
  CCorespondDisplay()
  {
    float n=999999;
    max_z=-n;min_z=n;
    margin=0;
    draw_by_layer=true;
  }
  std::multimap<float, CLineSetDisplay*> map_line_actors;
  std::map<int, CLineSetDisplay*> map_int_line_actors;
  void initial_actors(CShape*);
  void initial_actor(CShape*);
  void update_actors(float layerID,int& lineID);
  void solo_corespond_display(int contourID);
  float margin;
  bool draw_by_layer;
   float max_z,min_z;

   void check_level(float level)
   {
     if (level>max_z)
     {
       max_z=level;
     }
     if (level<min_z)
     {
       min_z=level;
     }
   }

   bool check_display(float id,float level)
   {
      float up=id+margin;
      float down=id;
      if (up>=max_z)
      {
        up=max_z-1;
      }
      if(down<min_z)
      {
        down=min_z;
      }
      if (level<=up&&level>=down)
      {
        return true;
      }
      else
      {
        return false;
      }
   }


};
