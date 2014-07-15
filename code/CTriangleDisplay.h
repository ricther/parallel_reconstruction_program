#pragma once
#include <vtkTriangle.h>
#include "CCorespond.h"
#include <map>
#include "CTriangleSetDisplay.h"
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>


class CShape;
class CTriangleDisplay
{
public:
  CTriangleDisplay();
  std::multimap<float, CTriangleSetDisplay*> map_layer_triangles_actors;
  void initial_actors(CShape*);
  void update_actors(float layerID,int lineID);
  float margin;

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
   vtkSmartPointer<vtkActor> m_actor;
   vtkSmartPointer<vtkActor> m_mesh_actor;
   vtkSmartPointer<vtkPolyData> m_polydata;
   vtkSmartPointer<vtkPolyDataMapper> m_polydata_mapper;
   vtkSmartPointer<vtkPoints>m_points;
   vtkSmartPointer<vtkCellArray>m_triangles;

   void set_up_vtk();
   int get_id(int ,int ,int);
   void organize_data(CCorrespond*);
   void set_up_data(CShape*);
   void initial_actor(CShape *);
   void smooth();
   private:
   int count;
   std::map<int,int> map_pindex_vtkindex;
};
