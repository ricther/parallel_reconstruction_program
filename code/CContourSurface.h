#pragma once
#include <vtkTriangle.h>
#include "CCorespond.h"
#include <map>
#include <vector>
#include "CTriangleSetDisplay.h"
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>

class CShape;
class CContourSurface
{
public:
  CContourSurface();
  void initial_actors(CShape*);
  void set_up_data_Delaunay2D(CShape *m_shape);
  void set_up_data_polygon(CShape* m_shape);
  //  void update_actors(float layerID,int lineID);

  std::vector< vtkSmartPointer<vtkActor> > vec_actor;
  vtkSmartPointer<vtkPoints>m_points;
  vtkSmartPointer<vtkCellArray>m_polygons;
  void set_up_data(CShape*);
private:
  int counter;
};
