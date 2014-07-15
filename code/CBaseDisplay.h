#pragma once
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
class CBaseDisplay
{
public:
  vtkSmartPointer<vtkActor> m_actor;
  vtkSmartPointer<vtkPolyData> m_polydata;
  vtkSmartPointer<vtkPolyDataMapper> m_polydata_mapper;
  vtkSmartPointer<vtkPoints>m_points;
  vtkSmartPointer<vtkCellArray>m_lines;
  vtkSmartPointer<vtkCellArray>m_vertices;
  CBaseDisplay()
  {
    m_lines=vtkSmartPointer<vtkCellArray>::New();
    m_vertices=vtkSmartPointer<vtkCellArray>::New();
    m_points=vtkSmartPointer<vtkPoints>::New();
    m_actor=vtkSmartPointer<vtkActor>::New();
    m_polydata=vtkSmartPointer<vtkPolyData>::New();
    m_polydata_mapper=vtkSmartPointer<vtkPolyDataMapper>::New();
    
    int n=999999;
    max_x=-n;min_x=n;
    max_y=-n;min_y=n;
    max_z=-n;min_z=n;  
  }
  double max_x,min_x,min_y,max_y,min_z,max_z;

  void check_edge(double x,double y,double z)
  {
    
    if(x>max_x)
    {
      max_x=x;
    }
    else if(x<min_x)
    {
      min_x=x;
    }

    if (y>max_y)
    {
      max_y=y;
    }
    else if (y<min_y)
    {
      min_y=y;
    }
    
    if (z>max_z)
    {
      max_z=z;
    }
    else if(z<min_z)
    {
      min_z=z;
    }
    
  }
};
