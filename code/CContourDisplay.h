#pragma once
#include "CBaseDisplay.h"
#include <map>
#include "CContour.h"
extern const bool show_vertex_on_edge;
class CLayer;
class CPoint;
class CContourDisplay:public CBaseDisplay
{
public:
CContourDisplay(std::map<float,CLayer*>& temp_map):map_Layer(temp_map){};
  vtkSmartPointer<vtkActor> initialActor(float layerID,int contourID);
  void initial_contour_actor(std::vector<CPoint*>::iterator,std::vector<CPoint*>::iterator,int,double);
  void update_contour_polydata(float LayerID,int id);
  void update_polydata(std::vector<CPoint*>::iterator,std::vector<CPoint*>::iterator,int,double);
  virtual void set_iterator(CContour* temp);
  float current_LayerID;
protected:
  std::vector<CPoint*>::iterator point_itr,point_etr;
  std::map<float,CLayer*>& map_Layer;
  int point_size;
};
