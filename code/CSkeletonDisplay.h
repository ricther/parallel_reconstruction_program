#pragma once

#include "CBaseDisplay.h"
#include <map>
#include "vtkUnsignedCharArray.h"
#include "vtkSphereSource.h"
#include "CSkeleton.h"

class CLayer;
class CSkeletonDisplay:public CBaseDisplay
{
public:
  vtkSmartPointer<vtkActor> initialActor(CSkeleton* m_skeleton,vtkSmartPointer<vtkRenderer>& renderer);
  void update_actor(float current_LayerID,std::map<float,CLayer*>& map_Layer);
  vtkSmartPointer<vtkActor> m_skeletal_point_actor;
private:
  vtkSmartPointer<vtkActor> draw_branch(CBranch* m_branch);
  std::map<float ,int> map_LayerID_vtkID;
  unsigned char green[3],red[3];
  vtkSmartPointer<vtkSphereSource> skeletal_point;
};










