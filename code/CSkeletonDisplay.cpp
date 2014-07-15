#include "CSkeletonDisplay.h"
#include "CLayer.h"
#include "CPoint.h"
#include "vtkPointData.h"
#include "vtkRenderer.h"

using namespace std;
vtkSmartPointer<vtkActor> CSkeletonDisplay::initialActor(CSkeleton* m_skeleton,vtkSmartPointer<vtkRenderer>& renderer)
{

  green[0]=120,green[1]=176,green[2]=91;
  red[0]=255,green[1]=0;red[2]=0;
  
  std::map<int,CBranch*>::iterator itr=m_skeleton->map_id_branch.begin(), etr=m_skeleton->map_id_branch.end();
  int size=m_skeleton->map_id_branch.size();
  for (; itr!=etr; ++itr)
  {
    vtkSmartPointer<vtkActor> temp_actor = draw_branch(itr->second);
    renderer->AddActor(temp_actor);
  }

  //  skeletal_point=vtkSmartPointer<vtkSphereSource>::New();
  
  // temp=map_Layer.begin()->second;
  // skeletal_point->SetCenter(temp->moment_one_point->x,temp->moment_one_point->y, temp->LayerID);
  // skeletal_point->SetRadius(3.0);
 
  // vtkSmartPointer<vtkPolyDataMapper> mapper = 
  //   vtkSmartPointer<vtkPolyDataMapper>::New();
  // mapper->SetInputConnection(skeletal_point->GetOutputPort());
 
  // m_skeletal_point_actor=vtkSmartPointer<vtkActor>::New();
  // m_skeletal_point_actor->SetMapper(mapper);
  // m_skeletal_point_actor->GetProperty()->SetColor(red[0],red[1],red[2]);
}

void CSkeletonDisplay::update_actor(float current_LayerID,std::map<float,CLayer*>&  map_Layer)
{
  CLayer* temp=map_Layer[current_LayerID];
  skeletal_point->SetCenter(temp->moment_one_point->x,temp->moment_one_point->y,temp->LayerID);
  skeletal_point->Modified();
  skeletal_point->Update();
}

vtkSmartPointer<vtkActor> CSkeletonDisplay::draw_branch(CBranch* m_branch)
{
  std::map<int,CSkeletalPoint*>::iterator itr=m_branch->map_id_spoint.begin(),etr=m_branch->map_id_spoint.end();
  int size=m_branch->map_id_spoint.size();
  vtkIdType * lineIndices = new vtkIdType[size];
  int count=0;

  vtkSmartPointer<vtkActor> m_actor;
  vtkSmartPointer<vtkPolyData> m_polydata;
  vtkSmartPointer<vtkPolyDataMapper> m_polydata_mapper;
  vtkSmartPointer<vtkPoints>m_points;
  vtkSmartPointer<vtkCellArray>m_lines;

  m_lines=vtkSmartPointer<vtkCellArray>::New();
  m_points=vtkSmartPointer<vtkPoints>::New();
  m_actor=vtkSmartPointer<vtkActor>::New();
  m_polydata=vtkSmartPointer<vtkPolyData>::New();
  m_polydata_mapper=vtkSmartPointer<vtkPolyDataMapper>::New();

  for (; itr!=etr; ++itr)
  {
    CSkeletalPoint * temp= itr->second;

    m_points->InsertPoint(static_cast<vtkIdType>(count),temp->node->x,temp->node->y,temp->node->z);
    lineIndices[count] = static_cast<vtkIdType>(count);
    map_LayerID_vtkID.insert(std::make_pair(temp->node->z,count));
    //    colors->InsertTupleValue(static_cast<vtkIdType>(count),green);
    ++count;
    check_edge(temp->node->x,temp->node->y,temp->node->z);
  }

  m_lines->InsertNextCell(size,lineIndices);
  delete [] lineIndices;

  m_polydata->SetPoints(m_points);

  m_polydata->SetLines(m_lines);
  //  m_polydata->GetPointData()->SetScalars(colors);
  m_polydata_mapper->SetInput(m_polydata);
  
  m_actor->SetMapper(m_polydata_mapper);
  m_actor->GetProperty()->SetPointSize(3);
  
  return m_actor;
}
