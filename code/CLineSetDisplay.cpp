#include "CLineSetDisplay.h"
#include "assert.h"
void CLineSetDisplay::initial_actors(CCorrespond* cor)
{
  std::map<int,Parapoint*>& tvec_corpoint = cor->map_cor;

  std::map<int,Parapoint*>::iterator itr,etr;
  itr= tvec_corpoint.begin();
  etr= tvec_corpoint.end();
  for(;itr!=etr;++itr)
  {
    Parapoint* n=itr->second;
    vtkSmartPointer<vtkLineSource> LineSource=vtkSmartPointer<vtkLineSource>::New();
    LineSource->SetPoint1(n->point1.x,n->point1.y,n->point1.z);
    LineSource->SetPoint2(n->point2.x,n->point2.y,n->point2.z);

    vec_line_source.push_back(LineSource);
    
    //Create a mapper and Actors
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(LineSource->GetOutputPort());

    vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    
    actor->SetMapper(mapper);

    vec_actor.push_back(actor);
  }
}

CLineSetDisplay::CLineSetDisplay()
{
  m_lines=vtkSmartPointer<vtkCellArray>::New();
  m_points=vtkSmartPointer<vtkPoints>::New();
  m_actor=vtkSmartPointer<vtkActor>::New();
  m_polydata=vtkSmartPointer<vtkPolyData>::New();
  m_polydata_mapper=vtkSmartPointer<vtkPolyDataMapper>::New();
}
void CLineSetDisplay::initial_actor(CCorrespond* cor)
{
  std::map<int,Parapoint*>& tvec_corpoint = cor->map_cor;

  std::map<int,Parapoint*>::iterator itr,etr;
  itr= tvec_corpoint.begin();
  etr= tvec_corpoint.end();
  int count=0;
  int counter=0;
  for(;itr!=etr;++itr)
  {
    counter++;
    if (counter%interval_of_line_layer!=0)
    {
      continue;
    }
    Parapoint* n=itr->second;
    m_points->InsertPoint(static_cast<vtkIdType>(count++),n->point1.x,n->point1.y,n->point1.z);
    m_points->InsertPoint(static_cast<vtkIdType>(count++),n->point2.x,n->point2.y,n->point2.z);
    m_lines->InsertNextCell(2);
    m_lines->InsertCellPoint(count-2);
    m_lines->InsertCellPoint(count-1);
  }
  m_polydata->SetPoints(m_points);
  m_polydata->SetLines(m_lines);
  m_polydata_mapper->SetInput(m_polydata);
  m_actor->SetMapper(m_polydata_mapper);
}
