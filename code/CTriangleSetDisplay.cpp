#include "CTriangleSetDisplay.h"
#include "vtkExtractEdges.h"
#include "vtkDecimatePro.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkTriangleFilter.h"
CTriangleSetDisplay::CTriangleSetDisplay()
{
  m_mesh_actor=vtkSmartPointer<vtkActor>::New();
}

void CTriangleSetDisplay::initial_actors(CCorrespond* cor)
{
  std::map<int,Parapoint*>& tvec_corpoint = cor->map_cor;

  std::map<int,Parapoint*>::iterator itr,etr;
  itr= tvec_corpoint.begin();
  etr= tvec_corpoint.end();
  int count=0;
  for(;itr!=etr;++itr)
  {
        m_points->InsertPoint(static_cast<vtkIdType>(count++),itr->second->point1.x,itr->second->point1.y,itr->second->point1.z);
        m_points->InsertPoint(static_cast<vtkIdType>(count++),itr->second->point2.x,itr->second->point2.y,itr->second->point2.z);

  }
  // the triangle sequnce n+0,n+1,n+2--n+1,n+3,n+2: so n must smaller than count-3 to assure that n+3 is exist.
  for (int i = 0; i < count; ++i)
  {
    vtkSmartPointer<vtkTriangle> triangle=vtkSmartPointer<vtkTriangle>::New();
    triangle->GetPointIds()->SetId(0,i%count);
    triangle->GetPointIds()->SetId(1,(i+1)%count);
    triangle->GetPointIds()->SetId(2,(i+2)%count);

    m_lines->InsertNextCell(triangle);
    triangle=vtkSmartPointer<vtkTriangle>::New();
    triangle->GetPointIds()->SetId(0,(i+1)%count);
    triangle->GetPointIds()->SetId(1,(i+3)%count);
    triangle->GetPointIds()->SetId(2,(i+2)%count);

    m_lines->InsertNextCell(triangle);
  }
  

  
  //Create a mapper and Actors
  m_polydata->SetPoints(m_points);
  m_polydata->SetPolys(m_lines);

  vtkSmartPointer<vtkTriangleFilter> m_triangles=vtkSmartPointer<vtkTriangleFilter>::New();
  // m_triangles->DebugOn();
  m_triangles->SetInput(m_polydata);
  m_triangles->Update();
  //m_triangles->PassVertsOff();
  
  // vtkSmartPointer<vtkDecimatePro> m_deci= vtkSmartPointer<vtkDecimatePro>::New();
  // m_deci->DebugOn();
  // m_deci->SetInputConnection(m_triangles->GetOutputPort());
  // m_deci->SetTargetReduction(0.9);
  // m_deci->PreserveTopologyOn();
  // m_deci->Update();
  // m_deci->Print(std::cout);

  vtkSmartPointer<vtkSmoothPolyDataFilter> smooth= vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
  //  smooth->DebugOn();
  smooth->SetInputConnection(m_triangles->GetOutputPort());
  // smooth-> SetBoundarySmoothing(1);
  smooth->SetNumberOfIterations(50);
  smooth->Update();
  //  smooth->Print(std::cout);

  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
      vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  //  smoother->DebugOn();
  smoother->SetInputConnection(smooth->GetOutputPort());
  smoother->SetNumberOfIterations(100);
  smoother->BoundarySmoothingOff();
  smoother->FeatureEdgeSmoothingOff();
  smoother->SetFeatureAngle(45.0);
  smoother->SetPassBand(.001);
  //smoother->NonManifoldSmoothingOn();
  // smoother->NormalizeCoordinatesOn();
  smoother->Update();
  //  smoother->Print(std::cout);
  
  vtkSmartPointer<vtkPolyDataNormals> normals=vtkSmartPointer<vtkPolyDataNormals>::New();
  //  normals->DebugOn();
  normals->SetInputConnection(smoother->GetOutputPort());
  normals->FlipNormalsOn();
  normals->Update();
  //  normals->Print(std::cout);
  m_polydata_mapper->SetInputConnection(normals->GetOutputPort());
  m_actor->SetMapper(m_polydata_mapper);
  //Create mesh
  vtkSmartPointer<vtkExtractEdges> edges = vtkSmartPointer<vtkExtractEdges>::New();
  edges->SetInput(m_polydata);
  vtkSmartPointer<vtkPolyDataMapper> mesh_maper=vtkSmartPointer<vtkPolyDataMapper>::New();
  mesh_maper->SetInput(edges->GetOutput());
  m_mesh_actor->SetMapper(mesh_maper);
}
