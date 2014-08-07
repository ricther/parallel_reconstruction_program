#include "CTriangleDisplay.h"
#include "CShape.h"
#include "CTriangleSetDisplay.h"
#include "vtkDecimatePro.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkTriangleFilter.h"
#include "CTriangle.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkPolyDataWriter.h"
#include "vtkNew.h"
#include "vtkExtractEdges.h"
#include "CVirtualContour.h"
using namespace std;
CTriangleDisplay::CTriangleDisplay()
{
  float n=999999;
  max_z=-n;min_z=n;
  margin=0;
  map_layer_triangles_actors.clear();
  count=0;
  m_triangles=vtkSmartPointer<vtkCellArray>::New();
  m_points=vtkSmartPointer<vtkPoints>::New();
  m_actor=vtkSmartPointer<vtkActor>::New();
  m_mesh_actor=vtkSmartPointer<vtkActor>::New();
  m_polydata=vtkSmartPointer<vtkPolyData>::New();
  m_polydata_mapper=vtkSmartPointer<vtkPolyDataMapper>::New();

}
void CTriangleDisplay:: initial_actors(CShape * m_shape)
{

  int size= m_shape->vec_Cor.size();
  for (int i = 0; i < size; ++i)
  {
    CTriangleSetDisplay* new_triangle_set=new CTriangleSetDisplay();
    new_triangle_set->initial_actors(m_shape->vec_Cor[i]); 
    map_layer_triangles_actors.insert(make_pair(m_shape->vec_Cor[i]->CorrespondLayer1,new_triangle_set));
    check_level(m_shape->vec_Cor[i]->CorrespondLayer1);
    check_level(m_shape->vec_Cor[i]->CorrespondLayer2);
  }
}

void CTriangleDisplay:: initial_actor(CShape* m_shape)
{
  int size= m_shape->vec_Cor.size();
  set_up_data(m_shape);
  for (int i = 0; i < size; ++i)
  {
    organize_data(m_shape->vec_Cor[i]);
    check_level(m_shape->vec_Cor[i]->CorrespondLayer1);
    check_level(m_shape->vec_Cor[i]->CorrespondLayer2);
  }
  set_up_vtk();
}

void CTriangleDisplay:: set_up_data(CShape* m_shape)
{
  std::vector<CPoint*>::iterator p_itr,p_etr;
  std::map<float,CLayer*>::iterator l_itr,l_etr;
  std::map<int,CContour*>::iterator c_itr,c_etr;
  l_itr=m_shape->map_Layer.begin();l_etr=m_shape->map_Layer.end();
  int count=0;
  for (; l_itr!=l_etr;++l_itr)
  {
    c_itr=l_itr->second->map_contour.begin();c_etr=l_itr->second->map_contour.end();
    for (; c_itr!=c_etr; ++c_itr)
    {
      p_itr=c_itr->second->vec_Points_Origin.begin();p_etr=c_itr->second->vec_Points_Origin.end();
      for (; p_itr!=p_etr; ++p_itr)
      {
        m_points->InsertPoint(static_cast<vtkIdType>(count++),(*p_itr)->x,(*p_itr)->y,(*p_itr)->z);
        map_pindex_vtkindex.insert(make_pair((*p_itr)->index,count-1));
      }

      int N=c_itr->second->vec_virtual_contour.size();
      for (int i = 0; i < N; ++i)
      {
        p_itr=c_itr->second->vec_virtual_contour[0]->vec_Points_project.begin();p_etr=c_itr->second->vec_virtual_contour[0]->vec_Points_project.end();
        for (; p_itr!=p_etr; ++p_itr)
        {
          m_points->InsertPoint(static_cast<vtkIdType>(count++),(*p_itr)->x,(*p_itr)->y,(*p_itr)->z);
          map_pindex_vtkindex.insert(make_pair((*p_itr)->index,count-1));
        }
      }
    
    }
    if (l_itr->second->vec_medial_points.size()>0)
    {
      int N=l_itr->second->vec_medial_points.size();
      for (int i = 0; i < N; ++i)
      {
        m_points->InsertPoint(static_cast<vtkIdType>(count++),l_itr->second->vec_medial_points[i]->x,l_itr->second->vec_medial_points[i]->y,l_itr->second->vec_medial_points[i]->z);
        map_pindex_vtkindex.insert(make_pair(l_itr->second->vec_medial_points[i]->index,count-1));
      }
    }
  }
  vtkNew<vtkPolyDataWriter> writer;
  vtkSmartPointer<vtkPolyData> temp_polydata = vtkSmartPointer<vtkPolyData>::New();
  temp_polydata->SetPoints(m_points);
  writer->SetFileName("./points.vtk");
  writer->SetInput(temp_polydata);
  //  writer->Write();

}
void CTriangleDisplay:: organize_data(CCorrespond* cor)
{
  // std::map<int,Parapoint*>& tvec_corpoint = cor->map_cor;
  // std::map<int,Parapoint*>::iterator itr,etr;
  // itr= tvec_corpoint.begin();
  // etr= tvec_corpoint.end();
  // int sub_count=count;
  // for(;itr!=etr;++itr)
  // {
  //       m_points->InsertPoint(static_cast<vtkIdType>(count++),itr->second->point1.x,itr->second->point1.y,itr->second->point1.z);
  //       m_points->InsertPoint(static_cast<vtkIdType>(count++),itr->second->point2.x,itr->second->point2.y,itr->second->point2.z);
    
  // }
  // // the triangle sequnce n+0,n+1,n+2--n+1,n+3,n+2: so n must smaller than count-3 to assure that n+3 is exist.
  // for (int i = sub_count; i < count; ++i)
  // {
  //   vtkSmartPointer<vtkTriangle> triangle=vtkSmartPointer<vtkTriangle>::New();
  //   triangle->GetPointIds()->SetId(0,i%count);
  //   triangle->GetPointIds()->SetId(1,get_id(i+1,count,sub_count));
  //   triangle->GetPointIds()->SetId(2,get_id(i+2,count,sub_count));

  //   m_triangles->InsertNextCell(triangle);
  //   triangle=vtkSmartPointer<vtkTriangle>::New();
  //   triangle->GetPointIds()->SetId(0,get_id(i+2,count,sub_count));
  //   triangle->GetPointIds()->SetId(1,get_id(i+1,count,sub_count));
  //   triangle->GetPointIds()->SetId(2,get_id(i+3,count,sub_count));

  //   m_triangles->InsertNextCell(triangle);
  // }
  cor->produce_triangle();
  std::vector<CTriangle*>::iterator itr=cor->vec_Triangle.begin(),etr=cor->vec_Triangle.end();
  for (; itr!=etr; ++itr)
  {
    vtkSmartPointer<vtkTriangle> triangle=vtkSmartPointer<vtkTriangle>::New();
    int a=map_pindex_vtkindex[(*itr)->a],b=map_pindex_vtkindex[(*itr)->b],c=map_pindex_vtkindex[(*itr)->c];
    triangle->GetPointIds()->SetId(0,a);
    triangle->GetPointIds()->SetId(1,b);
    triangle->GetPointIds()->SetId(2,c);
    m_triangles->InsertNextCell(triangle);
  }
}

void CTriangleDisplay:: set_up_vtk()
{
  m_polydata->SetPoints(m_points);
  m_polydata->SetPolys(m_triangles);
  smooth();
  //  m_polydata_mapper->SetInput(m_polydata);
  m_actor->SetMapper(m_polydata_mapper);
  m_actor->GetProperty()->SetColor(1.0,0.4,0);
  m_actor->GetProperty()->SetAmbient(0.5);
  m_actor->GetProperty()->SetDiffuse(0.5);
  m_actor->GetProperty()->SetSpecular(0.5);
  //  m_actor->GetProperty()->SetRepresentationToWireframe();

  //Create mesh
  vtkSmartPointer<vtkExtractEdges> edges = vtkSmartPointer<vtkExtractEdges>::New();
  edges->SetInput(m_polydata);
  vtkSmartPointer<vtkPolyDataMapper> mesh_maper=vtkSmartPointer<vtkPolyDataMapper>::New();
  mesh_maper->SetInput(edges->GetOutput());
  m_mesh_actor->SetMapper(mesh_maper);

}

void CTriangleDisplay:: smooth()
{
  vtkSmartPointer<vtkDataSetSurfaceFilter> surface= vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surface->SetInput(m_polydata);
  surface->Update();
  vtkNew<vtkPolyDataWriter> writer;
  writer->SetFileName("./surface.vtk");
  writer->SetInputConnection(surface->GetOutputPort());
  //  writer->Write();
  //  vtkSmartPointer<vtkTriangleFilter> m_triangles=vtkSmartPointer<vtkTriangleFilter>::New();
  // m_triangles->DebugOn();
  //  m_triangles->SetInput(m_polydata);
  //  m_triangles->Update();
  //m_triangles->PassVertsOff();
  
  // vtkSmartPointer<vtkDecimatePro> m_deci= vtkSmartPointer<vtkDecimatePro>::New();
  // m_deci->DebugOn();
  // m_deci->SetInputConnection(surface->GetOutputPort());
  // m_deci->SetTargetReduction(0.1);
  // m_deci->PreserveTopologyOn();
  // m_deci->Update();
  // m_deci->Print(std::cout);

  // vtkSmartPointer<vtkSmoothPolyDataFilter> smooth= vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
  // //  smooth->DebugOn();
  // smooth->SetInputConnection(m_deci ->GetOutputPort());
  // // smooth-> SetBoundarySmoothing(1);
  // smooth->SetNumberOfIterations(200);
  // smooth->Update();
  // //  smooth->Print(std::cout);

  // vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
  //     vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  // //  smoother->DebugOn();
  // smoother->SetInputConnection(smooth->GetOutputPort());
  // smoother->SetNumberOfIterations(100);
  // //  smoother->BoundarySmoothingOff();
  // smoother->FeatureEdgeSmoothingOff();
  // //  smoother->SetFeatureAngle(45.0);
  // //  smoother->SetPassBand(.001);
  // //smoother->NonManifoldSmoothingOn();
  // // smoother->NormalizeCoordinatesOn();
  // smoother->Update();

  // vtkSmartPointer<vtkPolyDataNormals> normals=vtkSmartPointer<vtkPolyDataNormals>::New();
  // //  normals->DebugOn();
  // normals->SetInputConnection(smoother->GetOutputPort());
  // normals->FlipNormalsOn();
  // normals->Update();

  m_polydata_mapper->SetInputConnection(surface->GetOutputPort());
}

int CTriangleDisplay:: get_id(int id,int count,int sub_count)
{
  if (id>=count)
  {
    id=id%count+sub_count;
  }
  return id;
}

void CTriangleDisplay::update_actors(float layerID,int lineID)
{
  std::multimap<float,CTriangleSetDisplay*>::iterator itr,etr;
  itr=map_layer_triangles_actors.begin();etr=map_layer_triangles_actors.end();
  for (;itr!=etr;++itr)
  {
    bool b=check_display(layerID,itr->first);
    if (b==true)
    {
      itr->second->set_actor_visibility(true,lineID);
    }
    else
    {
      itr->second->set_actor_visibility(false,lineID);
    }
  }
}
