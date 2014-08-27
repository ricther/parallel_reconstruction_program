#include "CContourSurface.h"
#include "COriginContour.h"
#include "vtkTriangleFilter.h"
#include "vtkDelaunay2D.h"
#include "CShape.h"
#include "vtkPolygon.h"

extern double contour_surface_color[3];
CContourSurface::CContourSurface()
{
    counter=0;
    m_polygons=  vtkSmartPointer<vtkCellArray>::New();
    m_points =  vtkSmartPointer<vtkPoints>::New();

}

void CContourSurface:: initial_actors(CShape * m_shape)
{
  //  set_up_data_Delaunay2D(m_shape);
  set_up_data_polygon(m_shape);

}


void CContourSurface:: set_up_data_polygon(CShape* m_shape)
{
  std::map<float,CLayer*>::iterator itr=m_shape->map_Layer.begin(),etr=m_shape->map_Layer.end();
  for (;itr!=etr;++itr)
  {
    std::map<int,int>::iterator itr_c,etr_c;
    itr_c=itr->second->map_contourID.begin();etr_c=itr->second->map_contourID.end();
    for (; itr_c!=etr_c; ++itr_c)
    {

      CContour* temp_contour= itr->second->map_contour[itr_c->second];
      int N = temp_contour->vec_Points_Origin.size();
      for (int i = 0; i < N; ++i)
      {
        m_points->InsertPoint(static_cast<vtkIdType>(counter++),temp_contour->vec_Points_Origin[i]->x,temp_contour->vec_Points_Origin[i]->y,temp_contour->vec_Points_Origin[i]->z);
      }
      vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
      polygon->GetPointIds()->SetNumberOfIds(N);
      for (int ii = 0; ii < N; ++ii)
      {
        polygon->GetPointIds()->SetId(ii,counter-N+ii);
      }
      m_polygons->InsertNextCell(polygon);
      //vtkSmartPointer<vtkIdList> id_list= vtkSmartPointer<vtkIdList>::New();
      //      polygon->Triangulate(id_list);

      //      vtkSmartPointer<vtkCellArray>m_lines=vtkSmartPointer<vtkCellArray>::New();
      // int size=id_list->GetNumberOfIds();
      // for (int n = 0; n < size; n+3)
      // {
      //   if (n+2>=size)
      //   {
      //     break;
      //   }
      //   vtkSmartPointer<vtkTriangle> triangle=vtkSmartPointer<vtkTriangle>::New();
      //   vtkIdType id=id_list->GetId(n);
      //   triangle->GetPointIds()->SetId(0,id);
      //   id=id_list->GetId(n+1);
      //   triangle->GetPointIds()->SetId(1,id);
      //   id=id_list->GetId(n+2);
      //   triangle->GetPointIds()->SetId(2,id);
      //   m_lines->InsertNextCell(triangle);
      // }
      
      
    }
  }
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(m_points);
  polydata->SetPolys(m_polygons);

  vtkSmartPointer<vtkTriangleFilter> triangleFilter= vtkSmartPointer<vtkTriangleFilter>::New();
  triangleFilter->SetInput(polydata);
  triangleFilter->Update();

  vtkSmartPointer<vtkPolyDataMapper> triangleMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  triangleMapper->SetInputConnection(triangleFilter->GetOutputPort());
      
  vtkSmartPointer<vtkActor> triangleActor= vtkSmartPointer<vtkActor>::New();
  triangleActor->SetMapper(triangleMapper);
  triangleActor->GetProperty()->SetColor(contour_surface_color);
  //  triangleActor->GetProperty()->SetOpacity(0.5);
  vec_actor.push_back(triangleActor);

}
void CContourSurface:: set_up_data_Delaunay2D(CShape *m_shape)
{

  std::map<float,CLayer*>::iterator itr=m_shape->map_Layer.begin(),etr=m_shape->map_Layer.end();
  for (;itr!=etr;++itr)
  {
    std::map<int,int>::iterator itr_c,etr_c;
    itr_c=itr->second->map_contourID.begin();etr_c=itr->second->map_contourID.end();
    for (; itr_c!=etr_c; ++itr_c)
    {
        COriginContour* new_contour=new COriginContour( m_shape->map_Layer);
        new_contour->initialActor(itr->first,itr_c->first);
        // vtkSmartPointer<vtkTriangleFilter> triangleFilter= vtkSmartPointer<vtkTriangleFilter>::New();
        // triangleFilter->SetInput(new_contour->m_polydata);
        // triangleFilter->Update();


        vtkSmartPointer<vtkDelaunay2D> Del = vtkSmartPointer<vtkDelaunay2D>::New();
        Del->SetInput(new_contour->m_polydata);
        Del->SetTolerance(0.01);
        
        vtkSmartPointer<vtkPolyDataMapper> triangleMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        triangleMapper->SetInputConnection(Del->GetOutputPort());

        vtkSmartPointer<vtkActor> triangleActor= vtkSmartPointer<vtkActor>::New();
        triangleActor->SetMapper(triangleMapper);
        vec_actor.push_back(triangleActor);
    }
  }
}

