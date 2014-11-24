#include "CContourDisplay.h"
#include "CPoint.h"
#include "CLayer.h"
#include "vtkVertex.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
extern double contour_line_color[3];
extern int linewidth;
vtkSmartPointer<vtkActor> CContourDisplay::initialActor(float layerID,int contourID)
{
  CLayer* temp = map_Layer[layerID];
  CContour* temp_contour= temp->map_contour[temp->map_contourID[contourID]];
  set_iterator(temp_contour,0);
  initial_contour_actor(point_itr,point_etr,point_size,1);
  current_LayerID=layerID;
  return m_actor;
}

void CContourDisplay::initial_contour_actor(std::vector<CPoint*>::iterator point_itr,std::vector<CPoint*>::iterator point_etr,int size,double level)
{
  vtkIdType * lineIndices = new vtkIdType[size+1];
  int count=0;
  
  for (;point_itr!=point_etr;++point_itr)
  {
    m_points->InsertPoint(static_cast<vtkIdType>(count),(*point_itr)->x,(*point_itr)->y,level*(*point_itr)->z);
    lineIndices[count] = static_cast<vtkIdType>(count);
      ++count;
  }
  m_vertices->InsertNextCell(size,lineIndices);
  
  lineIndices[count] = 0;
  m_lines->InsertNextCell(size+1,lineIndices);

  
  delete [] lineIndices;

  m_polydata->SetPoints(m_points);
  m_polydata->SetLines(m_lines);
  
  if( show_vertex_on_edge)
  {
    m_polydata->SetVerts(m_vertices);
    // double grey[3]={0.5,0.5,0.5};
    // vtkSmartPointer<vtkDoubleArray> colors= vtkSmartPointer<vtkDoubleArray>::New();
    // colors->SetName("colors");
    // colors->InsertNextTupleValue(grey);
    // m_polydata->GetPointData()->SetColor(colors);
  }

  m_polydata_mapper->SetInput(m_polydata);

  m_actor->SetMapper(m_polydata_mapper);
  m_actor->GetProperty()->SetLineWidth(linewidth);
  m_actor->GetProperty()->SetPointSize(5);
  m_actor->GetProperty()->SetColor(contour_line_color);
}

void CContourDisplay::update_contour_polydata(float layerID,int contourID)
{
  current_LayerID=layerID;
  CLayer* temp = map_Layer[layerID];
  CContour* temp_contour= temp->map_contour[temp->map_contourID[contourID]];
  set_iterator(temp_contour,0);
  int count=0;
  update_polydata(point_itr,point_etr,point_size,1);
}

void CContourDisplay::update_polydata(std::vector<CPoint*>::iterator point_itr,std::vector<CPoint*>::iterator point_etr,int size,double level)
{
  int count=0;
  vtkIdType * lineIndices = new vtkIdType[size+1];
  m_points->Reset();
  m_lines->Reset();
  m_vertices->Reset();

      
  for (;point_itr!=point_etr;++point_itr)
  {

    m_points->InsertPoint(static_cast<vtkIdType>(count),(*point_itr)->x,(*point_itr)->y,(*point_itr)->z*level);
    lineIndices[count] = static_cast<vtkIdType>(count);

    ++count;
  }
  m_vertices->InsertNextCell(size,lineIndices);
  lineIndices[count] = 0;
  m_lines->InsertNextCell(size+1,lineIndices);

  delete [] lineIndices;
  m_polydata->SetPoints(m_points);
  m_polydata->SetLines(m_lines);

  if( show_vertex_on_edge)
  {
    m_polydata->SetVerts(m_vertices);    
  }

  m_polydata->Modified();
  m_polydata->Update();  
}

void CContourDisplay::set_iterator(CContour* temp,int deformID)
{
  point_itr= temp->vec_points.begin();
  point_etr= temp->vec_points.end();
  point_size= temp->vec_points.size();
}
