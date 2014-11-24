#include"CShapeDisplay.h"
#include "CPoint.h"
#include <vtkProperty.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkTextProperty.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include "KeyPressInteractorStyle.h"
#include <vtkRenderWindow.h>
#include "CShape.h"
#include "CLayer.h"
#include "vtkProperty2D.h"
#include "CLineSetDisplay.h"

using namespace std;
CShapeDisplay::CShapeDisplay(CShape* temp,vtkSmartPointer<vtkRenderWindow> temp_window,vtkSmartPointer<vtkRenderWindowInteractor> temp_interactor)
{
  m_shape=temp;
  min_z=m_shape->min_z;max_z=m_shape->max_z;
  m_skeleton = new CSkeletonDisplay();
  m_skeleton_contour= new CContourDisplay(m_shape->map_Layer);
  m_correspond = new CCorespondDisplay();
  m_origin_contours=new CContourSetDisplay<COriginContour>(m_shape->map_Layer);
  m_normal_contours=new CContourSetDisplay<CNormalContour>(m_shape->map_Layer);
  m_transformed_contours=new CContourSetDisplay<CTransformedContour>(m_shape->map_Layer);
  m_triangle = new CTriangleDisplay();
  m_contour_surface= new CContourSurface();
  m_layer_text_actor= vtkSmartPointer<vtkTextActor>::New();
  
  m_render_window=temp_window;
  m_interactor=temp_interactor;
  b_draw_contour_for_skeleton=false;
  count_contour=0;
  current_LayerID=1;
  current_lineID=0;
  current_correspondID=0;
}

void CShapeDisplay::draw_correspond( vtkSmartPointer<vtkRenderer> renderer)
{
  // m_correspond->initial_actors(m_shape);
  // std::map<float,CLineSetDisplay*>::iterator itr,etr;
  // itr = m_correspond->map_line_actors.begin();etr = m_correspond->map_line_actors.end();
  // for (;itr!=etr; ++itr)
  // {
  //   std::vector< vtkSmartPointer<vtkActor> >::iterator itr_actor,etr_actor;
  //   itr_actor=(itr->second)->vec_actor.begin();
  //   etr_actor=(itr->second)->vec_actor.end();
  //   for (;itr_actor!=etr_actor; ++itr_actor)
  //   {
  //     renderer->AddActor(*itr_actor);
  //   }
  // }
  m_correspond->initial_actor(m_shape);
  m_correspond->draw_by_layer=true;
  std::map<float,CLineSetDisplay*>::iterator itr,etr;
  itr = m_correspond->map_line_actors.begin();etr = m_correspond->map_line_actors.end();
  for (;itr!=etr; ++itr)
  {
    renderer->AddActor(itr->second->m_actor);
  }
  //  draw_medial_axis(renderer);
  draw_text(renderer,string("correspond links"));
  draw_origin_points(renderer);
}

void CShapeDisplay::draw_text(vtkSmartPointer<vtkRenderer> renderer, string str)
{
    // Setup the text and add it to the window
  vtkSmartPointer<vtkTextActor> textActor = 
    vtkSmartPointer<vtkTextActor>::New();
  textActor->GetTextProperty()->SetFontSize ( 16 );
  textActor->SetPosition2 ( 0, 0 );
  renderer->AddActor2D ( textActor );
  textActor->SetInput (str.c_str());
  textActor->GetTextProperty()->SetColor ( 1.0,0.0,0.0 );
  line_space=textActor->GetTextProperty()->GetLineSpacing();
 }

#include <sstream>
void CShapeDisplay::draw_layerID(vtkSmartPointer<vtkRenderer> renderer ,float layerID)
{

  m_layer_text_actor->GetTextProperty()->SetFontSize (15);
  m_layer_text_actor->SetPosition(0,8*line_space);
  renderer->AddActor2D ( m_layer_text_actor );
  
  string str;
  str="current_layerID:";
  std::stringstream stream;
  stream<<layerID;

  str+=stream.str();
  m_layer_text_actor->SetInput (str.c_str());
  m_layer_text_actor->GetTextProperty()->SetColor ( 1.0,0.0,0.0 );
  m_layer_text_actor->Modified();
}

void CShapeDisplay::update_layer_text(float layerID)
{
  string str;
  str="current_layerID:";

  std::stringstream stream;
  stream<<layerID;
  
  str+=stream.str();
  m_layer_text_actor->SetInput (str.c_str());
  m_layer_text_actor->Modified();
}

void CShapeDisplay::draw_origin_points(vtkSmartPointer<vtkRenderer> renderer)
{
  m_origin_contours->initialActors();
  std::vector<COriginContour*>::iterator itr, etr;
  itr=m_origin_contours->vec_contour.begin();
  etr=m_origin_contours->vec_contour.end();
  for (; itr!=etr; ++itr)
  {
    renderer->AddActor((*itr)->m_actor);
  }
  vtkSmartPointer<KeyPressInteractorStyle> style= vtkSmartPointer<KeyPressInteractorStyle>::New(); 
  m_interactor->SetInteractorStyle(style);
  style->SetShape(this);
  //  draw_lattice(renderer);
  draw_text(renderer,"original points");
}

void CShapeDisplay::draw_normal_points(vtkSmartPointer<vtkRenderer> renderer)
{
  m_normal_contours->initialActors();
  std::vector<CNormalContour*>::iterator itr,etr;
  itr=m_normal_contours->vec_contour.begin();etr =m_normal_contours->vec_contour.end();
  for (;itr!=etr;++itr)
  {
    //    (*itr)->m_actor->GetProperty()->SetColor(1,1,0);
    renderer->AddActor((*itr)->m_actor);
  }
  //   draw_medial_axis(renderer);
  //  draw_lattice(renderer);
  vtkSmartPointer<KeyPressInteractorStyle> style= vtkSmartPointer<KeyPressInteractorStyle>::New(); 
  m_interactor->SetInteractorStyle(style);
  style->SetCurrentRenderer(renderer);
  style->SetShape(this);

  draw_text(renderer,"normal points");
  draw_layerID(renderer,1.0);
}

void CShapeDisplay::draw_lattice(vtkSmartPointer<vtkRenderer> renderer)
{
  CLayer* temp = m_shape->map_Layer[20];
  CContour* temp_contour=temp->map_contour[temp->map_contourID[1]];
  draw_lattice_content(renderer,temp_contour,1);
  draw_lattice_content(renderer,temp_contour,2);
}

void CShapeDisplay::draw_lattice_content(vtkSmartPointer<vtkRenderer> renderer,CContour* contour,int mode)//mode 1 draw x axis mode 2 draw y axis
{
  int N=MatrixRes;
  if(mode==1)
  {
    for (int i = 0; i < N; ++i)
    {
      vtkSmartPointer<vtkActor> actor= vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkPolyData> m_polydata= vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPolyDataMapper> mapper= vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkPoints> m_points= vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkCellArray> m_lines= vtkSmartPointer<vtkCellArray>::New();
      vtkIdType* lineIndices=new vtkIdType[MatrixRes];
      
      int count=0;
      for (int j = 0; j < N; ++j)
      {
       m_points->InsertPoint(static_cast<vtkIdType>(count),contour->vec_virtual_contour[0]->lattice_x[i][j],contour->vec_virtual_contour[0]->lattice_y[i][j],contour->LayerID);
        lineIndices[count] = static_cast<vtkIdType>(count);
        count++;
      }
      m_lines->InsertNextCell(MatrixRes,lineIndices);
      delete[] lineIndices;

      m_polydata->SetPoints(m_points);
      m_polydata->SetLines(m_lines);
      mapper->SetInput(m_polydata);
      actor->SetMapper(mapper);
      actor->GetProperty()->SetColor(0,0,0);
      renderer->AddActor(actor);
    }
  }
  else
  {
    for (int j = 0; j < N; ++j)
    {
      vtkSmartPointer<vtkActor> actor= vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkPolyData> m_polydata= vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPolyDataMapper> mapper= vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkPoints> m_points= vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkCellArray> m_lines= vtkSmartPointer<vtkCellArray>::New();
      vtkIdType* lineIndices=new vtkIdType[MatrixRes];
      
      int count=0;
      for (int i = 0; i < N; ++i)
      {
        m_points->InsertPoint(static_cast<vtkIdType>(count),contour->vec_virtual_contour[0]->lattice_x[i][j],contour->vec_virtual_contour[0]->lattice_y[i][j],contour->LayerID);
        lineIndices[count] = static_cast<vtkIdType>(count);
        count++;
      }
      m_lines->InsertNextCell(MatrixRes,lineIndices);
      delete[] lineIndices;

      m_polydata->SetPoints(m_points);
      m_polydata->SetLines(m_lines);
      mapper->SetInput(m_polydata);
      actor->SetMapper(mapper);
      actor->GetProperty()->SetColor(0,0,0);
      renderer->AddActor(actor);
    }
  }
}

void CShapeDisplay::draw_medial_axis(vtkSmartPointer<vtkRenderer> renderer)
{
  CVirtualContour* temp = m_shape->contour_arrangement->vec_contour_couple[1]->vcontour2;
  vector<CPoint*> vec=temp->vec_medial_points;
  int size= vec.size();
  for (int i = 0; i < size; ++i)
  {
    vtkSmartPointer<vtkSphereSource> point=vtkSmartPointer<vtkSphereSource>::New();
    point->SetCenter(vec[i]->x,vec[i]->y,vec[i]->z);
    point->SetRadius(0.5);
 
    vtkSmartPointer<vtkPolyDataMapper> mapper = 
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(point->GetOutputPort());
 
    vtkSmartPointer<vtkActor>point_actor=vtkSmartPointer<vtkActor>::New();
    point_actor->SetMapper(mapper);
    point_actor->GetProperty()->SetColor(255,0,0);
    renderer->AddActor(point_actor);
  }
}

void CShapeDisplay::draw_transformed_points(vtkSmartPointer<vtkRenderer> renderer)
{
  m_transformed_contours->initialActors();
  std::vector<CTransformedContour*>::iterator itr, etr;
  itr=m_transformed_contours->vec_contour.begin();
  etr=m_transformed_contours->vec_contour.end();
  for (; itr!=etr; ++itr)
  {
    renderer->AddActor((*itr)->m_actor);
  }
  //  draw_medial_axis(renderer);
  //  draw_lattice(renderer);
  draw_text(renderer,"transformed_points|current_layer_in_red");
  draw_layerID(renderer,1.0);
}

void CShapeDisplay::draw_compared_points(vtkSmartPointer<vtkRenderer> renderer)
{
  CLayer* temp;

  // itr= m_shape->map_Layer.begin(); etr = m_shape->map_Layer.end();

  // for (;itr!=etr;++itr)
  // {
  //   temp = itr->second;
  //   int size = temp->vec_new_points.size();
  //   if (size==0)
  //   {
  //     continue;
  //   }
  //   init_contour_represetation(renderer,temp->vec_Points_Origin,size,temp->LayerID*10);
  //   init_contour_represetation(renderer,temp->vec_new_points,size,temp->LayerID*15);
  // }
  draw_text(renderer,"same contour compare"); 
}

void CShapeDisplay::init_contour_represetation(vtkSmartPointer<vtkRenderer> renderer,vector<CPoint*> vec_p,int scale,double level)
{
  std::vector<CPoint*>::iterator point_itr,point_etr;
  
    
  point_itr=vec_p.begin();point_etr=vec_p.end();
  int size=vec_p.size();
  vtkSmartPointer<vtkActor> actor=initialActor(point_itr,point_etr,size,level);
  renderer->AddActor(actor);
}

vtkSmartPointer<vtkActor> CShapeDisplay::initialActor(std::vector<CPoint*>::iterator point_itr,std::vector<CPoint*>::iterator point_etr,int size,double level)
{
  vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkIdType * lineIndices = new vtkIdType[size+1];
  int count=0;
  for (;point_itr!=point_etr;++point_itr)
  {
    if(level==1)
      level=(*point_itr)->z;
    points->InsertPoint(static_cast<vtkIdType>(count),(*point_itr)->x,(*point_itr)->y,level );
    lineIndices[count] = static_cast<vtkIdType>(count);
    ++count;
  }
  lineIndices[count] = 0;
  lines->InsertNextCell(size+1,lineIndices);
  delete [] lineIndices;
  pd->SetPoints(points);
  pd->SetLines(lines);

  vtkSmartPointer<vtkPolyDataMapper> mapper=vtkSmartPointer<vtkPolyDataMapper>::New();

  mapper->SetInput(pd);
  vtkSmartPointer<vtkActor>actor= vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  return actor;
}

void CShapeDisplay::draw_skeleton(vtkSmartPointer<vtkRenderer> renderer,vtkSmartPointer<vtkRenderWindowInteractor> interactor)
{
  b_draw_contour_for_skeleton=false;//turn the switch for the draw_contour_for_skeleton

  m_skeleton->initialActor(m_shape->m_skeleton,renderer);

  if (b_draw_contour_for_skeleton)
  {
    m_skeleton_contour->initialActor(1,0);    
  }


  renderer->AddActor(m_skeleton_contour->m_actor);
  
  //  renderer->AddActor(m_skeleton->m_actor);

  //  renderer->AddActor(m_skeleton->m_skeletal_point_actor);
  draw_text(renderer,"skeleton");
  draw_layerID(renderer,1.0);
}
void CShapeDisplay::draw_triangle(vtkSmartPointer<vtkRenderer> renderer,vtkSmartPointer<vtkRenderWindowInteractor> interactor)
{
  // m_triangle->initial_actors(m_shape);
  // std::map<float,CTriangleSetDisplay*>::iterator itr,etr;
  // itr = m_triangle->map_layer_triangles_actors.begin();etr = m_triangle->map_layer_triangles_actors.end();
  // for (;itr!=etr; ++itr)
  // {//m_actor draw the surface
  //   renderer->AddActor(itr->second->m_actor);
  //   //m_mesh_actor draw the mesh
  //   //    renderer->AddActor(itr->second->m_mesh_actor);
  // }

  
  m_triangle->initial_actor(m_shape);
  renderer->AddActor(m_triangle->m_actor);
  //  renderer->AddActor(m_triangle->m_mesh_actor);

  if(show_contour_surface)
  {
    m_contour_surface->initial_actors(m_shape);
    int N= m_contour_surface->vec_actor.size();
    for (int i = 0; i < N; ++i)
    {
      renderer->AddActor(m_contour_surface->vec_actor[i]);
    }
  }
  draw_text(renderer,string("triangle surface"));

}
void CShapeDisplay:: key_down()
{
  current_LayerID=m_shape->get_next_layer(current_LayerID,0);
  if (current_LayerID<min_z)
  {
    current_LayerID=min_z;
  }
  //  m_skeleton->update_actor(current_LayerID,m_shape->map_Layer);
  if (b_draw_contour_for_skeleton)
  {
    m_skeleton_contour->update_contour_polydata(current_LayerID,0);    
  }

  m_origin_contours->update_actor(current_LayerID);
  m_normal_contours->update_actor(current_LayerID);
  m_triangle->update_actors(current_LayerID,0);
  m_correspond->update_actors(current_LayerID,current_lineID);
  m_transformed_contours->update_actor(current_LayerID);
  update_layer_text(current_LayerID);
  m_render_window->Render();
}

void CShapeDisplay:: key_up()
{
  current_LayerID=m_shape->get_next_layer(current_LayerID,1);
  if (current_LayerID>max_z)
  {
    current_LayerID=max_z;
  }
  //  m_skeleton->update_actor(current_LayerID,m_shape->map_Layer);
  if (b_draw_contour_for_skeleton)
  {
    m_skeleton_contour->update_contour_polydata(current_LayerID,0);    
  }

  m_origin_contours->update_actor(current_LayerID);
  m_normal_contours->update_actor(current_LayerID);
  m_triangle->update_actors(current_LayerID,0);
  m_correspond->update_actors(current_LayerID,current_lineID);
  m_transformed_contours->update_actor(current_LayerID);
  update_layer_text(current_LayerID);
  m_render_window->Render();
}

void CShapeDisplay::key_minus()
{
  float temp = m_origin_contours->margin;
  temp=temp*2;
  if(temp<0)
  {
    temp=min_z;
  }
  m_origin_contours->margin=temp;
  m_normal_contours->margin=temp;
  m_transformed_contours->margin=temp;
  m_origin_contours->update_actor(current_LayerID);
  m_normal_contours->update_actor(current_LayerID);
  m_transformed_contours->update_actor(current_LayerID);
  m_render_window->Render();
}
void CShapeDisplay::key_equal()
{
  float temp = m_origin_contours->margin;
  temp=temp*2;
  if(temp>max_z)
  {
    temp=max_z;
  }
  m_origin_contours->margin=temp;
  m_normal_contours->margin=temp;
  m_transformed_contours->margin=temp;
  m_origin_contours->update_actor(current_LayerID);
  m_normal_contours->update_actor(current_LayerID);
  m_transformed_contours->update_actor(current_LayerID);
  m_render_window->Render();
}
void CShapeDisplay::key_r()
{
  m_origin_contours->update_actor(current_LayerID);
  m_normal_contours->update_actor(current_LayerID);
  m_transformed_contours->update_actor(current_LayerID);
  m_render_window->Render();
}
void CShapeDisplay::key_f()
{
  m_origin_contours->show_all_actor();
  m_normal_contours->show_all_actor();
  m_transformed_contours->show_all_actor();
  m_render_window->Render();
}

void CShapeDisplay::key_bracketleft()
{
  --current_lineID; 
  if (current_lineID<0)
  {
    current_lineID=0;
  }
  m_correspond->update_actors(current_LayerID,current_lineID);
  m_render_window->Render();
}

void CShapeDisplay::key_bracketright()
{
  ++current_lineID;
  m_correspond->update_actors(current_LayerID,current_lineID);
  m_render_window->Render();

    
}

void CShapeDisplay::key_o()
{
  m_correspond->solo_corespond_display(current_correspondID++);
  if (current_correspondID>m_shape->vec_Cor.size())
  {
    current_correspondID=m_shape->vec_Cor.size();
  }
  m_render_window->Render();
}

void CShapeDisplay::key_p()
{
  m_correspond->solo_corespond_display(current_correspondID--);
  if (current_correspondID<0)
  {
    current_correspondID=0;
  }
  m_render_window->Render();
}
extern vtkSmartPointer<vtkRenderer> Renderers[5];
extern double xmins[5];
extern double xmaxs[5];
extern double ymins[5];
extern double ymaxs[5];

void CShapeDisplay::key_change_viewport(std::string key)
{
  update_viewport(key);
  m_render_window->Render();
}

void CShapeDisplay::update_viewport(std::string key)
{
  int N=5;
  int index= atoi(key.c_str());
  if (index!=0)
  {
    index=index -1;
    for (int i = 0; i < N; ++i)
    {
      if (i == index)
      {
        Renderers[i]->SetViewport(0,0,1,1);
      }
      else
      {
        Renderers[i]->SetViewport(0,0,0,0);
      }
    }
  }
  else if(index==0)
  {
    for (int i = 0; i < N; ++i)
    {
      Renderers[i]->SetViewport(xmins[i],ymins[i],xmaxs[i],ymaxs[i]);
    }
  }
}
