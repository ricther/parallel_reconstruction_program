#include"CShape.h"
#include"CLayer.h"
#include"CRegistration.h"
#include"CMap.h"
#include"CPoint.h"
#include"CSkeleton.h"
#include<omp.h>
#include "CContourArrangement.h"
using namespace std;
/** 
 * initial m_skeleton;m_registration
 * 
 */
CShape::CShape()
{
  float n=999999;
  max_x=max_y=max_z=-n;
  min_x=min_y=min_z=n;
  sum_x=sum_y=0;
  m_skeleton=new CSkeleton();
  m_registration= new CRegistration(this);
  original_max_x=-9999999;
  original_max_y=-9999999;
  original_min_x=9999999;
  original_min_y=9999999;
  contour_arrangement= new CContourArrangement(this);
  m_total_points = vtkSmartPointer<vtkPoints>::New();
  vtk_points_counter=0;
  map_CPointsIndex_vtkIndex.clear();
  total_points_num=0;
}

/** 
 * intial the shape, invoke the read()
 * 
 * @param filename 
 */
void CShape::initial(string filename)
{
  m_filename=filename;
  read();   
}

/** 
 * 
 * 
 * @param temp_window 
 * @param temp_interactor 
 */
void CShape::initial_display(vtkSmartPointer<vtkRenderWindow> temp_window,vtkSmartPointer<vtkRenderWindowInteractor> temp_interactor)
{
  m_display= new CShapeDisplay(this,temp_window,temp_interactor);
}

#include<fstream>
/** 
 * read the image data, and invoke the layer.read() and contour.read()
 * fill the map_layer
 */
void CShape::read()
{
  //read the structure like that: level contourID x y
    fstream fin(m_filename.c_str(),ios_base::in);
  if(fin.is_open())
  {
    bool result=true;
    float levelID=-1;
    while(result&&fin.good())
    {
      int pos=fin.tellg();
      fin>>levelID;
      fin.seekg(pos);
      
      CLayer * temp_layer;
      if (map_Layer.find(levelID)==map_Layer.end())
      {
        temp_layer=new CLayer(levelID);
        result=temp_layer->read_layer_multi_contour_with_z(fin);
        if (temp_layer->contour_Num<=0)
        {
          cout<<"layer:"<<temp_layer->LayerID<<"deleted"<<endl;
          delete temp_layer;
          continue;
        }
        //        check_edge(temp_layer->center_point->x,temp_layer->center_point->y,temp_layer->center_point->z);
        check_edge(temp_layer->max_x,temp_layer->max_y,temp_layer->center_point->z);
        check_edge(temp_layer->min_x,temp_layer->min_y,temp_layer->center_point->z);
        check_edge_before_scale(temp_layer->original_max_x,temp_layer->original_max_y);
        check_edge_before_scale(temp_layer->original_min_x,temp_layer->original_min_y);

        map_Layer.insert(make_pair(temp_layer->LayerID,temp_layer));
      }
      else
      {
        temp_layer=map_Layer.find(levelID)->second;
        result=temp_layer->read_layer_multi_contour_with_z(fin);
        //        check_edge(temp_layer->center_point->x,temp_layer->center_point->y,temp_layer->center_point->z);
        check_edge(temp_layer->max_x,temp_layer->max_y,temp_layer->center_point->z);
        check_edge(temp_layer->min_x,temp_layer->min_y,temp_layer->center_point->z);

        check_edge_before_scale(temp_layer->original_max_x,temp_layer->original_max_y);
        check_edge_before_scale(temp_layer->original_min_x,temp_layer->original_min_y);

      }

    }
  }
  get_center();
  calculate_one_moment();
  check_point_scale();
}

void CShape::check_point_scale()
{
  float x= original_max_x-original_min_x;
  float y= original_max_y-original_min_y;
  float v= x>y?x:y;
  float scale = NumRows/v;
  cout<<"\n********************************************************************\nthe longest distance in x direction is:" <<x<<"\t longest distance in y direction is:" <<y <<" the specified point scale is : "<<point_scale <<"\t the recommend point scale is:"<<scale<<"\n********************************************************************\n";
}
/** 
 * build skeleton;
 * normalize the shape position;
 * 
 */
void CShape::Setup()
{
  //  m_skeleton->build_skeleton_new(map_Layer);
  std::map<float,CLayer*>::iterator itr=map_Layer.begin(),etr=map_Layer.end();
  int layer_count=0;
  int contour_count=0;
  for(;itr!=etr;++itr)
  {
    vec_layerID.push_back(itr->first);
    layer_count++;
    (itr->second)->setup(center_point);
    contour_count+=(itr->second)->map_contour.size();
  }
    cout<<"\n the max and min point value after scale. max_x ="<<CContour::static_max_x<<"\t min_x="<< CContour::static_min_x<<"\t max_y="<< CContour::static_max_y<<"\t min_y="<< CContour::static_min_y<<"\n";

  std::cout<<"**********Layer NUM:"<<layer_count<<"\t"<<"Contour NUM:"<<contour_count<<"\n";
  if (write_shape_points_to_data)
  {
      write_vtk_points();
      std::cout<<"**********Write points data to vtk_points.vtk, total_points_number:"<<total_points_num<<"***************\n";

  }
}

void CShape::Setup_use_openmp()
{
  // m_skeleton->build_skeleton_new(map_Layer);
  std::map<float,CLayer*>::iterator itr=map_Layer.begin(),etr=map_Layer.end();
  for(;itr!=etr;++itr)
  {
    vec_layerID.push_back(itr->first);
  }

  int size=vec_layerID.size();
  int layer_count=0;
  int contour_count=0;
#pragma omp parallel for num_threads(thread_num)
  for(int i=0;i<size;++i)
  {
    layer_count++;
    map_Layer[vec_layerID[i]]->setup(moment_one_point);
    contour_count+=map_Layer[vec_layerID[i]]->map_contour.size();
  }

    cout<<"\n the max and min point value after scale. max_x ="<<CContour::static_max_x<<"\t min_x="<< CContour::static_min_x<<"\t max_y="<< CContour::static_max_y<<"\t min_y="<< CContour::static_min_y<<"\n";

  std::cout<<"**********Layer NUM:"<<layer_count<<"\t"<<"Contour NUM:"<<contour_count<<"\n";
}


void CShape::Registration()
{
  if (use_open_mp)
  {
    m_registration->Register_use_openmp();
  }
  else
  {
    m_registration->Register();
  }

}

void CShape:: check_edge(float tempx,float tempy,float tempz)
{

    if (tempx>max_x)
    {
      max_x=tempx;
    }
    if (tempx<min_x)
    {
      min_x=tempx;
    }



    if (tempy>max_y)
    {
      max_y=tempy;
    }
    if(tempy<min_y)
    {
      min_y=tempy;
    }



    if (tempz>max_z)
    {
      max_z=tempz;
    }
    if (tempz<min_z)
    {
      min_z=tempz;
    }

    sum_x+=tempx;sum_y+=tempy;
}

void CShape::check_edge_before_scale(float tempx,float tempy)
{
    if (tempx>original_max_x)
    {
      original_max_x=tempx;
    }
    if (tempx<original_min_x)
    {
      original_min_x=tempx;
    }

    if (tempy>original_max_y)
    {
      original_max_y=tempy;
    }
    if(tempy<original_min_y)
    {
      original_min_y=tempy;
    }

}



void CShape::get_center()
{
  center_point.x=(max_x+min_x)/2;
  center_point.y=(max_y+min_y)/2;
  center_point.z=0;//because when do the normalize, we don't want change the z level;
  cout<<"max_x"<<max_x<<"\t min_x"<<min_x<<"\t max_y"<<max_y<<"\t min_y"<<min_y<<"center_point.x:"<<center_point.x<<"center_point.y:"<<center_point.y<<"\n";
}

void CShape::calculate_one_moment()
{
  moment_one_point.x=sum_x/map_Layer.size();
  moment_one_point.y=sum_y/map_Layer.size();
}
//direction = 0 down 1 up
float CShape::get_next_layer(float now_layerID,int direction)
{
  std::map<float,CLayer*>::iterator itr = map_Layer.find(now_layerID);
  if (direction==0)
  {
    if (itr !=map_Layer.begin())
    {
      itr--;
      return itr->first;
    }
    else
    {
      return itr->first;
    }
  }
  else
  {
    std::map<float,CLayer*>::iterator nitr=itr;
    nitr++;
    if (nitr != map_Layer.end())
    {
      return nitr->first;
    }
    else
    {
      return itr->first;
    }
  }
}

void CShape::initial_vtk_points()
{
  /// intial the points in the CContour
  std::map<float,CLayer*>::iterator itr_l=map_Layer.begin(),etr_l=map_Layer.end();
  std::map<int,CContour*>::iterator itr_c,etr_c;
  for (;itr_l!=etr_l; ++itr_l)
  {
    CLayer* temp_layer= itr_l->second;
    itr_c=temp_layer->map_contour.begin();
    etr_c=temp_layer->map_contour.end();
    for (; itr_c!=etr_c; ++itr_c)
    {
      CContour* temp_contour= itr_c->second;
      insert_vtk_points(temp_contour->vec_points);
      insert_vtk_points(temp_contour->vec_Points_Origin);
    }
  }
  /// intial the points in the CVirtual

  int size = contour_arrangement->vec_contour_couple.size();
  for (int i = 0; i < size; ++i)
  {
    CVirtualContour* vcontour1 = contour_arrangement->vec_contour_couple[i]->vcontour1;
    CVirtualContour* vcontour2 = contour_arrangement->vec_contour_couple[i]->vcontour2;

    insert_vtk_points(vcontour1->vec_Points_project);
    insert_vtk_points(vcontour1->vec_medial_points);

    insert_vtk_points(vcontour2->vec_Points_project);
    insert_vtk_points(vcontour2->vec_medial_points);

  }
  
}

void CShape:: insert_vtk_points(std::vector<CPoint*>& vec_points)
{
  if (vec_points.size()==0)
  {
    return;
  }
  std::vector<CPoint*>::iterator itr_p, etr_p;
  itr_p=vec_points.begin();
  etr_p=vec_points.end();
  for (; itr_p!=etr_p; ++itr_p)
  {
    insert_map_CPointsIndex_vtkIndex((*itr_p),get_vtk_points_index());
  }
}

#include "vtkPolyDataWriter.h"
#include "vtkNew.h"
void CShape::write_vtk_points()
{
  vtkSmartPointer<vtkPoints> m_vtk_points= vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> m_vertices= vtkSmartPointer<vtkCellArray>::New();

  std::map<float,CLayer*>::iterator itr_l=map_Layer.begin(),etr_l=map_Layer.end();
  std::map<int,CContour*>::iterator itr_c,etr_c;

  for (;itr_l!=etr_l; ++itr_l)
  {
    CLayer* temp_layer= itr_l->second;
    itr_c=temp_layer->map_contour.begin();
    etr_c=temp_layer->map_contour.end();
    for (; itr_c!=etr_c; ++itr_c)
    {
      CContour* temp_contour= itr_c->second;
      total_points_num+=temp_contour->vec_points.size();
    }
  }

  itr_l=map_Layer.begin(),etr_l=map_Layer.end();

  vtkIdType* Indices = new vtkIdType[total_points_num];
  int counter=0;
  for (;itr_l!=etr_l; ++itr_l)
  {
    CLayer* temp_layer= itr_l->second;
    itr_c=temp_layer->map_contour.begin();
    etr_c=temp_layer->map_contour.end();
    for (; itr_c!=etr_c; ++itr_c)
    {
      CContour* temp_contour= itr_c->second;
      std::vector<CPoint*> vec_points=temp_contour->vec_points;
      std::vector<CPoint*>::iterator itr_p, etr_p;
      itr_p=vec_points.begin();
      etr_p=vec_points.end();
      for (; itr_p!=etr_p; ++itr_p)
      {
        m_vtk_points->InsertPoint(static_cast<vtkIdType>(counter),(*itr_p)->x,(*itr_p)->y,(*itr_p)->z);
        Indices[counter]=static_cast<vtkIdType>(counter);
        ++counter;
      }
    }
  }
  m_vertices->InsertNextCell(total_points_num,Indices);
  vtkSmartPointer<vtkPolyData> m_polydata= vtkSmartPointer<vtkPolyData>::New();
  m_polydata->SetPoints(m_vtk_points);
  m_polydata->SetVerts(m_vertices);
  vtkNew<vtkPolyDataWriter> writer;
  writer->SetFileName("./points.vtk");
  writer->SetInput(m_polydata);
  writer->Write();

}
