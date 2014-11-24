#include"CContour.h"
#include"CPoint.h"
#include "assert.h"
#include "CLayer.h"
#include "math.h"
using namespace std;

float  CContour::static_max_x=-999999;
float  CContour::static_min_x=999999;
float  CContour::static_max_y=-999999;
float  CContour::static_min_y=999999;
int CContour::static_contour_index=0;
CContour::CContour(const float ID,CLayer* layer)
{
  m_layer=layer;
  contourID=ID;
  filename="";
  max_x=-9999999;
  max_y=-9999999;
  min_x=9999999;
  min_y=9999999;
  original_max_x=-9999999;
  original_max_y=-9999999;
  original_min_x=9999999;
  original_min_y=9999999;
  center_point= new CPoint();
  moment_one_point= new CPoint();
  length=0;original_length=0;
  last_x=last_y=9999999;
  original_last_x=original_last_y=9999999;
  //  use_as_higher_contour_count=0;
  //  use_as_lower_contour_count=0;
  use_counter=0;
  sum_x=0;
  sum_y=0;
  m_Map= new CMap(this);

  contour_index=static_contour_index++;
}

void CContour:: operator=(CContour &temp)
{
  assert(false);
  contourID = temp.contourID;
  m_Map = temp.m_Map;
  sum_x=sum_y=0;
}
#include "stdio.h"
#include "iostream"
void CContour:: check_edge(float tempx,float tempy )
{
    if (tempx>max_x)
    {
      max_x=tempx;
    }
    if (tempx<min_x)
    {
      min_x=tempx;
    }
    sum_x+=tempx;

    if (tempy>max_y)
    {
      max_y=tempy;
    }
    if(tempy<min_y)
    {
      min_y=tempy;
    }
    sum_y+=tempy;
    if (last_x==9999999&&last_y==9999999)
    {
      last_x=tempx;last_y=tempy;
    }
    float temp_value=(tempx-last_x)*(tempx-last_x)+(tempy-last_y)*(tempy-last_y);
    length+=sqrt(temp_value);
    // std::cout<<length<<"\t"<<tempx<<"\t"<<last_x<<"\t"<<tempy<<"\t"<<last_y<<"\n";
    last_x=tempx;last_y=tempy;
}

void CContour::check_edge_before_scale(float tempx,float tempy)
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

    if (original_last_x==9999999&&original_last_y==9999999)
    {
      original_last_x=tempx;original_last_y=tempy;
    }
    float temp_value=(tempx-original_last_x)*(tempx-original_last_x)+(tempy-original_last_y)*(tempy-original_last_y);
    original_length+=sqrt(temp_value);

    original_last_x=tempx;original_last_y=tempy;
}

void CContour::get_center()
{
  center_point->x=(max_x+min_x)/2;
  center_point->y=(max_y+min_y)/2;
  center_point->z=LayerID;//because when do the normalize, we don't want change the z level;
}

void CContour::calculate_one_moment()
{
  moment_one_point->x=sum_x/(PointNum*1.0);
  moment_one_point->y=sum_y/(PointNum*1.0);
  moment_one_point->z=LayerID;
}
bool CContour::read_single_layer_without_z(fstream& fin)
{
  int oldz=0,newz=-9999;

  float tempvalue;
  while(fin.good())
  {
    CPoint *temp=new CPoint();
    fin>>temp->x;
    fin>>temp->y;
    check_edge(temp->x,temp->y);
    temp->z=LayerID*10;//TODO
    if(fin.good())
    {
      vec_points.push_back(temp);
    }
  }
  PointNum=vec_points.size();
  get_center();
  calculate_one_moment();
  return true;
}


bool CContour::read_contour_with_z(fstream& fin)
{
  float levelID=-1;
  int temp_contourID=-1;
  int pos=-1;
  if(fin.good())
    {
      pos=fin.tellg();
      fin>>levelID;
      fin>>temp_contourID;
      if (temp_contourID!=contourID)
      {
        return false;
      }
      LayerID=levelID;
      fin.seekg(pos);
    }
  else
    {
      return false;
    }
  int count=0;
  while(fin.good())
    {
      count++;
      CPoint *temp=new CPoint();
      int pos=fin.tellg();
      fin>>temp->z;
      fin>>temp_contourID;
      fin>>temp->x;
      fin>>temp->y;
      check_edge_before_scale(temp->x,temp->y);
      //      cout<<temp->x<<"\t"<<temp->y<<"\t"<<temp->z<<"\n";
      temp->x=temp->x*point_scale;
      temp->y=temp->y*point_scale;
      //      cout<<temp->x<<"\t"<<temp->y<<"\t"<<temp->z<<"\n";
    

      if(temp_contourID!=contourID||temp->z!=LayerID)
	{
          fin.seekg(pos);
	  delete temp;
	  break;
	}
      else if(count%interval_of_point==0)
	{
          temp->z=temp->z*layer_interval_scale;
          temp->index=temp->get_index();
	  vec_points.push_back(temp);
          check_edge(temp->x,temp->y);
	}
      if(fin.eof())
      {
        break;
      }
    }
  PointNum=vec_points.size();
  check_point_seq();
  get_center();
  calculate_one_moment();
  return true;
}

void CContour::check_point_seq()
{
  std::vector<CPoint*>::iterator itr,etr,nitr;
  itr=vec_points.begin();etr=vec_points.end();
  float avarge = length/PointNum;
  (*itr)->used_times++;
  vec_points_orderd.push_back(*itr);
  for (; itr!=etr; ++itr)
  {
    nitr=itr;
    nitr++;
    if (nitr!=etr)
    {
      float dis = ((*itr)->x-(*nitr)->x)*((*itr)->x-(*nitr)->x) +  ((*itr)->y-(*nitr)->y)*((*itr)->y-(*nitr)->y);
      dis = sqrt(dis);
      if (dis>point_distance_scale*avarge&&point_distance_scale!=0)
      {
        CPoint* point=get_the_nearest_unused_point(*itr);
        point->used_times++;
        vec_points_orderd.push_back(point);
      }
      else
      {
        (*nitr)->used_times++;
        vec_points_orderd.push_back(*nitr);
      }
    }
  }
}

CPoint*CContour:: get_the_nearest_unused_point(CPoint* point)
{
  std::vector<CPoint*>::iterator itr,etr,nitr;
  itr=vec_points.begin();etr=vec_points.end();
  float max=99999999;
  CPoint* first_nearest=NULL,*second_nearest=NULL;
  for (;itr!=etr ; ++itr)
  {
    if (point!=*itr)
    {
      float dis= ((*itr)->x-point->x) * ((*itr)->x-point->x) + ((*itr)->y-point->y)*((*itr)->y-point->y);
      dis=sqrt(dis);
      if (dis<max)
      {
        if ((*itr)->used_times==0)
        {
          first_nearest=*itr;
          max=dis;
        }
      }
    }
  }
  first_nearest->used_times++;
  return first_nearest;
}

#include "CFileDebug.h"
void CContour::InitMap()
{
  // normalize();
  m_Map= new CMap(this);
  m_Map->setup();
  // if(LayerID==1&&contourID==1)
  // {
  //   CFileDebug m_file2("synthetic3_1_1_dismap");
  //   m_file2.Output(m_Map->DistancsMap);
  // }
  m_Map->gradient();
}

void CContour::normalize(CPoint shape_center_point)//may be use the moment point
{
  std::vector<CPoint*>::iterator itr,etr;
  itr=vec_points_orderd.begin();
  etr=vec_points_orderd.end();
  int n=99999999;

  int map_center_x=NumRows/2;
  int map_center_y=NumCols/2;
  for (;itr!=etr;++itr)
  {
    CPoint* temp = new CPoint();
    temp->x=(*itr)->x-shape_center_point.x + map_center_x;
    temp->y=(*itr)->y-shape_center_point.y + map_center_y;
    temp->z=(*itr)->z;
    // if (temp->x>NumRows||temp->y>NumRows)
    // {
    //   cout<<temp->x<<"\t:"<<temp->y<<"\n";
    // }

    temp->index=temp->get_index();
    vec_Points_Origin.push_back(temp);
    if (temp->x>static_max_x)
    {
      static_max_x=temp->x;
    }
    if (temp->x<static_min_x )
    {
      static_min_x=temp->x;
    }
    if (temp->y>static_max_y)
    {
      static_max_y=temp->y;
    }
    if (temp->y<static_min_y)
    {
      static_min_y=temp->y;
    }
  }
  // center_point->x=center_point->x-shape_center_point.x+map_center_x;
  // center_point->y=center_point->y-shape_center_point.y+map_center_y;
  // center_point->z=LayerID;

  // moment_one_point->x=moment_one_point->x-shape_center_point.x+map_center_x;
  // moment_one_point->y=moment_one_point->y-shape_center_point.y+map_center_y;
  // moment_one_point->z=LayerID;

}

void CContour::reset()
{
  assert(false);
  //use_as_lower_contour_count=0;
  //use_as_higher_contour_count=0;
}

void CContour::smooth()
{
  if (use_contour_smooth==false)
  {
    return;
  }
  float smooth_factor=30.0;
  int size=vec_Points_Origin.size();
  for (int i = 0; i < size; ++i)
  {
    CPoint temp1 = *vec_Points_Origin[i];
    for(int j=-smooth_factor;j <= smooth_factor;j++ )
    {
      int index=(i+j)%size;
      if (index<0)
      {
        index=index+size;
      }
      CPoint temp2 = *vec_Points_Origin[index];
      temp1 = (temp1 + temp2);
    }
    CPoint temp = temp1/(smooth_factor*2+2);
    *vec_Points_Origin[i]= temp;
  }
  
}

// #include"CFileDebug.h"
// void CContour::calculate_medial_map(float** medial_axis)
// {
//   m_medial_map = new CMedialMap(this,medial_axis);
//   m_medial_map->setup();
//   m_medial_map->gradient();
//   //CFileDebug m_debugfile6("./medial_DistancsMap_20");
//   //m_debugfile6.Output(m_medial_map->DistancsMap);
//   //CFileDebug file("./only_medial_axis");
//   //file.Output(medial_axis);
//   // CFileDebug m_debugfile7("./medial_sign_20");
//   //  m_debugfile7.Output_sign(m_medial_map->SignMap);
  


// }

// void CContour::swap_map_medialmap()
// {
//   assert(m_Map);
//   assert(m_medial_map);
//   m_temp_map=m_Map;
//   m_Map=(CMap*)m_medial_map;
// }

// void CContour::swap_medialmap_map()
// {
//   assert(m_Map);
//   assert(m_medial_map);
//   m_medial_map=(CMedialMap*)m_Map;
//   m_Map=m_temp_map;
//   delete m_medial_map;
// }



