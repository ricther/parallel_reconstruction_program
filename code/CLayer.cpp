#include"CLayer.h"
#include"CPoint.h"
#include "CMap.h"
#include "CContour.h"
#include "assert.h"
#include "Memory.h"
#include "stdio.h"
#include "iostream"
using namespace std;


int CLayer::s_layer_index=0;
CLayer::CLayer(const float ID)
{
  LayerID=ID;
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
  contour_Num=0;
  total_length=0;
  layer_index=s_layer_index++;
}

void CLayer:: operator=(CLayer &temp)
{
  LayerID = temp.LayerID;
  map_contourID = temp.map_contourID;
  map_contour=temp.map_contour;
  //just want a pause if call this function
  assert(false);
  sum_x=sum_y=0;
}

void CLayer:: check_length(float length)
{
  total_length+=length;
}
void CLayer:: check_edge(float tempx,float tempy )
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
}

void CLayer::check_edge_before_scale(float tempx,float tempy)
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


void CLayer::get_center()
{
  center_point->x=(max_x+min_x)/2;
  center_point->y=(max_y+min_y)/2;
  center_point->z=LayerID;//because when do the normalize, we don't want change the z level;
}

void CLayer::calculate_one_moment()
{
  moment_one_point->x=sum_x/(map_contour.size()*1.0);
  moment_one_point->y=sum_y/(map_contour.size()*1.0);
  moment_one_point->z=LayerID;
}

bool CLayer:: read_layer_multi_contour_with_z(std::fstream& fin)
{

  int pos=-1;
  float levelID;
  int contourID;
  bool result=true;

  if (fin.good())
  {
    pos=fin.tellg();
    fin>>levelID;
    fin>>contourID;
    fin.seekg(pos);
  }
  else
  {
    return false;
  }
  
  if (fin.good()&&LayerID!=levelID)
  {
    assert(false);
    return true;// buffer has content need read more
  }


  while(fin.good()&&result)
  {
    
    if (fin.good())
    {
      pos=fin.tellg();
      fin>>levelID;
      fin>>contourID;
      fin.seekg(pos);
    }
    else
    {
      return false;
    }
    
    if (levelID!=LayerID)
    {
      return true;
    }

    CContour* new_contour;
    if (map_contour.find(contourID)==map_contour.end())//handle the disorder situiation
    {
      new_contour= new CContour(contourID,this);
      result=new_contour->read_contour_with_z(fin);
      cout<<"level:"<<new_contour->LayerID<<"\t length:"<<new_contour->length<<"\t point num:"<<new_contour->PointNum;
      if (new_contour->length<contour_length_threshold&&new_contour->PointNum<contour_pointnum_threshold)
      {
        cout<<"deleted"<<endl;
        delete new_contour;
        continue;
      }else
      {
        cout<<"\n";
      }


      map_contourID.insert(make_pair(contour_Num,contourID));
      map_contour.insert(make_pair(contourID,new_contour));
      check_edge(new_contour->max_x,new_contour->max_y);
      check_edge(new_contour->min_x,new_contour->min_y);
      check_edge_before_scale(new_contour->original_max_x,new_contour->original_max_y);
      check_edge_before_scale(new_contour->original_min_x,new_contour->original_min_y);
      get_center();
      check_length(new_contour->length);
      calculate_one_moment();
      ++contour_Num;
    }
    else
    {
      new_contour=map_contour[contourID];
      result=new_contour->read_contour_with_z(fin);
      check_edge(new_contour->max_x,new_contour->max_y);
      check_edge(new_contour->min_x,new_contour->min_y);
      check_edge_before_scale(new_contour->original_max_x,new_contour->original_max_y);
      check_edge_before_scale(new_contour->original_min_x,new_contour->original_min_y);
      get_center();
      calculate_one_moment();
    }
  }
  return true;
}

void CLayer::setup(CPoint shape_center_point)
{
  std::map<int,CContour*>::iterator itr,etr;
  itr=map_contour.begin();etr=map_contour.end();
  for (; itr!=etr; ++itr)
  {
    (itr->second)->normalize(shape_center_point);
    (itr->second)->smooth();
  }
}

void CLayer::reset()
{
  std::map<int,CContour*>::iterator itr,etr;
  itr=map_contour.begin();etr=map_contour.end();
  for (; itr!=etr; ++itr)
  {
    (itr->second)->reset();
  }
}
