#include"CShape.h"
#include"CLayer.h"
#include"CRegistration.h"
#include"CMap.h"
#include"CPoint.h"
#include"CSkeleton.h"
using namespace std;
CShape::CShape()
{
  float n=999999;
  max_x=max_y=max_z=-n;
  min_x=min_y=min_z=n;
  sum_x=sum_y=0;
  m_skeleton=new CSkeleton();
}
void CShape::initial(vector<string>&t_vec_filename)
{
  m_registration= new CRegistration(this);
  vec_filename=t_vec_filename;
  int size=vec_filename.size();
  if (size>1)
  {
    // for (int i = 0; i < size; ++i)
    // {
    //   string filename=vec_filename[i];
    //   fstream fin(vec_filename[i].c_str(),ios_base::in);
    //   if (fin.is_open())
    //   {
    //     CLayer * newlayer = new CLayer(i+1);//in order to keep consistent with the single file,the first level is 1, not zero
    //     newlayer->filename=filename;
    //     newlayer->read_layer_single_contour_without_z(fin);
    //     check_edge(newlayer->center_point->x,newlayer->center_point->y,newlayer->center_point->z);
    //     map_Layer.insert(make_pair(newlayer->LayerID,newlayer));
    //   }
    // }
    // get_center();
    // calculate_one_moment();
  }
  else if(size==1)
  {
    read();   
  }
 
}

void CShape::initial_display(vtkSmartPointer<vtkRenderWindow> temp_window,vtkSmartPointer<vtkRenderWindowInteractor> temp_interactor)
{
  m_display= new CShapeDisplay(this,temp_window,temp_interactor);
}

#include<fstream>
void CShape::read()
{
  //read the structure like that: level contourID x y
  string filename = vec_filename[0];
  fstream fin(filename.c_str(),ios_base::in);
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
        check_edge(temp_layer->center_point->x,temp_layer->center_point->y,temp_layer->center_point->z);
        map_Layer.insert(make_pair(temp_layer->LayerID,temp_layer));
      }
      else
      {
        temp_layer=map_Layer.find(levelID)->second;
        result=temp_layer->read_layer_multi_contour_with_z(fin);
        check_edge(temp_layer->center_point->x,temp_layer->center_point->y,temp_layer->center_point->z);
      }

    }
  }
  get_center();
  calculate_one_moment();
}

void CShape::Setup()
{
  m_skeleton->build_skeleton_new(map_Layer);
  std::map<float,CLayer*>::iterator itr=map_Layer.begin(),etr=map_Layer.end();
  int layer_count=0;
  int contour_count=0;
  for(;itr!=etr;++itr)
  {
    layer_count++;
    (itr->second)->setup(moment_one_point);
    contour_count+=(itr->second)->map_contour.size();
    //    (itr->second)->normalize(center_point);
    //    (itr->second)->normalize(moment_one_point);
    //    (itr->second)->InitMap();
    //(itr->second)->m_Map->gradient();
  }
  std::cout<<"**********Layer NUM:"<<layer_count<<"\t"<<"Contour NUM"<<contour_count<<"\n";
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

void CShape::get_center()
{
  center_point.x=(max_x+min_x)/2;
  center_point.y=(max_y+min_y)/2;
  center_point.z=0;//because when do the normalize, we don't want change the z level;
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
