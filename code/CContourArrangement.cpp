#include "CContourArrangement.h"
#include "math.h"
#include "map"
#include "CPoint.h"
#include "assert.h"
#include "stdio.h"
#include "iostream"

contour_couple::contour_couple()
{
  contour1=NULL;contour2=NULL;
  b_intersection=false;
  distance=99999999999;
  boundary_coverage_rate=0.0;
  b_couple=false;
}

bool contour_couple:: operator == (const contour_couple& right)
{
  if(this->contour1==right.contour1)
  {
    if(this->contour2==right.contour2)
    {
      return true;
    }
  }
  if (this->contour1==right.contour2)
  {
    if (this->contour2==right.contour1)
    {
      return true;
    }
  }
  return false;
}


void CContourArrangement::setup()
{
  std::vector<contour_couple*> vec_couple1,vec_couple2;
  get_couple(up_layer,down_layer,vec_couple1);
  get_couple(down_layer,up_layer,vec_couple2);
  get_union_couple(vec_couple1,vec_couple2);
  vec_couple1.clear();vec_couple2.clear();
}

void CContourArrangement::get_union_couple(  std::vector<contour_couple*>& couple1,  std::vector<contour_couple*>& couple2)
{
  vec_contour_couple=couple1;
  std::vector<contour_couple*>::iterator itr=couple2.begin(),etr=couple2.end();
  for (; itr!=etr; ++itr)
  {
    std::vector<contour_couple*>::iterator itr2=vec_contour_couple.begin(),etr2=vec_contour_couple.end();
    bool find_same=false;
    for (; itr2 !=etr2; ++itr2)
    {
      if(*(*itr)==*(*itr2))
      {
        find_same=true;
        break;
      }
    }
    if (find_same==false)
    {
      vec_contour_couple.push_back(*itr);
    }
  }
}
//make couple depend on the intersection rate and distance, did't decided the compute arrangement.
void CContourArrangement:: get_couple(CLayer* first_layer, CLayer* second_layer,  std::vector<contour_couple*>& vec_couple )
{
    std::map<int,CContour*>::iterator u_itr,u_etr,d_itr,d_etr;
    u_itr=first_layer->map_contour.begin();u_etr=first_layer->map_contour.end();

  float intersection_rate;
  float distance;
  CContour* temp1,*temp2;
  for (;u_itr!=u_etr ; u_itr++)
  {
    contour_couple* new_couple=new contour_couple();
    new_couple->contour1=u_itr->second;
    d_itr=second_layer->map_contour.begin();d_etr=second_layer->map_contour.end();
    for (;d_itr!=d_etr ;d_itr++ )
    {
      temp1=u_itr->second;temp2=d_itr->second;
      intersection_rate=check_intersection_level_set(temp1,temp2);
      distance=cal_distance(u_itr->second,d_itr->second);
      if (intersection_rate>boundary_coverage_rate_threshold)
      {
        if (intersection_rate>new_couple->boundary_coverage_rate)
        {
          new_couple->boundary_coverage_rate=intersection_rate;
          new_couple->contour2=d_itr->second;
          new_couple->b_couple=true;
        }
        // if (distance<new_couple->distance)
        // {
        //   new_couple->distance=distance;
        //   new_couple->contour2=d_itr->second;
        //   new_couple->b_couple=true;
        // }
      }
      else if (intersection_rate<boundary_coverage_rate_threshold)
      {
        if (distance<contour_distance_threshold)
        {
          if (distance<new_couple->distance)
          {
            std::cout<<first_layer->LayerID<<"-"<<second_layer->LayerID<<"distance:"<<distance<<std::endl;
            new_couple->distance=distance;
            new_couple->contour2=d_itr->second;
            new_couple->b_couple=true;
          }
        }
      }
    }
    if (new_couple->b_couple==false)
    {
      delete new_couple;
    }
    else
    {
      vec_couple.push_back(new_couple);
    }
  }
}

float CContourArrangement::check_intersection(CContour* u, CContour*d)
{
  if (u->max_x>=d->min_x&&u->min_x<=d->max_x)
  {
    if (u->max_y>=d->min_y&&u->min_y<=d->max_y)
    {
     float rate= check_boundary_box_coverage(u,d);
     return rate;
    }
    else
    {
      return 0;
    }
  }
  return 0;

  // if (u->max_x>d->min_x&&u->min_x<d->max_x)
  // {
  //   if (u->max_y>d->min_y&&u->min_y<d->max_y)
  //   {
  //     return true;
  //   }
  //   else
  //   {
  //     return false;
  //   }
  // }
  // return false;
}

float CContourArrangement::check_boundary_box_coverage(CContour* u, CContour*d)
{
  float dis_x=  check_boundary_box_distance(u->min_x,u->max_x,d->min_x,d->max_x);
  float dis_y=  check_boundary_box_distance(u->min_y,u->max_y,d->min_y,d->max_y);
  float area_intersection=dis_x*dis_y;
  float area_u=(u->max_x-u->min_x)*(u->max_y-u->min_y);
  float area_d=(d->max_x-d->min_x)*(d->max_y-d->min_y);
  float rate_u=area_intersection/area_u;
  float rate_d=area_intersection/area_d;
  float rate=rate_u*rate_d;
  std::cout<<u->LayerID<<"-"<<d->LayerID<<"rate:"<<rate<<std::endl;
  //  float rate=rate_u>rate_d?rate_u:rate_d;
  return rate;
}

float  CContourArrangement::check_boundary_box_distance(float num1,float num2,float num3, float num4)//just calculate the sequnce, num1<num2 , num3 < num4
{
  float seq1,seq2,seq3,seq4,temp;
  if (num1<=num3)
  {
    seq1=num1;
    if (num2>=num4)
    {
      seq2=num3;seq3=num4;seq4=num2;
    }
    else if (num2<=num4&&num3<=num2)
    {
      seq2=num3;seq3=num2;seq4=num4;
    }
    else if (num3>=num2)//in this case no intersection
    {
      assert(false);
    }
  }
  else if (num1>num3)
  {
    seq1=num3;
    if (num2<=num4)
    {
      seq2=num1;seq3=num2;seq4=num4;
    }
    else if (num2>=num4&&num1<=num4)
    {
      seq2=num1;seq3=num4;seq4=num2;
    }
    else if (num1>=num4)
    {
      assert(false);//no intersection
    }
  }
  return seq3-seq2;
}

float CContourArrangement::cal_distance(CContour* u, CContour*d)
{
  float dis=(u->moment_one_point->x-d->moment_one_point->x)*(u->moment_one_point->x-d->moment_one_point->x)+(u->moment_one_point->y-d->moment_one_point->y)*(u->moment_one_point->y-d->moment_one_point->y);
  dis=sqrt(dis);
  return dis;
}

float CContourArrangement::check_intersection_level_set(CContour* u, CContour*d)
{
  if (u->max_x>=d->min_x&&u->min_x<=d->max_x)
  {
    if (u->max_y>=d->min_y&&u->min_y<=d->max_y)
    {
     float rate= check_coverage_level_set(u,d);
     return rate;
    }
    else
    {
      return 0;
    }
  }
  return 0;

}

float CContourArrangement::check_coverage_level_set(CContour* u, CContour*d)
{
  float area_u=u->m_Map->area;
  float area_d=d->m_Map->area;
  float area_intersection=check_intersection_use_level_set(u,d);
  if (area_intersection==0)
  {
    return 0;
  }
  float rate_u=area_intersection/area_u;
  float rate_d=area_intersection/area_d;
  float rate=rate_u*rate_d;
  std::cout<<u->LayerID<<"-"<<d->LayerID<<"rate:"<<rate<<std::endl;
  //  float rate=rate_u>rate_d?rate_u:rate_d;
  return rate;
  
}

float CContourArrangement::check_intersection_use_level_set(CContour* u, CContour* d)
{
  int temp=0;
  for (int i = 0; i < NumRows; ++i)
  {
    for (int j = 0; j < NumCols; ++j)
    {
      if (u->m_Map->DistancsMap[j][i]<0&&d->m_Map->DistancsMap[j][i]<0)
      {
        ++temp;
      }
    }
  }
  return temp*1.0;
}
