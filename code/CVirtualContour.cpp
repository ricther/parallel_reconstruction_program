
#include "CVirtualContour.h"
#include "assert.h"
#include "CFileDebug.h"
#include <stdlib.h>
#include "CPoint.h"
#include "CMedialMap.h"

void CVirtualContour::calculate_medial_axis(float layerID)
{
  if (m_contour->map_neighbour.size()==0)
  {
    return;
  }
  std::multimap<float,CContour*>::iterator itr,etr,nitr;
  
  std::pair<std::multimap<float,CContour*>::iterator,std::multimap<float,CContour*>::iterator>pair_itr = m_contour->map_neighbour.equal_range(layerID);

  medial_axis_count++;
  medial_axis=make_2D_float_array(NumRows,NumCols);
  for (int i = 0; i < NumRows; ++i)
  {
    for (int j = 0; j < NumCols; ++j)
    {
      medial_axis[i][j]=255;
    }
  }

  itr=pair_itr.first;etr=pair_itr.second;

  /// there is no key in the multimap
  if (itr==etr)
  {
    return;
  }
  for (; itr!=etr; ++itr)
  {
   
    //    char intStr[10];
    //    sprintf(intStr,"%d",a);
    //    string str = string(intStr);
    
    //    CFileDebug file(str);
    //    file.Output(itr->second->m_Map->DistancsMap);
    //    a++;

    /// there is only one key in the map;
    nitr=itr;nitr++;
    if (nitr==etr)
    {
      break;
    }
    for (; nitr!=etr; ++nitr)
    {
      for (int i = 0; i < NumRows; ++i)
      {
        for (int j = 0; j < NumCols; ++j)
        {
          float first=(int) itr->second->m_Map->get_distance_map(i,j);
          float second=(int) nitr->second->m_Map->get_distance_map(i,j);
          float margin=1;
          if(first-second>=-margin&&first-second<=margin)
          {
            medial_axis[i][j]=0;
          }
          else if (medial_axis[i][j]==0)
          {
            if (itr->second->m_Map->get_distance_map(i,j)<=distance_medialaxis_contour||nitr->second->m_Map->get_distance_map(i,j)<=distance_medialaxis_contour)
            {
              medial_axis[i][j]=255;              
            }
          }
        }
      }
    }
  }
  for (int i = 0; i < NumRows; ++i)
  {
    for (int j = 0; j < NumCols; ++j)
    {
      if( medial_axis[i][j]==0)
      {
        CPoint * new_point=new CPoint();
        //because the x,y in the distance map are exchanged.
        new_point->x=j;new_point->y=i;new_point->z=LayerID*layer_interval_scale;
        new_point->index=new_point->get_index();
        vec_medial_points.push_back(new_point);
      }
    }
  }
}

#include"CFileDebug.h"
void CVirtualContour::calculate_medial_map(float** medial_axis)
{
  m_medial_map = new CMedialMap(this->m_contour,medial_axis);
  m_medial_map->setup();
  m_medial_map->gradient();
  //to-do to test clear the medial_axis may have error,need more test 2014/10/28
  free_2D_float_array(medial_axis);
  //CFileDebug m_debugfile6("./medial_DistancsMap_20");
  //m_debugfile6.Output(m_medial_map->DistancsMap);
  //CFileDebug file("./only_medial_axis");
  //file.Output(medial_axis);
  // CFileDebug m_debugfile7("./medial_sign_20");
  //  m_debugfile7.Output_sign(m_medial_map->SignMap);
  


}

void CVirtualContour::swap_map_medialmap()
{
  assert(m_Map);
  assert(m_medial_map);
  m_temp_map=m_Map;
  m_Map=(CMap*)m_medial_map;
}

void CVirtualContour::swap_medialmap_map()
{
  assert(m_Map);
  assert(m_medial_map);
  m_medial_map=(CMedialMap*)m_Map;
  m_Map=m_temp_map;
  delete m_medial_map;
}


