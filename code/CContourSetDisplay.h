#pragma once
#include "CContourSetDisplay.h"
#include <vector>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <map>
#include "CLayer.h"

template <class Contour>
class CContourSetDisplay
{
public:
  CContourSetDisplay(std::map<float,CLayer*>& temp_map):map_Layer(temp_map)
  {
     float n=999999;
     max_z=-n;min_z=n;
     float o=map_Layer.begin()->first;
     float m=map_Layer.rbegin()->first;
     margin=(m-o)/(map_Layer.size()-1);
     margin=margin+margin*0.1;
  }
  float margin;
private:

   float max_z,min_z;

   void check_level(float level)
   {
     if (level>max_z)
     {
       max_z=level;
     }
     if (level<min_z)
     {
       min_z=level;
     }
   }

   bool check_display(float id,float level)
   {
      float up=id+margin;
      float down=id;
      if (up>=max_z)
      {
        up=max_z;
      }
      if(down<min_z)
      {
        down=min_z;
      }
      /* if (up<down) */
      /* { */
      /*   up=down; */
      /* } */
      if (level<=up&&level>=down)
      {
        return true;
      }
      else
      {
        return false;
      }
   }

public:
  std::map<float,CLayer*>& map_Layer;
  std::vector<Contour*> vec_contour;
  void initialActors()
  {
    std::map<float,CLayer*>::iterator itr,etr;
    itr = map_Layer.begin();etr= map_Layer.end();
    // layer loop
    int count=-1;
    for (;itr!=etr;++itr)
    {
      count++;
      CLayer* temp = itr->second;
      check_level(temp->LayerID);
      //for paper input graph
      /* if (count%8!=0) */
      /* { */
      /*   continue; */
      /* } */
      std::map<int,int> map_contourID=temp->map_contourID;

      std::map<int,int>::iterator itr_c,etr_c;
      itr_c=map_contourID.begin();etr_c=map_contourID.end();
      //contour loop
      for (; itr_c!=etr_c; ++itr_c)
      {
        Contour* new_contour= new Contour(map_Layer);
        new_contour->initialActor(itr->first,itr_c->first);
        vec_contour.push_back(new_contour);
      
      }
    }
  }
  
  void update_actor(float layerID)
  {
    typename std::vector<Contour*>::iterator itr, etr;

    itr=vec_contour.begin(),etr=vec_contour.end();

    double color=1,scale=0;
    for (; itr!=etr; ++itr)
    {
     
      bool b=check_display(layerID,(*itr)->current_LayerID);
      if (b==true)
      {
           (*itr)->m_actor->SetVisibility(true);
           scale=((*itr)->current_LayerID-layerID)/(margin);
           if (scale == 0)
           {
             (*itr)->m_actor->GetProperty()->SetColor(1,0,0);
           }
           else
           {
             (*itr)->m_actor->GetProperty()->SetColor(0.2,0.6,0);
             // (*itr)->m_actor->GetProperty()->SetColor(0.8*scale,0.4*scale,0);
           }
           (*itr)->m_actor->GetProperty()->SetLineWidth(3);
           //                        (*itr)->m_actor->GetProperty()->SetColor(color,color-scale*color,color-scale*color);        

      }
      else
      {
          (*itr)->m_actor->SetVisibility(false);
          (*itr)->m_actor->GetProperty()->SetColor(color,color,color);
      }
      (*itr)->m_actor->Modified();
 
    }
  }

  void show_all_actor()
  {
    typename std::vector<Contour*>::iterator itr, etr;

    itr=vec_contour.begin(),etr=vec_contour.end();
    for (; itr!=etr; ++itr)
    {
      (*itr)->m_actor->SetVisibility(true);
      (*itr)->m_actor->Modified();
    }
  }

};

