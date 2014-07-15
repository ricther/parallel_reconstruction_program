#include "CCorespondDisplay.h"
#include "CShape.h"
void CCorespondDisplay:: initial_actors(CShape * m_shape)
{

  int size= m_shape->vec_Cor.size();
  for (int i = 0; i < size; ++i)
  {
    CLineSetDisplay* new_line_set=new CLineSetDisplay();
    new_line_set->initial_actors(m_shape->vec_Cor[i]); 
    map_line_actors.insert(std::make_pair(m_shape->vec_Cor[i]->CorrespondLayer1,new_line_set));
    check_level(m_shape->vec_Cor[i]->CorrespondLayer1);
    check_level(m_shape->vec_Cor[i]->CorrespondLayer2);
  }
}

void CCorespondDisplay:: initial_actor(CShape * m_shape)
{

  int size= m_shape->vec_Cor.size();
  for (int i = 0; i < size; ++i)
  {
    CLineSetDisplay* new_line_set=new CLineSetDisplay();
    new_line_set->initial_actor(m_shape->vec_Cor[i]); 
    map_line_actors.insert(std::make_pair(m_shape->vec_Cor[i]->CorrespondLayer1,new_line_set));
    check_level(m_shape->vec_Cor[i]->CorrespondLayer1);
    check_level(m_shape->vec_Cor[i]->CorrespondLayer2);
  }
}

void CCorespondDisplay::update_actors(float layerID,int& lineID)
{
  std::map<float,CLineSetDisplay*>::iterator itr,etr;
  itr=map_line_actors.begin();etr=map_line_actors.end();
  for (;itr!=etr;++itr)
  {
    if (draw_by_layer)
    {
      bool b=check_display(layerID,itr->first);
      if (b==true)
      {
        itr->second->m_actor->SetVisibility(true);

      }
      else
      {
        itr->second->m_actor->SetVisibility(false);

      }
      itr->second->m_actor->Modified();
    }
    else
    {
      bool b=check_display(layerID,itr->first);
      if (b==true)
      {
        itr->second->set_actor_visibility(true,lineID);
      }
      else
      {
        itr->second->set_actor_visibility(false,lineID);
      }
    }
  }
}








