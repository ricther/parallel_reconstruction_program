#pragma once
//athour xu liu 1/1/2014
//the class of CMedialMap derive from CMap

#include "CMap.h"

class CMedialMap:public CMap
{
  public:
  //map is the medial map of the layer
  CMedialMap(CContour* temp_contour,double** map):CMap(temp_contour)
  {
    medial_axis=map;
  }
  virtual void fillupoutside();
  double** medial_axis;
};
