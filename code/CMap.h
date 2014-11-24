#pragma once
#include <stdlib.h>
#include <assert.h>
// this is the class about the bit map
class CContour;
class CPoint;

extern const int MatrixRes;
extern const int NumRows;
extern const int NumCols;

class CMap
{
public:
  CMap(CContour*);
  ~CMap();
  const int NUMROWS;
  const int NUMCOLS;
  unsigned char** SignMap;

  CContour* m_contour;
  void setup();
  bool mallocMap_float(float**&);
  bool mallocMap_char(unsigned char**&);
  void freeMap(void**);
  void gradient();
  int area;
protected:
  float **gx,**gy;
  float** DistancsMap;//if the image is too large, here can use double
public:
  float**& get_distance_map()
  {
    if(DistancsMap==NULL)
    {
      setup();
    }
    return DistancsMap;
  }
  
  float**& get_gx()
  {
    if (gx==NULL)
    {
      gradient();
    }
    return gx;
  }
  
  float**& get_gy()
  {
    if (gy==NULL)
    {
      gradient();
    }
    return gy;
  }

  float get_distance_map(int i,int j)
  {
    if(DistancsMap==NULL)
    {
      setup();
    }
    return DistancsMap[i][j];
  }
  float get_gx(int i,int j)
  {
    if (gx==NULL)
    {
      if (DistancsMap==NULL)
      {
        assert(false);
      }
      gradient();
    }
    return gx[i][j];
  }
  float get_gy(int i , int j)
  {
    if (gy==NULL)
    {
      gradient();
    }
    return gy[i][j];
  }

protected:
  void interp_closed();
  int  digi_img_line(CPoint*,CPoint*,int*,int*);
  void insertNewPoints(int Num, int*, int*);
  void insertPoint(CPoint*);
  void drawcontour();
  virtual void fillupoutside();
  void caldistancesmap();
};
