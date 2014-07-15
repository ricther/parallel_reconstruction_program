#pragma once
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
  double** DistancsMap;
  CContour* m_contour;
  void setup();
  bool mallocMap_double(double**&);
  bool mallocMap_char(unsigned char**&);
  void freeMap(void**);
  void gradient();
  double**gx,**gy;
  int area;
protected:
  void interp_closed();
  int  digi_img_line(CPoint*,CPoint*,int*,int*);
  void insertNewPoints(int Num, int*, int*);
  void insertPoint(CPoint*);
  void drawcontour();
  virtual void fillupoutside();
  void caldistancesmap();
};
