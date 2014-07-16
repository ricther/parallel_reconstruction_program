#include "CMap.h"
#include "CLayer.h"
#include <stdlib.h>
#include <string.h>
#include "CPoint.h"
#include"Regionfill.h"
#include "CFileDebug.h"
#include "CContour.h"
#include <iostream>
using namespace std;
CMap::CMap(CContour* temp_contour):NUMROWS(NumRows),NUMCOLS(NumCols)
{
  m_contour=temp_contour;
  SignMap=NULL;
  DistancsMap=NULL;
  area=0;
}

CMap::~CMap()
{
  if(SignMap!=NULL)
  freeMap((void**)SignMap);
  
  if(DistancsMap!=NULL)
  freeMap((void**)DistancsMap);

  freeMap((void**)gx);
  freeMap((void**)gy);
}

void CMap::setup()
{
  string filestr=m_contour->filename;
  //string tempstr=filestr+"_before_region_fill";
  // CFileDebug m_debugfile(tempstr.c_str());
  //tempstr=filestr+"_after_region_fill";
  //CFileDebug m_debugfile2(tempstr.c_str());
  //tempstr=filestr+"_sign";

  mallocMap_char(SignMap);
  mallocMap_double(DistancsMap);
  interp_closed();
  //tempstr=filestr+"_interpolate_points";
  //CFileDebug m_debugfile4(tempstr.c_str());
  //m_debugfile4.Output_points(m_contour->vec_Points_Inter);
  drawcontour();
  // if(m_contour->LayerID==1&&m_contour->contourID==1)
  // {
  //   CFileDebug m_debugfile3("synthetic2_sign_1_1");
  //   m_debugfile3.Output_sign(SignMap);
  // }
  //m_debugfile.Output(DistancsMap);
  fillupoutside();
  // string tempstr = filestr+"_after_fillupoutside";
  // CFileDebug m_debugfile5(tempstr.c_str());/
  // m_debugfile5.Output_sign(SignMap);
  // tempstr= filestr+"_after_cal_dist";
  // CFileDebug m_debugfile6(tempstr.c_str());
  caldistancesmap();
  //  m_debugfile6.Output(DistancsMap);
  // m_debugfile2.Output(DistancsMap);
}

bool CMap::mallocMap_double(float **& array)
{
  array=0;
  float *arr_ptr;
  arr_ptr=(float*) calloc(NUMCOLS*NUMROWS,sizeof(float));//calloc can initialize the memory with 0
  //  memset(arr_ptr,NUMCOLS*NUMROWS,sizeof(double)*NUMCOLS*NUMROWS);
  //  std::cout<<"mallocMap_double_size:"<<NUMCOLS*NUMROWS*sizeof(double)<<"double_size:"<<sizeof(double)<<"\n";
  if (arr_ptr==NULL)
  {
    return false;
  }
  array =(float**) calloc(NUMROWS,sizeof(float*));
  if (array==NULL)
  {
    return false;
  }
  for (int i = 0; i < NUMROWS; ++i)
  {
     array[i]=&arr_ptr[i*NUMCOLS];
     //    array[i]=(double*)malloc(sizeof(double)*NUMCOLS);
  }
  return true;
}

bool CMap::mallocMap_char(unsigned char **& array)
{
  array=0;
  unsigned char *arr_ptr;
  arr_ptr=(unsigned char*) calloc(NUMCOLS*NUMROWS,sizeof(unsigned char));//calloc can initialize the memory with 0
  //  memset(arr_ptr,NUMCOLS*NUMROWS,sizeof(unsigned char)*NUMCOLS*NUMROWS);
  if (arr_ptr==NULL)
  {
    return false;
  }
  array =(unsigned char**) calloc(NUMROWS,sizeof(unsigned char*));
  if (array==NULL)
  {
    return false;
  }
  for (int i = 0; i < NUMROWS; ++i)
  {
     array[i]=&arr_ptr[i*NUMCOLS];
     //    array[i]=(unsigned char*)malloc(sizeof(unsigned char)*NUMCOLS);
  }
  return true;
}

void CMap::freeMap(void** array)
{
  free((char*)array[0]);
  free((char*)array);
}

void CMap::interp_closed()
{
  if (m_contour->PointNum>0)
  {
    int size=m_contour->PointNum;
    int newrow[2000], newcol[2000];
    int temprow[200],tempcol[200];
    int numpoints = 0;
    int newpoints;
    int i=0;
    for ( i= 0; i < size-1; i++)
    {	
      newpoints = 0;
      numpoints++;
      newpoints = digi_img_line(m_contour->vec_Points_Origin[i],m_contour->vec_Points_Origin[i+1],temprow,tempcol);
      insertPoint( m_contour->vec_Points_Origin[i]);
      insertNewPoints(newpoints,temprow,tempcol);
      numpoints = numpoints+newpoints;
    }
               
    insertPoint(m_contour->vec_Points_Origin.back()); 
    newpoints = 0;
    newpoints = digi_img_line(m_contour->vec_Points_Origin.back(),m_contour->vec_Points_Origin.front(),temprow,tempcol);
    numpoints++;
    insertNewPoints(newpoints,temprow,tempcol);
    numpoints = numpoints+newpoints;

    insertPoint(m_contour->vec_Points_Origin.front());
  
    numpoints++;        
  }
  if (m_contour->LayerID==1&&m_contour->contourID==1)
  {
    CFileDebug m_debugfile3("synthetic2_points_1_1");
    m_debugfile3.Output_points(m_contour->vec_points);
  }
}
    
/***********************/
// Digital line drawning sample points by double to int, and interploate point in void pixelint
int CMap::digi_img_line(CPoint* point1,CPoint* point2,int *xint, int *yint)
{
  int x1=point1->x , y1=point1->y,x2=point2->x,y2=point2->y;
  int d,x,y,ax,ay,sx,sy,dx,dy;
  int next=0;
  dx=x2-x1;  
  ax=abs(dx)<<1;  //why multiply 2 question..
  if (dx>=0) sx=1; else sx=-1;
  dy=y2-y1;  
  ay=abs(dy)<<1;   //why multiply 2 question..
  if (dy>=0) sy=1; else sy=-1;
  x=x1;
  y=y1;
  
  if (ax>ay) 
  {
    d=ay-(ax>>1); //why ay minus half ax question..
    for (;;) 
    {
      xint[next] = x;
      yint[next] = y;
      next=next+1;
      if (x==x2) return(next-1);
      if (d>=0) {y+=sy; d-=ax;}
      x+=sx;
      d+=ay;
    }
  }
  else 
  {			
    d=ax-(ay>>1);
    for (;;) 
    {
      xint[next] = x;
      yint[next] = y;
      next=next+1;
      if (y==y2) return(next-1);
      if (d>=0) {x+=sx; d-=ay;}
      y+=sy;
      d+=ax;
    }
  }
}

void CMap::insertNewPoints(int Num, int* temprow, int*tempcol)
{
  for (int j = 0; j < Num; j++)
  {
    CPoint* temp = new CPoint(temprow[j],tempcol[j],m_contour->LayerID);
    m_contour->vec_Points_Inter.push_back(temp);
  }
}

void CMap::insertPoint(CPoint* temp)
{
  m_contour->vec_Points_Inter.push_back(temp);
}

///////////////////////////////////////////////
///the edge pixel in the signmap will be 255
///the edge pixel in the distances map will nbe 0
void CMap::drawcontour()
{
  memset(DistancsMap[0],'w',sizeof(float)*NUMCOLS*NUMROWS);// initialize the DistancsMap with a big number
  memset(SignMap[0],0,sizeof(unsigned char)*NUMCOLS*NUMROWS);// initialize the SignMap with 0;
  for(int i=0;i<m_contour->vec_Points_Inter.size();++i)
  {
    int x=m_contour->vec_Points_Inter[i]->x;
    int y=m_contour->vec_Points_Inter[i]->y;
    SignMap[x][y]=255;
    DistancsMap[x][y]=0;
  }
}
//the outside pixel value is 128, the edge value is 128 and the pixel value inside the contour is 0.
void CMap::fillupoutside()
{
  //regionfill can fill the pixel which is outside the contour but with value 255
    Regionfill::regionfill(SignMap,NUMROWS,NUMCOLS,1,1,255);
    //relabel the contour and the outside with value 128
    for (int i = 0; i < NUMROWS; ++i)
    {
      for (int j = 0; j < NUMCOLS; ++j)
      {
        if(SignMap[i][j]==255)
          SignMap[i][j]=128;
        else
          SignMap[i][j]=0;
      }
    }
    //relabel the contour with value 255
    for (int i = 0; i < m_contour->vec_Points_Inter.size(); ++i)
    {
      SignMap[(int)m_contour->vec_Points_Inter[i]->x][(int)m_contour->vec_Points_Inter[i]->y]=255;
    }
}

void CMap::caldistancesmap()
{
  double M_SQRT_2 = 1.4142135623730950488016887220097;
  int ONE=1;
  int k=0,l=0;
  float t;

  for(int i=0;i<NUMROWS;i++)
  {
    for(int j=0;j<NUMCOLS;j++)
    {
      t=DistancsMap[i][j];

      k=i-1;l=j-1;
      if ((k>=0)&&(l>=0)&&(t>DistancsMap[k][l]+M_SQRT_2)) 
        t=DistancsMap[k][l]+M_SQRT_2;
      k=i-1;l=j;
      if ((k>=0)&&(t>DistancsMap[k][l]+ONE)) 
        t=DistancsMap[k][l]+ONE;
      k=i-1;l=j+1;
      if ((k>=0)&&(l<NUMCOLS)&&(t>DistancsMap[k][l]+M_SQRT_2))
        t=DistancsMap[k][l]+M_SQRT_2;
      k=i;l=j-1;
      if ((l>=0)&&(t>DistancsMap[k][l]+ONE)) 
        t=DistancsMap[k][l]+ONE;
      DistancsMap[i][j]=t;
    }
  }
  
   for(int i=NUMROWS-1;i>=0;i--)
  {
    for(int j=NUMCOLS-1;j>=0;j--)
    {
      t=DistancsMap[i][j];
      float test=t;
      k=i+1;l=j-1;if ((k<NUMROWS)&&(l>=0)&&(t>DistancsMap[k][l]+M_SQRT_2)) t=DistancsMap[k][l]+M_SQRT_2;
      k=i+1;l=j;if ((k<NUMROWS)&&(t>DistancsMap[k][l]+ONE)) t=DistancsMap[k][l]+ONE;
      k=i+1;l=j+1;if ((k<NUMROWS)&&(l<NUMCOLS)&&(t>DistancsMap[k][l]+M_SQRT_2)) t=DistancsMap[k][l]+M_SQRT_2;
      k=i;l=j+1;if ((l<NUMCOLS)&&(t>DistancsMap[k][l]+ONE)) t=DistancsMap[k][l]+ONE;
      DistancsMap[i][j]=t;
    }
  }

  for (int i = 0; i < NUMROWS-1; ++i)
  {
    for (int j = 0; j < NUMCOLS-1; ++j)
    {
      if (SignMap[i][j]==0)
      {
        DistancsMap[i][j]=-DistancsMap[i][j];
        ++area;
      }
    }
  }
  float Temp[NUMROWS][NUMCOLS];
  for (int i = 0; i < NUMROWS; ++i)
  {
    for (int j = 0; j < NUMCOLS; ++j)
    {
      Temp[i][j]=DistancsMap[i][j];
    }
  }
  for (int i = 0; i < NUMROWS; ++i)
  {
    for (int j = 0; j < NUMCOLS; ++j)
    {
      DistancsMap[j][i]=Temp[i][j];
    }
  }
  //////////////////////////////////////////////////////
  // string file_name;
  // file_name=m_contour->filename+"_dist_output";
  // CFileDebug m_file(file_name.c_str());
  // for (int i = 0; i < NUMROWS; ++i)
  // {
  //   for (int j = 0; j < NUMCOLS; ++j)
  //   {
  //     m_file.Output_distdata(i,j,DistancsMap[i][j]);
  //   }
  // }
  ////////////////////////////////////////////////////
}

void CMap::gradient()
{
  mallocMap_double(gx);
  mallocMap_double(gy);
  int r1 = NUMROWS -1;
  int c1 = NUMCOLS -1;

  // string filestr=m_contour->filename;
  // filestr+="_gradient_gx_new";
  // CFileDebug file_gx(filestr);
  //filestr=m_contour->filename+"_gradient_gy_new";
  // CFileDebug file_gy(filestr);
  for (int i = 1; i < r1; ++i)
  {
    for (int j = 1; j < c1; ++j)
    {
      gx[i][j] = (DistancsMap[i][j+1] - DistancsMap[i][j-1]) / 2.0;
      gy[i][j] = (DistancsMap[i+1][j] - DistancsMap[i-1][j]) / 2.0;
      //file_gx.Output_distdata(i,j,gx[i][j]);
      //file_gy.Output_distdata(i,j,gy[i][j]);
    }
  }
  //borders
  for (int i = 1; i < r1; ++i)
  {
    //left
    gx[i][0] = DistancsMap[i][1]-DistancsMap[i][0];
    gy[i][0] = (DistancsMap[i+1][0] - DistancsMap[i-1][0]) /2.0;
    //file_gx.Output_distdata(i,0,gx[i][0]);
    //file_gy.Output_distdata(i,0,gy[i][0]);
    //right
    gx[i][r1] = DistancsMap[i][r1] - DistancsMap[i][r1-1];
    gy[i][r1] = (DistancsMap[i+1][r1] -DistancsMap[i-1][r1]) /2.0;
    // file_gx.Output_distdata(i,r1,gx[i][r1]);
    //file_gy.Output_distdata(i,r1,gy[i][r1]);
  }

  
  for (int j = 1; j < r1; j++)
  {
    // Top 
    gx[0][j] = (DistancsMap[0][j+1]-DistancsMap[0][j-1]) / 2.0;
    gy[0][j] = DistancsMap[1][j]-DistancsMap[0][j];

    //file_gx.Output_distdata(0,j,gx[0][j]);
    // file_gy.Output_distdata(0,j,gy[0][j]);
    // Bottom
    gx[c1][j] = (DistancsMap[c1][j+1]-DistancsMap[c1][j-1]) / 2.0;
    gy[c1][j] = DistancsMap[c1][j]-DistancsMap[c1-1][j];

    //file_gx.Output_distdata(c1,j,gx[c1][j]);
    //file_gy.Output_distdata(c1,j,gy[c1][j]);
  }
  //  file_gx.close();
  // file_gy.close();
  // corners
  gx[0][0] = DistancsMap[0][1]- DistancsMap[0][0];
  gy[0][0] = DistancsMap[1][0]- DistancsMap[0][0];

  gx[0][c1] = DistancsMap[0][c1] - DistancsMap[0][c1 -1];
  gy[0][c1] = DistancsMap[1][c1] - DistancsMap[0][c1];

  gx[r1][0] = DistancsMap[r1][1] - DistancsMap[r1][0];
  gy[r1][0] = DistancsMap[r1][0] - DistancsMap[r1-1][0];
  
  gx[r1][c1] = DistancsMap[r1][c1] - DistancsMap[r1][c1 -1];
  gy[r1][c1] = DistancsMap[r1][c1] - DistancsMap[r1-1][c1];
  //  freeMap((void**)SignMap);
  //SignMap=NULL;
}
