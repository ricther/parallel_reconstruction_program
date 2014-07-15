#include "CFileDebug.h"
#include "CPoint.h"
bool ifprint=true;
using namespace std;
CFileDebug::CFileDebug(string FileName)
{
  int num =-1;
  filename=FileName;
  num=filename.find("/");
  if(num!=-1)
  {
    filename.erase(0,num+1);
  }
  filename="../Debug/"+filename;
  if (is_noprint())
  {
    return;
  }
  m_file.open(filename.c_str(),ios_base::out);
  bool test=m_file.is_open();
}

void CFileDebug::Output(double **& map)
{
  if(is_noprint())
  {
    return;
  }
  int row=NumRows;
  int col=NumCols;
  int irow=0,icol=0;
  int temp=0;
  for(;irow<row;++irow)
  {
    for(icol=0;icol<col;++icol)
    {
      if(map[irow][icol]>10)
      {
        temp=1;
      }
      else if (map[irow][icol]<0)
      {
        temp=3;
      }
      else
      {
        temp=map[irow][icol];
      }

      m_file<<temp;
    }
    m_file<<std::endl;
  }
  m_file.close();

}

void CFileDebug::Output_sign(unsigned char **& map)
{
  if (is_noprint())
  {
    return;
  }
  int row=NumRows;
  int col=NumCols;
  int irow=0,icol=0;
  int temp=0;
  for(;irow<row;++irow)
  {
    for(icol=0;icol<col;++icol)
    {
      if(map[irow][icol]<10)
      {
        temp=0;
      }
      else
      {
        temp=1;
      }

      m_file<<temp;
    }
    m_file<<std::endl;
  }
  m_file.close();  
}

void CFileDebug::Output_points(std::vector<CPoint*> vec)
{
  if (is_noprint())
  {
    return;
  }
  std::vector<CPoint*>::iterator itr,etr;
  itr=vec.begin();etr=vec.end();
  for (;itr!=etr; ++itr)
  {
    m_file<<(int)(*itr)->x;
    m_file<<"\t";
    m_file<<(int)(*itr)->y<<"\r\n";
  }
  m_file.close();
  
}

void CFileDebug::Output_distdata(int i, int j,double value)
{
  if (is_noprint())
  {
    return;
  }
  m_file<<i<<"\t"<<j<<"\t"<<value<<"\r\n";
}

void CFileDebug::close()
{
  m_file.close();
}

bool CFileDebug::is_noprint()
{
  if (ifprint==false)
  {
    return true;
  }
  else
  {
    return false;
  }
}
void CFileDebug::Output_intensity_data(vector<CPoint*>vec_point,vector<double> vec)
{
  std::vector<double>::iterator itr,etr;
  std::vector<CPoint*>::iterator pitr,petr;
  pitr=vec_point.begin();petr=vec_point.end();
  itr = vec.begin();etr= vec.end();
  for(;itr!=etr&&(pitr!=petr);++itr,++pitr)
  {
    m_file<<(int)(*pitr)->x<<"\t"<<(int)(*pitr)->y<<"\t"<<*itr<<"\r\n";
  }
  m_file.close();
}
