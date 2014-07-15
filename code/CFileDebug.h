#pragma once
#include <fstream>
#include <string>
#include <vector>
class CPoint;
extern const int NumRows;
extern const int NumCols;

class CFileDebug
{
  public:
  CFileDebug(std::string FileName);
  void Output(double**& map);
  void Output_sign(unsigned char **&map);
  void Output_points(std::vector<CPoint*>);
  void Output_distdata(int,int,double);
  void Output_intensity_data(std::vector<CPoint*>,std::vector<double>);
  void close();

  private:
  std::string filename;
  std::ofstream m_file;
  bool is_noprint();
};
