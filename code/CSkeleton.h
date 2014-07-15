#pragma once
#include "CPoint.h"
#include "CLayer.h"
#include <map>
#include <queue>
class CSkeletalPoint
{
  public:
  CSkeletalPoint(int);
  //  CPoint* pre_node;
  //  CPoint* next_node;
  CPoint* node;
  int ID;
  double color;
  int use_as_lower_skeletal_counter;
  int use_as_higher_skeletal_counter;
};

class CBranch
{
  public:
  CBranch(int,std::map<int,CSkeletalPoint*>&,std::map<float,CLayer*>::iterator,std::map<float,CLayer*>::iterator);
  int nodeID_counter;
  int branchID;
  std::map<int,CSkeletalPoint*> map_id_spoint;
  void build_branch(int sign,CPoint*, std::map<float,CLayer*>::iterator,std::map<float,CLayer*>::iterator etr);
  CSkeletalPoint* check_skeletalpoint(CPoint*);
  CPoint* find_next_skeletalpoint(int ,bool ,CPoint* , CLayer*);
  CSkeletalPoint* new_skeletalpoint(CPoint* point);
  private:
  std::map<int,CSkeletalPoint*>::iterator itr,etr;
  std::map<int,CSkeletalPoint*>& map_id_skeletal;
  std::map<float,CLayer*>::iterator map_begin,map_end;
}
;

struct skeletal_branch
{
  std::map<float,CLayer*>::iterator itr,etr;
};

class CSkeleton
{
  public:
  CSkeleton();
  int branchID_counter;
  int skeletal_point_counter;
  std::map<int,CBranch*> map_id_branch;
  void split_branch(int,int,int);
  void build_skeleton(std::map<float,CLayer*>&);
  void build_skeleton_new(std::map<float ,CLayer*>&);
  void build_skeleton_from_queue(std::queue<skeletal_branch>&);
  void initial_skeletal_points(std::map<float,CLayer*>&);
  std::map<int,CSkeletalPoint*> map_id_skeletal;
  std::map<float,CLayer*>::iterator map_begin,map_end;
  private:
  int check_contour_num(int,int);
  void insert_queue(std::map<float,CLayer*>::iterator itr, std::map<float,CLayer*> :: iterator etr, std::queue<skeletal_branch>& queue);
}
;
