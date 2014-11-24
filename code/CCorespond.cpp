#include "CCorespond.h"
#include "CTriangle.h"
#include "iostream"
void CCorrespond::produce_triangle()
{
  std::map<int,Parapoint*>::iterator itr=map_cor.begin(),etr=map_cor.end();
  int count=0;
  std::map<int,int> map_int_index;
  for (; itr!=etr; ++itr)
  {
    map_int_index.insert(std::make_pair(count++,itr->second->point1.index));
    map_int_index.insert(std::make_pair(count++,itr->second->point2.index));
  }
  for (int i = 0; i < count; ++i)
  {
    //    if (i%2==0)
    {
          CTriangle* tri= new CTriangle();
          tri->a=map_int_index[i];
          tri->b=map_int_index[(i+1)%count];
          tri->c=map_int_index[(i+2)%count];
          vec_Triangle.push_back(tri);
    //    std::cout<<"tri:"<<tri->a<<"\t"<<tri->b<<"\t"<<tri->c<<"\n";
          tri = new CTriangle();
          tri->a=map_int_index[(i+2)%count];
          tri->b=map_int_index[(i+1)%count];
          tri->c=map_int_index[(i+3)%count];
          vec_Triangle.push_back(tri);
    //    std::cout<<"tri:"<<tri->a<<"\t"<<tri->b<<"\t"<<tri->c<<"\n";

    }
    // else{
    //       CTriangle* tri= new CTriangle();
    //       tri->a=map_int_index[i];
    //       tri->b=map_int_index[(i+1)%count];
    //       tri->c=map_int_index[(i+3)%count];
    //       vec_Triangle.push_back(tri);
    // //    std::cout<<"tri:"<<tri->a<<"\t"<<tri->b<<"\t"<<tri->c<<"\n";
    //       tri = new CTriangle();
    //       tri->a=map_int_index[(i+0)%count];
    //       tri->b=map_int_index[(i+3)%count];
    //       tri->c=map_int_index[(i+2)%count];
    //       vec_Triangle.push_back(tri);
    // //    std::cout<<"tri:"<<tri->a<<"\t"<<tri->b<<"\t"<<tri->c<<"\n";

    // }
      
  }
  
}
