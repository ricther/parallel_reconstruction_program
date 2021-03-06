#pragma once
#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "CShape.h"
class CDataReader
{
public:
  CDataReader(const char* str):file(str)
  {
    doc.parse<0>(file.data());
    initial_node=doc.first_node("initial");
    data_source_node=doc.first_node("datasource");
    debug_node=doc.first_node("debug");
  }
  rapidxml::xml_node<>* initial_node, *data_source_node, *debug_node;

  int get_value_int(const char* str, int mode =1)
  {
    rapidxml::xml_node<>* node;
    if (mode==1)
    {
      node=initial_node;
    }
    else if (mode==2)
    {
      node=debug_node;
    }
    char* t=node->first_node(str)->value();
    int i=atoi(t);
    return i;
  }

  float get_value_float(const char* str, int mode=1)
  {
    rapidxml::xml_node<>* node;
    if (mode==1)
    {
      node=initial_node;
    }
    else if (mode==2)
    {
      node=debug_node;
    }

    char* t=node->first_node(str)->value();
    float f=atof(t);
    return f;
  }

  bool get_value_bool(const char* str, int mode=1)
  {
    rapidxml::xml_node<>* node;
    if (mode==1)
    {
      node=initial_node;
    }
    else if (mode==2)
    {
      node=debug_node;
    }

    char* t=node->first_node(str)->value();
    if (strcmp(t,"true")==0||strcmp(t,"TRUE")==0||strcmp(t,"1")==0)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  char* get_data_source()
  {

    char* t=data_source_node->first_node("file_name")->value();
    return t;
  }

  void read_file(CShape* m_source)
  {
    m_source->initial(get_data_source());
  }


private:
  rapidxml::file<char> file;
  rapidxml::xml_document<> doc;

}
;
