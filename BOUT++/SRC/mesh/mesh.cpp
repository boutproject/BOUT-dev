
#include "mesh.h"

static Mesh::wtime_comms = 0.0;

/**************************************************************************
 * Data sources
 **************************************************************************/

int Mesh::add_source(GridDataSource &source)
{
  return add_source(&source);
}

int Mesh::add_source(GridDataSource *source)
{
  source_list.push_front(source);
  return 0;
}

int Mesh::load(const GridDataSource &source)
{
  std::list<GridDataSource*> old_list;
  
  old_list = source_list;
  source_list.clear();
  
  if(add_source(source))
    return 1;
  
  return load();
}

GridDataSource* Mesh::find_source(const char *name)
{
  for(std::list<GridDataSource*>::iterator it = source_list.begin(); 
      it != source_list.end(); it++) {
    
    // Query this source
    if((*it) != NULL)
      if((*it)->hasVar(name)) {
	return *it;
      }
  }
  return NULL;
}

/**************************************************************************
 * Data get routines
 **************************************************************************/

int Mesh::get(Vector2D &var, const char *name)
{
  return get(var, string(name));
}

int Mesh::get(Vector3D &var, const char *name)
{
  return get(var, string(name));
}

int Mesh::get(Vector2D &var, const string &name)
{
#ifdef CHECK
  msg_stack.push("Loading 2D vector: Mesh::get(Vector2D, %s)", name.c_str());
#endif

  if(var.covariant) {
    output << "\tReading covariant vector " << name << endl;
    
    get(var.x, name+"_x");
    get(var.y, name+"_y");
    get(var.z, name+"_z");
    
  }else {
    output << "\tReading contravariant vector " << name << endl;
    
    get(var.x, name+"x");
    get(var.y, name+"y");
    get(var.z, name+"z");
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif  

  return 0;
}

int Mesh::get(Vector3D &var, const string &name)
{
#ifdef CHECK
  msg_stack.push("Loading 3D vector: Mesh::get(Vector3D, %s)", name.c_str());
#endif

  if(var.covariant) {
    output << "\tReading covariant vector " << name << endl;
    
    get(var.x, name+"_x");
    get(var.y, name+"_y");
    get(var.z, name+"_z");
    
  }else {
    output << "\tReading contravariant vector " << name << endl;
    
    get(var.x, name+"x");
    get(var.y, name+"y");
    get(var.z, name+"z");
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif  

  return 0;
}

/**************************************************************************
 * Communications
 **************************************************************************/

int Mesh::communicate(FieldData &f)
{
  FieldGroup group;
  group.add(f);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2)
{
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2, FieldData &f3)
{
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  group.add(f3);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4)
{
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  group.add(f3);
  group.add(f4);
  return communicate(group);
}

comm_handle Mesh::send(FieldData &f)
{
  FieldGroup group;
  group.add(f);
  return send(group);
}
