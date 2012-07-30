#ifndef __FIELDGROUP_H__
#define __FIELDGROUP_H__

#include "field_data.hxx"

#include <vector>
using std::vector;

/// Group together fields
class FieldGroup {
 public:
  FieldGroup() {}
  FieldGroup(FieldData &f);
  FieldGroup(FieldData &f1, FieldData &f2);
  FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3);
  FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4);
  FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4, FieldData &f5);
  FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4, FieldData &f5, FieldData &f6);
  FieldGroup(const FieldGroup &other) {fvec = other.fvec;}
  
  FieldGroup& operator=(const FieldGroup &other) {fvec = other.fvec; return *this;}
  
  void add(FieldData &f) {fvec.push_back(&f);}
  void add(FieldData &f1, FieldData &f2) {
    fvec.push_back(&f1); fvec.push_back(&f2);}
  void add(FieldData &f1, FieldData &f2, FieldData &f3) {
    fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3);}
  void add(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4) {
    fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3); fvec.push_back(&f4);}
  void add(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4, FieldData &f5) {
    fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3); 
    fvec.push_back(&f4); fvec.push_back(&f5);}
  void add(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4, FieldData &f5, FieldData &f6) {
    fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3); 
    fvec.push_back(&f4); fvec.push_back(&f5); fvec.push_back(&f6);}
  
  const vector<FieldData*> get() const {return fvec;} 
  void clear() {fvec.clear();}
 private:
  vector<FieldData*> fvec; // Vector of fields
};

#endif // __FIELDGROUP_H__
