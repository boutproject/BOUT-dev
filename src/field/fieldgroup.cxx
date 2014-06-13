
#include <bout/fieldgroup.hxx>

FieldGroup::FieldGroup(FieldData &f) {
  fvec.push_back(&f);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2) {
  fvec.push_back(&f1); fvec.push_back(&f2);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3) {
  fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4) {
  fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3); fvec.push_back(&f4);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4, FieldData &f5) {
  fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3); fvec.push_back(&f4); fvec.push_back(&f5);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4, FieldData &f5, FieldData &f6) {
  fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3); fvec.push_back(&f4);
  fvec.push_back(&f5); fvec.push_back(&f6);
}

FieldGroup& FieldGroup::operator+(const FieldGroup &other){
  vector<FieldData*> vec;
  vec=other.get();
  for(int i=0;i<vec.size();i++){
    this->add(*(vec[i]));
  };
  return *this;
};

FieldGroup& FieldGroup::operator+=(const FieldGroup &other){
  return (*this)+other;
};
