
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

FieldGroup FieldGroup::operator+(const FieldGroup &other){
  FieldGroup temp=(*this); //Temporary field group -- Initialise temporary to hold contents of this

  //Now add contents of other
  for(int i=0;i<other.fvec.size();i++){
    temp.add(*(other.fvec[i]));
  };

  //Return copy of temp
  return temp;
};

FieldGroup& FieldGroup::operator+=(const FieldGroup &other){
  (*this)=(*this)+other;
  return *this;
};

void FieldGroup::makeUnique(){
  //Need to sort vector before making unique
  std::sort(fvec.begin(), fvec.end());
  //Remove duplicate entries (doesn't resize vector though)
  //auto last = std::unique(fvec.begin(), fvec.end());
  //Nicer to do the above but can't use auto until we have c++11 by default
  vector<FieldData*>::iterator last = std::unique(fvec.begin(), fvec.end());
  //Resizes vector to remove memory no longer required
  fvec.erase(last, fvec.end());
};
