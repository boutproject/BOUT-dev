/*
  Provides a convenient and relatively efficient way to loop
  over a set of ranges
  
  Examples
  --------
  
  RangeIterator it(1,4);
  for(it.first(); !it.isDone(); it.next())
    cout << *it;
  
  or 

  for(it.first(); it != RangeIterator::end(); it++)
    cout << *it;
  
  Call to first() is optional, and resets to start of ranges

  for(RangeIterator it(1,4); !it.isDone(); it++)
    cout << *it;
  
*/

#ifndef __RANGE_H__
#define __RANGE_H__

class RangeIterator {
 public:
  /// Can be given a single range
  RangeIterator(int start, int end, RangeIterator* join=0);
  RangeIterator(int start, int end, const RangeIterator& join);
  
  void first();
  void next();
  bool isDone() const;
  
  int operator*() {return ind;}
  RangeIterator& operator++() { next(); return *this; }
  RangeIterator& operator++(int) { next(); return *this; }
  
  bool operator==(const RangeIterator &x) const {
    return cur == x.cur;
  }
  bool operator!=(const RangeIterator &x) const {
    return cur != x.cur;
  }

  static RangeIterator end() { return RangeIterator(1,0);}

  int ind; // The index
 private:
  int is, ie;
  RangeIterator *n; // Next range. Doesn't change after creation
  RangeIterator *cur;  // Currently iterating. Changes during iteration
  int curend; // End of current range
};

#endif // __RANGE_H__
