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
  RangeIterator() = default;
  RangeIterator(int start, int end, RangeIterator *join = nullptr);
  RangeIterator(int start, int end, const RangeIterator &join);
  RangeIterator(const RangeIterator &r);
  ~RangeIterator();

  void first();
  void next();
  bool isDone() const;

  int operator*() { return ind; }
  RangeIterator &operator++() {
    next();
    return *this;
  }
  RangeIterator operator++(int) {
    RangeIterator original(*this);
    next();
    return original;
  }

  bool operator==(const RangeIterator &x) const { return cur == x.cur; }
  bool operator!=(const RangeIterator &x) const { return cur != x.cur; }

  bool intersects(const RangeIterator &other, bool all = true) const;
  bool intersects(int ind, bool all = true) const;

  RangeIterator &operator=(const RangeIterator &r);
  RangeIterator &operator+=(const RangeIterator &r); // Add ranges
  RangeIterator &operator-=(const RangeIterator &r); // Remove ranges

  static RangeIterator end() { return RangeIterator(1, 0); }

  int ind; // The index

  int min() const { return is; }
  int max() const { return ie; }
  RangeIterator *nextRange() const { return n; };

private:
  int is{1}, ie{0};
  RangeIterator* n{nullptr};   // Next range. Doesn't change after creation
  RangeIterator* cur{nullptr}; // Currently iterating. Changes during iteration
  int curend;               // End of current range
  bool delete_next = false; // Flag to delete this->n if we created it
};

#endif // __RANGE_H__
