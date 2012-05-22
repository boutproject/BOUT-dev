


#include <bout/sys/range.hxx>

RangeIterator::RangeIterator(int start, int end, RangeIterator* join) 
  : is(start), ie(end), n(join) {
  
  cur = this;
  if(start > end) {
    // Null range
    cur = n;
  }
  
  if(cur != 0) {
    ind = cur->is;
    curend = cur->ie;
  }
}

RangeIterator::RangeIterator(int start, int end, const RangeIterator& join)
  : is(start), ie(end) {
  
  cur = this;

  n = new RangeIterator(join);
  
  if(start > end) {
    // Null range
    cur = n;
  }
  
  if(cur != 0) {
    ind = cur->is;
    curend = cur->ie;
  }
}

void RangeIterator::first() {
  cur = this;
  ind = is;
  curend = ie;
  
  if(is > ie) {
    // Null range, skip to next
    cur = cur->n;
    if(cur != 0) {
      ind = cur->is;
      curend = cur->ie;
    }
  }
}

void RangeIterator::next() {
  ind++;
  if(ind > curend) {
    // End of this range
    
    cur = cur->n;
    if(cur != 0) {
      // Another range
      ind = cur->is;
      curend = cur->ie;
    }
  }
}

bool RangeIterator::isDone() const {
  return cur == 0;
}
