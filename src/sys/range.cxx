


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

RangeIterator::RangeIterator(const RangeIterator& r) {
  ind    = r.ind;
  is     = r.is;
  ie     = r.ie;
  n      = r.n;
  cur    = r.cur;
  if(cur == &r)
    cur = this;
  curend = r.curend;
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
  if(isDone())
    return;
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

bool RangeIterator::intersects(const RangeIterator &other, bool all) const {
  if((other.is <= ie) && (other.ie >= is))
    return true;
  if(all && (n != 0))
    return n->intersects(other, all);
  return false;
}

bool RangeIterator::intersects(int ind, bool all) const {
  if( (is <= ind) && (ie >= ind) )
    return true;
  if(all && (n != 0))
    return n->intersects(ind, all);
  return false;
}

RangeIterator& RangeIterator::operator=(const RangeIterator &r) {
  ind    = r.ind;
  is     = r.is;
  ie     = r.ie;
  n      = r.n;
  cur    = r.cur;
  if(cur == &r)
    cur = this;
  curend = r.curend;

  return *this;
}

RangeIterator& RangeIterator::operator+=(const RangeIterator &r) {
  // For now just put at the end
  RangeIterator *it = this;
  while(it->n != 0) {
    it = it->n;
  }
  it->n = new RangeIterator(r);
  return *this;
}

RangeIterator& RangeIterator::operator-=(const RangeIterator &r) {
  // Find any ranges which overlap
  RangeIterator *it = this;
  do{
    const RangeIterator *itr = &r;
    do {
      // Check if it and itr overlap
      if( it->intersects(*itr, false) ){
        // Overlap

        if( (itr->is <= it->is) && (itr->ie >= it->ie) ) {
          // Total overlap
          is = 1; ie = 0; // Make invalid
        }else if( (itr->is > it->is) && (itr->ie >= it->ie) ) {
          // Removing upper end
          it->ie = itr->is-1;
        }else if( (itr->is <= it->is) && (itr->ie < it->ie) ) {
          // Removing lower end
          it->is = itr->ie+1;
        }else {
          // Removing a chunk from the middle
          it->n = new RangeIterator(itr->ie+1, it->ie, it->n); // Upper piece
          it->ie = itr->is-1; // Truncate lower piece
        }
      }
      itr = itr->n;
    }while(itr != 0);

    // Check if this range is still valid
    if(is > ie) {
      // Invalid range
      if(it->n != 0) {
        // Copy from next one
        RangeIterator *tmp = it->n;
        *it = *it->n;
        // and delete the redundant object
        delete tmp;
        continue; // Check this range
      }
    }
    it = it->n;
  }while(it != 0);
  return *this;
}
