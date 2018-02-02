#include <bout/region.hxx>

template<typename T>
Region<T>::Region(int xstart, int xend, int ystart, int yend,
		  int zstart, int zend, int ny, int nz){
  indices = createRegionIndices(xstart, xend, ystart, yend, zstart, zend, ny, nz);
  blocks  = getContiguousBlocks();
};

template<typename T>
Region<T>::Region(regionIndices& indices) : indices(indices){
  blocks = getContiguousBlocks();
}

/// Helper function to create a regionIndices, given the start and end
/// points in x, y, z, and the total y, z lengths
template<typename T>
typename Region<T>::regionIndices Region<T>::createRegionIndices(int xstart, int xend, int ystart, int yend,
								 int zstart, int zend, int ny, int nz) {
  
  int len = (xend - xstart + 1) * (yend - ystart + 1) * (zend - zstart + 1);
  regionIndices region(len);
  int j = 0;
  int x = xstart;
  int y = ystart;
  int z = zstart;
  
  bool done = false;
  j = -1;
  while (!done) {
    j++;
    region[j] = (x * ny + y) * nz + z;
    if (x == xend && y == yend && z == zend) {
      done = true;
    }
    ++z;
    if (z > zend) {
      z = zstart;
      ++y;
      if (y > yend) {
	y = ystart;
	++x;
      }
    }
  }
  return region;
}
    
    
/* 
 * Returns a vector of all contiguous blocks contained in the passed region.
 * Limits the maximum size of any contiguous block to maxBlockSize.
 * A contiguous block is described by the inclusive start and the exclusive end
 * of the contiguous block.
 */
template<typename T>
typename Region<T>::contiguousBlocks Region<T>::getContiguousBlocks() const{
  const int lastPoint = indices.size();
  contiguousBlocks result;
  int index = 0;
  
  while (index < lastPoint){
    const int startIndex = indices[index];
    int count = 1; //We will always have at least startPair in the block so count starts at 1
    
    //Consider if the next point should be added to this block
    for(index++; count<MAXREGIONBLOCKSIZE; index++){
      if((indices[index]-indices[index-1])==1){
	count++;
      }else{//Reached the end of this block so break
	break;
      }
    }
    
    //Add pair to output, denotes inclusive start and exclusive end
    result.push_back({startIndex, indices[index-1]});
  }
  
  return result;
}
