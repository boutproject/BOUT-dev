/************************************************************
 * Partition a rectangular domain into rectangles of equal area,
 * whilst minimising the length of the boundary
 *
 * Algorithm from:
 * 
 * Michelle L Wachs "Minimum Cost Partitions of a Rectangle"
 *        Discrete Mathematics 36 (1981) 305-324
 *
 ************************************************************/

#ifndef __PARTITION_H__
#define __PARTITION_H__

#include "domain.hxx"

/// Partition a single domain into n pieces
void partition(Domain* d, int n);

/// Partition domain and its neighbours into n pieces
void partitionAll(Domain* d, int n);

#endif // __PARTITION_H__
