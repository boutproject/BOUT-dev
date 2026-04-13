#pragma once

#ifndef BOUT_PETSC_JACOBIAN_H
#define BOUT_PETSC_JACOBIAN_H

#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include <petscmat.h>

/// Insert the nonzero pattern of @p sub into the variable block
/// (@p out_var, @p in_var) of the Jacobian @p Jfd.
///
/// @p Jfd is a square matrix of size (n_evolving * nvars) where nvars is
/// inferred as Jfd_global_size / sub_global_size.  Each nonzero (r, c) in
/// @p sub produces an entry at (r * nvars + out_var, c * nvars + in_var)
/// in @p Jfd.
///
/// @p Jfd must already be preallocated.  Entries are inserted with
/// INSERT_VALUES; MatAssemblyBegin/End must be called by the caller after
/// all insertions are complete.
///
/// @param Jfd     The Jacobian matrix to populate. Must be preallocated.
/// @param sub     Evolving-cell submatrix providing the nonzero pattern.
/// @param out_var Row variable index in [0, nvars).
/// @param in_var  Column variable index in [0, nvars).
void addOperatorSparsity(Mat Jfd, Mat sub, int out_var, int in_var);

/// @brief Insert the nonzero pattern of @p sub into every variable block of
///        @p Jfd.
///
/// Equivalent to calling addOperatorSparsity(Jfd, sub, out_var, in_var) for
/// every combination of @p out_var and @p in_var in [0, nvars), where
/// @c nvars is inferred as @c Jfd_global / @c sub_global.
///
/// @param Jfd     The Jacobian matrix to populate. Must be preallocated.
/// @param sub     Evolving-cell submatrix providing the nonzero pattern.
void addOperatorSparsity(Mat Jfd, Mat sub);

#endif // BOUT_HAS_PETSC

#endif // BOUT_PETSC_JACOBIAN_H
