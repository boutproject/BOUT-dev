/* Host-only SUNDIALS ranges sharing BOUT++'s thread-local depth and palette.
 * Used by instrument_sundials.py; requires ELF weak symbols and GCC/Clang's
 * cleanup attribute. An uninstrumented application leaves these ranges inert.
 */
#ifndef SUNDIALS_BOUT_NVTX_H
#define SUNDIALS_BOUT_NVTX_H

#if !defined(__GNUC__) && !defined(__clang__)
#error "SUNDIALS BOUT NVTX scopes require GCC or Clang cleanup support"
#endif

extern void bout_nvtx_range_push(const char* name) __attribute__((weak));
extern void bout_nvtx_range_pop(void) __attribute__((weak));

static inline int sundials_bout_nvtx_begin(const char* name)
{
  if (!bout_nvtx_range_push || !bout_nvtx_range_pop) { return 0; }
  bout_nvtx_range_push(name);
  return 1;
}

static inline void sundials_bout_nvtx_end(int* active)
{
  if (*active) { bout_nvtx_range_pop(); }
}

/* Cleanup also runs on return, break, continue, and goto out of a scope. */
#define SUNDIALS_BOUT_NVTX_SCOPE(variable, name) \
  int variable __attribute__((cleanup(sundials_bout_nvtx_end), unused)) = \
    sundials_bout_nvtx_begin(name)

#endif
