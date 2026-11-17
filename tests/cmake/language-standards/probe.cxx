#ifdef _MSVC_LANG
#define PROBE_CXX_STANDARD _MSVC_LANG
#else
#define PROBE_CXX_STANDARD __cplusplus
#endif

#if defined(CHECK_CXX23) && PROBE_CXX_STANDARD < 202100L
#error "The caller's C++23 standard was not preserved"
#endif
#if PROBE_CXX_STANDARD < 202002L
#error "The BOUT++ C++20 minimum was not applied"
#endif
int main() { return 0; }
