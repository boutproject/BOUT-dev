#include "bout/build_defines.hxx"

#if BOUT_HAS_CVODE

#include "gtest/gtest.h"

#include "../../../src/solver/impls/cvode/cvode.hxx"

#include <sundials/sundials_linearsolver.h>

using bout::classifyCvodeKrylovFlag;
using bout::CvodeKrylovOutcome;

// Classification of Krylov return flags for solver:cvode_failure_counters, which
// must agree with cvLsSolve in SUNDIALS cvode_ls.c

TEST(CvodeKrylovFlagTest, SuccessIsConverged) {
  EXPECT_EQ(classifyCvodeKrylovFlag(0, 0), CvodeKrylovOutcome::converged);
  EXPECT_EQ(classifyCvodeKrylovFlag(0, 2), CvodeKrylovOutcome::converged);
}

TEST(CvodeKrylovFlagTest, ReducedResidualAcceptedOnlyInFirstNewtonIteration) {
  EXPECT_EQ(classifyCvodeKrylovFlag(SUNLS_RES_REDUCED, 0), CvodeKrylovOutcome::accepted);
  EXPECT_EQ(classifyCvodeKrylovFlag(SUNLS_RES_REDUCED, 1),
            CvodeKrylovOutcome::newton_failure);
  EXPECT_EQ(classifyCvodeKrylovFlag(SUNLS_RES_REDUCED, 2),
            CvodeKrylovOutcome::newton_failure);
}

TEST(CvodeKrylovFlagTest, RecoverableFailuresFailNewtonIteration) {
  for (const int flag : {SUNLS_CONV_FAIL, SUNLS_ATIMES_FAIL_REC, SUNLS_PSOLVE_FAIL_REC,
                         SUNLS_PACKAGE_FAIL_REC, SUNLS_QRFACT_FAIL, SUNLS_LUFACT_FAIL}) {
    EXPECT_EQ(classifyCvodeKrylovFlag(flag, 0), CvodeKrylovOutcome::newton_failure)
        << "flag " << flag;
    EXPECT_EQ(classifyCvodeKrylovFlag(flag, 1), CvodeKrylovOutcome::newton_failure)
        << "flag " << flag;
  }
}

TEST(CvodeKrylovFlagTest, OtherFlagsAreUnrecoverable) {
  for (const int flag : {SUNLS_ATIMES_FAIL_UNREC, SUNLS_PSOLVE_FAIL_UNREC, SUNLS_GS_FAIL,
                         SUNLS_QRSOL_FAIL, -1}) {
    EXPECT_EQ(classifyCvodeKrylovFlag(flag, 0), CvodeKrylovOutcome::unrecoverable)
        << "flag " << flag;
    EXPECT_EQ(classifyCvodeKrylovFlag(flag, 1), CvodeKrylovOutcome::unrecoverable)
        << "flag " << flag;
  }
}

#endif
