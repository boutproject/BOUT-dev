#ifndef BOUT_GIT_METADATA_H
#define BOUT_GIT_METADATA_H

#include <string_view>

namespace bout {
namespace version {

/// True if git metadata could be collected for this build.
auto git_metadata_available() -> bool;

/// True if the source tree had tracked uncommitted changes at build time.
auto git_dirty() -> bool;

/// Diff of tracked changes in the source tree at build time.
auto git_diff() -> std::string_view;

} // namespace version
} // namespace bout

#endif // BOUT_GIT_METADATA_H
