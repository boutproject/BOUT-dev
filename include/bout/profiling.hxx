#pragma once

#include <cstdint>

#if defined(USE_NVTX) && defined(USE_ROCTX)
#error "Enable only one of USE_NVTX and USE_ROCTX"
#endif

#if defined(USE_NVTX)
#include <nvtx3/nvToolsExt.h>
#elif defined(USE_ROCTX)
#include <rocprofiler-sdk-roctx/roctx.h>
#endif

namespace bout::profiling {

// Shared by all host ranges on the current thread, including library callbacks.
inline unsigned& rangeDepth() {
  static thread_local unsigned depth = 0;
  return depth;
}

inline uint32_t depthColor(unsigned depth) {
  constexpr uint32_t colors[] = {0xFF0000FF, 0xFF00FFFF, 0xFF00FF00, 0xFFFFFF00,
                                 0xFFFFA500, 0xFFFF00FF, 0xFF800080};
  return colors[depth % 7];
}

inline void nvtxPushColor(const char* name, uint32_t argb) {
#if defined(USE_NVTX)
  nvtxEventAttributes_t ev{};
  ev.version = NVTX_VERSION;
  ev.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  ev.messageType = NVTX_MESSAGE_TYPE_ASCII;
  ev.message.ascii = name;
  ev.colorType = NVTX_COLOR_ARGB;
  ev.color = argb;
  nvtxRangePushEx(&ev);
  ++rangeDepth();
#elif defined(USE_ROCTX)
  (void)argb; // ROCTX ranges do not support colors.
  roctxRangePush(name);
#else
  (void)name;
  (void)argb;
#endif
}

inline void nvtxPop() {
#if defined(USE_NVTX)
  nvtxRangePop();
  --rangeDepth();
#elif defined(USE_ROCTX)
  roctxRangePop();
#endif
}

// Host scopes measure submission time without forcing GPU synchronization.
class ProfileRange {
public:
  explicit ProfileRange(const char* name) {
#if defined(USE_NVTX)
    nvtxPushColor(name, depthColor(rangeDepth()));
#elif defined(USE_ROCTX)
    roctxRangePush(name);
#else
    (void)name;
#endif
  }
  ~ProfileRange() { nvtxPop(); }
  ProfileRange(const ProfileRange&) = delete;
  ProfileRange& operator=(const ProfileRange&) = delete;
};

namespace nvtxColor {
constexpr uint32_t Red = 0xFFFF0000;
constexpr uint32_t Green = 0xFF00FF00;
constexpr uint32_t Blue = 0xFF0000FF;
constexpr uint32_t Yellow = 0xFFFFFF00;
constexpr uint32_t Cyan = 0xFF00FFFF;
constexpr uint32_t Magenta = 0xFFFF00FF;

constexpr uint32_t Orange = 0xFFFFA500;
constexpr uint32_t Purple = 0xFF800080;
constexpr uint32_t Teal = 0xFF008080;

constexpr uint32_t Gray = 0xFF808080;
constexpr uint32_t LightGray = 0xFFD3D3D3;
constexpr uint32_t DarkGray = 0xFF404040;

constexpr uint32_t Pink = 0xFFFF69B4;
constexpr uint32_t Lime = 0xFF32CD32;
} // namespace nvtxColor

} // namespace bout::profiling
