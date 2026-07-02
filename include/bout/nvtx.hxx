#ifndef BOUT_NVTX_H
#define BOUT_NVTX_H

#include <nvtx3/nvToolsExt.h>
#include <cstdint>

namespace bout::profiling {

    inline void nvtxPushColor(const char* name, uint32_t argb)
    {
    nvtxEventAttributes_t ev{};
    ev.version = NVTX_VERSION;
    ev.size    = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    ev.messageType = NVTX_MESSAGE_TYPE_ASCII;
    ev.message.ascii = name;
    ev.colorType = NVTX_COLOR_ARGB;
    ev.color = argb;

    nvtxRangePushEx(&ev);
    }

    inline void nvtxPop()
    {
    nvtxRangePop();
    }

    namespace nvtxColor
    {
        constexpr uint32_t Red        = 0xFFFF0000;
        constexpr uint32_t Green      = 0xFF00FF00;
        constexpr uint32_t Blue       = 0xFF0000FF;
        constexpr uint32_t Yellow     = 0xFFFFFF00;
        constexpr uint32_t Cyan       = 0xFF00FFFF;
        constexpr uint32_t Magenta    = 0xFFFF00FF;

        constexpr uint32_t Orange     = 0xFFFFA500;
        constexpr uint32_t Purple     = 0xFF800080;
        constexpr uint32_t Teal       = 0xFF008080;

        constexpr uint32_t Gray       = 0xFF808080;
        constexpr uint32_t LightGray  = 0xFFD3D3D3;
        constexpr uint32_t DarkGray   = 0xFF404040;

        constexpr uint32_t Pink       = 0xFFFF69B4;
        constexpr uint32_t Lime       = 0xFF32CD32;
    }

} // namespace bout::profiling

#endif // BOUT_NVTX_H