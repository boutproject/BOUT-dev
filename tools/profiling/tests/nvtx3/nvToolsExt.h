#pragma once
#include <cstdint>
#define NVTX_VERSION 3
#define NVTX_EVENT_ATTRIB_STRUCT_SIZE sizeof(nvtxEventAttributes_t)
#define NVTX_MESSAGE_TYPE_ASCII 1
#define NVTX_COLOR_ARGB 1
struct nvtxEventAttributes_t { int version, size, messageType, colorType; union { const char* ascii; } message; uint32_t color; };
int nvtxRangePushEx(const nvtxEventAttributes_t*);
int nvtxRangePop();
