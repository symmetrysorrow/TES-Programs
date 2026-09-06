#include <hip/hip_runtime.h>

#include <cstdio>

int main() {
  int count = 0;
  hipError_t status = hipGetDeviceCount(&count);
  if (status != hipSuccess) {
    std::fprintf(stderr, "hipGetDeviceCount failed: %s\n", hipGetErrorString(status));
    return 2;
  }
  std::printf("device_count=%d\n", count);
  for (int device = 0; device < count; ++device) {
    hipDeviceProp_t props{};
    status = hipGetDeviceProperties(&props, device);
    if (status != hipSuccess) {
      std::fprintf(stderr, "hipGetDeviceProperties(%d) failed: %s\n", device,
                   hipGetErrorString(status));
      return 3;
    }
    std::printf("device=%d name=%s gcnArchName=%s totalGlobalMem=%llu\n",
                device, props.name, props.gcnArchName,
                static_cast<unsigned long long>(props.totalGlobalMem));
    status = hipSetDevice(device);
    if (status != hipSuccess) {
      std::fprintf(stderr, "hipSetDevice(%d) failed: %s\n", device,
                   hipGetErrorString(status));
      return 4;
    }
    void *ptr = nullptr;
    status = hipMalloc(&ptr, 4096);
    if (status != hipSuccess) {
      std::fprintf(stderr, "hipMalloc failed: %s\n", hipGetErrorString(status));
      return 5;
    }
    status = hipMemset(ptr, 0, 4096);
    if (status != hipSuccess) {
      std::fprintf(stderr, "hipMemset failed: %s\n", hipGetErrorString(status));
      (void)hipFree(ptr);
      return 6;
    }
    status = hipDeviceSynchronize();
    (void)hipFree(ptr);
    if (status != hipSuccess) {
      std::fprintf(stderr, "hipDeviceSynchronize failed: %s\n", hipGetErrorString(status));
      return 7;
    }
    std::puts("device_memory_probe=PASS");
  }
  return count > 0 ? 0 : 8;
}
