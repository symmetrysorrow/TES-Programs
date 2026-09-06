#include <hipsparse/hipsparse.h>

#include <cstdio>

int main() {
  hipsparseHandle_t handle = nullptr;
  hipsparseStatus_t status = hipsparseCreate(&handle);
  if (status != HIPSPARSE_STATUS_SUCCESS) {
    std::fprintf(stderr, "hipsparseCreate failed: %d\n", static_cast<int>(status));
    return 2;
  }
  int version = 0;
  status = hipsparseGetVersion(handle, &version);
  if (status != HIPSPARSE_STATUS_SUCCESS) {
    std::fprintf(stderr, "hipsparseGetVersion failed: %d\n", static_cast<int>(status));
    hipsparseDestroy(handle);
    return 3;
  }
  std::printf("hipsparse_handle=PASS version=%d\n", version);
  status = hipsparseDestroy(handle);
  return status == HIPSPARSE_STATUS_SUCCESS ? 0 : 4;
}
