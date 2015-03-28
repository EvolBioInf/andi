#ifdef __cplusplus

std::vector<int> psufsort(const std::string& T);

#else

#ifdef __cplusplus
extern "C" {
#endif

int c_psufsort(const char *str, int* SA);

#ifdef __cplusplus
}
#endif

#endif
