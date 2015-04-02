#include <string>
#include <vector>
#include <cstring>

std::vector<int> psufsort(const std::string& T);

extern "C" int c_psufsort(const char *str, int* SA){
	if( !str || !SA){
		return 1;
	}
	auto T = std::string(str);
	auto temp = psufsort(T);
	memcpy(SA, temp.data()+1, T.size() * sizeof(int));
	return 0;
}
