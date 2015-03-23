#include <string>
#include <vector>
#include <utility>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>
#include <global.h>

void mk_sort (std::vector<int>& SA, const std::string& T, size_t l, size_t r, size_t depth);
void insertion_sort (std::vector<int>& SA, const std::string& T, size_t l, size_t r, size_t depth);
void TSQS (std::vector<int>& SA, const std::string& T, size_t l, size_t r, size_t depth);
void mk_buildin (std::vector<int>& SA, const std::string& T, size_t l, size_t r, size_t depth);

std::vector<int> psufsort(const std::string& T);

class Bucket {

public:
	size_t start, size;
	constexpr Bucket() noexcept : start(0), size(0) {};
};

class PSufSort
{
	const std::string& T;
	std::vector<int>& SA;
	size_t threshold;
public:
	PSufSort(const std::string& _T, std::vector<int>& _SA, size_t size) : T(_T), SA(_SA) {
		threshold = std::log(size);
	}
	~PSufSort() {};

	// All intervals are semi-open: [l,r)
	void sort(size_t l, size_t r, size_t depth, size_t calls);
	void sort_tsqs(size_t l, size_t r, size_t depth, size_t calls);
	void sort_insert(size_t l, size_t r, size_t depth, size_t calls);
	void sort_heap(size_t l, size_t	r, size_t depth, size_t calls);

	void build_heap( int* rSA, size_t n, size_t depth);
	void heapify( int* rSA, size_t heap_size, size_t i, size_t depth);

	char median3(size_t a, size_t b, size_t c, size_t depth);

	void swap_range(size_t a, size_t b, size_t n);
	inline char char_at( size_t sai, size_t depth);
	inline const char *str_from( size_t sai, size_t depth);
};

std::vector<int> psufsort(const std::string& T){
	auto n = T.size();
	auto SA = std::vector<int>(n+1);

	auto L = std::vector<Bucket>(256);
	auto bucket_SM = std::vector<Bucket>(256*256); // S-
	auto bucket_SS = std::vector<Bucket>(256*256); // S*

	auto SM = [&](size_t X, size_t Y) -> Bucket& {
		return bucket_SM[(X<<8) + Y];
	};
	auto SS = [&](size_t X, size_t Y) -> Bucket& {
		return bucket_SS[(X<<8) + Y];
	};

	// classify suffixes and compute the bucket sizes
	ssize_t i = n -1;
	while( i >= 0){
		if( T[i] >= T[i+1]){
			L[T[i]].size++;
			i--;
			continue;
		}

		SS(T[i], T[i+1]).size++;
		i--;

		while( i >= 0 && T[i] <= T[i+1]){
			SM(T[i], T[i+1]).size++;
			i--;
		}
	}

	// Deal with the null byte
	SS(0,0).size = 1;
	SA[0] = n;

	// compute bucket starting points
	int j = 0;
	for(i=0; i<256; i++){
		L[i].start = j;
		j += L[i].size;
		for(auto k =i; k< 256; k++){
			SS(i, k).start = j;
			j += SS(i, k).size;
			SM(i, k).start = j;
			j += SM(i, k).size;
		}
	}

	// fill in S* buckets
	i = n -1;
	while( i >= 0){
		auto c = T[i];

		if( c >= T[i+1]){ // skip L
			i--;
			continue;
		}

		SA[SS(c,T[i+1]).start] = i;
		SS(c,T[i+1]).start++;
		i--;

		while( i >= 0 && T[i] <= T[i+1]){ // skip S-
			i--;
		}
	}

	// correct the `++` from the previous loop. 
	for(i=0; i<256*256; i++){
		bucket_SS[i].start -= bucket_SS[i].size;
	}

	// sort all S* suffixes
	#pragma omp parallel for shared(SA,T) schedule(dynamic, 1) num_threads( (FLAGS & F_LOW_MEMORY) ? THREADS : 1)
	for(i=0; i<256*256; i++){
		const auto buc = bucket_SS[i];
		if( buc.size > 1){
			int b = buc.start;
			int e = b + buc.size;
			auto sorter = PSufSort( T, SA, buc.size);

			// sort
			sorter.sort( b, e, 2, 0);
		}
	}

	// induced insert all S-
	for(i=n; i >= 0;i--){
		int j = SA[i];
		if( j == 0) continue;

		auto a = T[j-1];
		auto b = T[j];
		if( a <= b){
			SA[SM(a, b).start + SM(a, b).size - 1] = j-1;
			SM(a, b).size--;
		}
	}

	// induced insert all Ls
	for(i=0; i<n+1; i++){
		int j = SA[i];
		if( j == 0) continue;

		auto a = T[j-1];
		if( SA[L[a].start] != 0) continue;
		if( a >= T[j]){
			SA[L[a].start] = j -1;
			L[a].start++;
		}
	}

	return SA; // move doesnt move
}

inline char PSufSort::char_at( size_t sai, size_t depth){
	return T[sai + depth];
}

inline const char *PSufSort::str_from( size_t sai, size_t depth){
	return &T[sai + depth];
}

void PSufSort::swap_range(size_t a, size_t b, size_t n){
	auto& SA = this->SA;

	for(auto i=0; i< n; i++){
		std::swap(SA[a++], SA[b++]);
	}
}

void PSufSort::sort (size_t l, size_t r, size_t depth, size_t calls) {
	if(l >= r){
		return;
	}

	auto m = r - l;
	if( m < 2 ){
		return;
	}

	if (m <= 16){
		sort_insert(l, r, depth, calls);
		return;
	}

	if( calls < threshold){
		sort_tsqs(l, r, depth, calls);
	} else {
		sort_heap(l, r, depth, calls);
	}
}


void PSufSort::sort_insert (size_t l, size_t r, size_t depth, size_t unused){
	auto cmp_from = [&]( size_t a, size_t b){
		auto ta = T.data()+ a +depth;
		auto tb = T.data()+ b +depth;

		return strcmp(ta,tb);
	};

	for(auto j = l+1; j < r; j++){
		auto X = SA[j];

		auto i = j;
		for(; i > l && cmp_from( SA[i-1], X) > 0 ; i--){
			SA[i] = SA[i-1];
		}

		SA[i] = X;
	}
}

char PSufSort::median3(size_t b, size_t a, size_t c, size_t depth){
	auto key = [&](size_t i){
		return char_at(SA[i], depth);
	};

	if( key(a) > key(b) ){ std::swap(SA[a], SA[b]); }
	if( key(b) > key(c) ){ std::swap(SA[b], SA[c]); }
	if( key(a) > key(b) ){ std::swap(SA[a], SA[b]); }

	return key(b);
}

void PSufSort::sort_tsqs (size_t l, size_t r, size_t depth, size_t calls){
	auto K = median3(l, (r-l)/2 + l, r-1, depth); // pick K

	auto a = l;
	auto b = l;
	auto c = r-1;
	auto d = r-1;

	while(true) {
		for(; b <= c && char_at(SA[b], depth) <= K; b++){
			if( char_at(SA[b], depth) == K){
				std::swap(SA[a], SA[b]);
				a++;
			}
		}

		for(; b <= c && char_at(SA[c], depth) >= K; c--){
			if( char_at(SA[c], depth) == K){
				std::swap(SA[c],SA[d]);
				d--;
			}
		}

		if( b > c) break;

		std::swap(SA[b],SA[c]);
		b++, c--;
	}

	auto m = std::min(a-l, b-a);
	swap_range(l, b-m, m);

	m = std::min(d-c, r-d-1);
	swap_range(b, r-m, m);

	auto i = l + b - a;
	auto j = r - d + c;

	this->sort(l, i, depth, calls + 1);
	this->sort(i, j, depth+1, calls + 1);
	this->sort(j, r, depth, calls + 1);
}

constexpr inline size_t LEFT(size_t i) noexcept {
	return (i << 1) + 1;
}

constexpr inline size_t RIGHT(size_t i) noexcept {
	return (i << 1) + 2;
}

constexpr inline size_t PARENT( size_t i) noexcept {
	return (i-1) >> 1;
}

void PSufSort::build_heap( int* rSA, size_t n, size_t depth){
	auto heap_size = n;
	for( ssize_t i= PARENT(n-1); i>=0 ; i--){
		heapify(rSA, heap_size, i, depth);
	}
}

void PSufSort::heapify( int* rSA, size_t heap_size, size_t i, size_t depth){ // aka. siftDown
	auto key = [&](size_t j){
		return T.data() + j + depth;
	};

	auto l = LEFT(i);
	auto r = RIGHT(i);
	auto largest = i;

	if( l < heap_size && strcmp( key(rSA[l]) , key(rSA[i])) > 0 ){
		largest = l;
	}
	if( r < heap_size && strcmp( key(rSA[r]) , key(rSA[largest])) > 0){
		largest = r;
	}
	if( largest != i){
		std::swap(rSA[i], rSA[largest]);
		heapify(rSA, heap_size, largest, depth);
	}
}

void PSufSort::sort_heap(size_t left, size_t right, size_t depth, size_t calls){
	auto rSA = SA.data() + left;
	auto n = right - left;
	auto heap_size = n;
	
	build_heap(rSA, n, depth);

	for( auto i = n-1; i ; i--){
		std::swap(rSA[0],rSA[i]);
		heap_size--;
		heapify(rSA, heap_size, 0, depth);
	}
}

