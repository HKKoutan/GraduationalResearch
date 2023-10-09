#ifndef __util_RANDOMBITS__
#define __util_RANDOMBITS__

#include <random>
#include <bitset>
#include <vector>
#include <cstdint>

namespace util {

template<std::size_t L>
class RandomBits {
	std::mt19937_64 mt;
public:
	explicit RandomBits(std::int64_t seed): mt(seed){}
	RandomBits() noexcept: RandomBits(0i64){}
	void generate(std::bitset<L> &bits);
};

template<std::size_t L>
void RandomBits<L>::generate(std::bitset<L> &bits){
	constexpr auto L64 = L/64u;
	for(std::size_t i=0u, iend=L64; i<iend; ++i){
		std::bitset<64u> rand(mt());
		for(std::size_t j=0u, jend=64u, pos=i*64u; j<jend; ++j) bits[pos+j]=rand[j];
	}
	if constexpr(L%64u!=0u){
		std::bitset<L%64u> rand(mt());
		for(std::size_t j=0u, jend=L%64u, pos=L64*64u; j<jend; ++j) bits[pos+j]=rand[j];
	}
}

}

#endif