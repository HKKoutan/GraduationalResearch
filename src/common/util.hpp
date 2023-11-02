#ifndef INCLUDE_GUARD_common_util
#define INCLUDE_GUARD_common_util

#include <iostream>
#include <bitset>
#include <vector>
#include <cstdint>
#include <chrono>
#include <random>

namespace util {

class Timekeep {
	std::vector<std::chrono::milliseconds> times; //記録した時間
	std::chrono::system_clock::time_point prev; //計測開始時刻

public:
	inline void start() noexcept{//計測開始
		prev = std::chrono::system_clock::now();
	}
	inline void stop(){//計測終了
		times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - prev));
		std::cout<<"time"<<std::endl;
		for(auto &t: times) std::cout<<t.count()<<"ms"<<std::endl;
	}
	inline void split(){//今の時間を記録して新たに計測
		times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - prev));
		prev = std::chrono::system_clock::now();
	}
};

class RandomBits {
	std::mt19937_64 mt;
public:
	explicit RandomBits(std::int64_t seed): mt(seed){}
	RandomBits() noexcept: RandomBits(0i64){}
	template<std::size_t L>
	void generate(std::bitset<L> &bits);
};

template<std::size_t L>
void RandomBits::generate(std::bitset<L> &bits){
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