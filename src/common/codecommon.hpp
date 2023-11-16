#ifndef INCLUDE_GUARD_common_codecommon
#define INCLUDE_GUARD_common_codecommon

#include <array>
#include <bitset>
#include <concepts>

namespace code {

template<std::size_t R, std::floating_point T, std::size_t L>
auto estimate_crop(std::array<T,L> LLR){
	std::bitset<R> est;
	for(std::size_t i=0; i<R; ++i) est[i] = LLR[i]<static_cast<T>(0);
	return est;
}

template<std::floating_point T, std::size_t L>
auto estimate(std::array<T,L> LLR){
	std::bitset<L> est;
	for(std::size_t i=0; i<L; ++i) est[i] = LLR[i]<static_cast<T>(0);
	return est;
}

template<typename T, std::size_t L1, std::size_t L2>
auto concatenate(const std::array<T,L1> arr1, const std::array<T,L2> arr2){
	std::array<T,L1+L2> arr;
	for(std::size_t i=0;i<L1;++i) arr[i]=arr1[i];
	for(std::size_t i=0;i<L2;++i) arr[i+L1]=arr2[i];
	return arr;
}

}

#endif