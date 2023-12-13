#ifndef INCLUDE_GUARD_common_codecommon
#define INCLUDE_GUARD_common_codecommon

#include <array>
#include <bitset>
#include <concepts>

namespace code {

template<std::integral T, std::floating_point U>
__global__ void estimate_crop(T *est, U *LLR){
	std::size_t i = BlockIdx.x*blockDim.x+threadIdx.x;
	est[i] = LLR[i]<0;
}

// template<typename T, std::size_t L1, std::size_t L2>
// auto concatenate(const std::array<T,L1> arr1, const std::array<T,L2> arr2){
// 	std::array<T,L1+L2> arr;
// 	for(std::size_t i=0;i<L1;++i) arr[i]=arr1[i];
// 	for(std::size_t i=0;i<L2;++i) arr[i+L1]=arr2[i];
// 	return arr;
// }

}

#endif