#ifndef INCLUDE_GUARD_common_codecommon
#define INCLUDE_GUARD_common_codecommon

#include <array>
#include <bitset>
#include <concepts>

namespace code {

template<std::floating_point T, std::size_t L>
auto estimate(std::array<T,L> LLR){
	std::bitset<L> est;
	for(std::size_t i=0, i<L, ++i) est[i] = LLR[i]<static_cast<T>(0);
	return est;
}

}

#endif