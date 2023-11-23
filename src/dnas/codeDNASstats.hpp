#ifndef INCLUDE_GUARD_dnas_codeDNASstats
#define INCLUDE_GUARD_dnas_codeDNASstats

#include <array>
#include "DNASnttype.hpp"

namespace code::DNAS {

////////////////////////////////////////////////////////////////
//                                                            //
//                         statistics                         //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint8_t ATGC, std::size_t L>
auto countAT(const std::array<nucleotide_t<ATGC>,L> &c){
	std::uint64_t qty_AT=0;
	for(const auto &i: c) qty_AT+=i.is_AT();
	return qty_AT;
}

template<std::size_t BS = 0, std::uint8_t ATGC, std::size_t L>
auto countBlockGCmaxDeviation(const std::array<nucleotide_t<ATGC>,L> &c){
	constexpr std::size_t block_size = BS==0?L:BS;
	std::uint64_t maxDeviation = 0;

	for(std::size_t i=0u, iend=L/block_size; i<iend; ++i){
		section block(i*block_size, block_size);
		std::int64_t blockGC = 0;
		for(std::size_t j=block.head, jend=block.tail; j<jend; ++j){
			blockGC += c[j].is_GC();
		}
		std::uint64_t deviation = static_cast<std::uint64_t>(std::abs(static_cast<std::int64_t>(block_size>>1)-blockGC));
		if(deviation>maxDeviation) maxDeviation=deviation;
	}
	return maxDeviation;
}

template<std::uint8_t ATGC, std::size_t L>
auto countRunlength(const std::array<nucleotide_t<ATGC>,L> &c){
	std::size_t rlmax=0, rl=0;
	nucleotide_t<ATGC> prev;
	for(const auto &i: c){
		if(i!=prev){
			rlmax = (rl>rlmax)?rl:rlmax;
			rl=0;
		}else{
			++rl;
		}
	}
	return (rl>rlmax)?rl:rlmax;
}

template<std::uint8_t ATGC, std::size_t L>
auto countError(const std::array<nucleotide_t<ATGC>,L> &c1, const std::array<nucleotide_t<ATGC>,L> &c2){
	std::uint64_t errorcount = 0;
	for(std::size_t i=0; i<L; ++i){
		if(c1[i]!=c2[i]) ++errorcount;
	}
	return errorcount;
}

template<std::uint8_t ATGC, std::size_t L>
auto countDifferentialError(const std::array<nucleotide_t<ATGC>,L> &c1, const std::array<nucleotide_t<ATGC>,L> &c2){
	std::uint64_t errorcount = 0;
	nucleotide_t<ATGC> prev1=0, prev2=0;
	for(std::size_t i=0; i<L; ++i){
		auto current1 = c1[i], current2 = c2[i];
		if(current1-prev1!=current2-prev2) ++errorcount;
		prev1 = current1;
		prev2 = current2;
	}
	return errorcount;
}

}

#endif