#ifndef INCLUDE_GUARD_dnas_codeDNASstats
#define INCLUDE_GUARD_dnas_codeDNASstats

#include <array>
#include "DNASnttype.hpp"

namespace code::DNAS {

template<std::uint8_t ATGC, std::size_t L>
auto countAT(const std::array<nucleotide_t<ATGC>,L> &c);
template<std::uint8_t ATGC, std::size_t L>
auto countRunlength(const std::array<nucleotide_t<ATGC>,L> &c);
template<std::uint8_t ATGC, std::size_t L>
auto countError(const std::array<nucleotide_t<ATGC>,L> &c1, const std::array<nucleotide_t<ATGC>,L> &c2);
template<std::uint8_t ATGC, std::size_t L>
auto countDifferentialError(const std::array<nucleotide_t<ATGC>,L> &c1, const std::array<nucleotide_t<ATGC>,L> &c2);

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