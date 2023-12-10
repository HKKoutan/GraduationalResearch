#ifndef INCLUDE_GUARD_dnas_codeDNASadapter
#define INCLUDE_GUARD_dnas_codeDNASadapter

#include <array>
#include <bitset>
#include <cmath>
#include <limits>
#include "DNASnttype.hpp"

namespace code::DNAS {

struct differential {
	template<std::uint8_t ATGC, std::size_t S>
	static auto encode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::uint8_t ATGC, std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state = 0);
	template<std::uint8_t ATGC, std::floating_point T, std::size_t S>
	static auto decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state = {1,0,0,0});
};

template<std::uint8_t ATGC> struct convert;

template<>
class convert<0x1B> {
	static constexpr std::uint8_t ATGC = 0x1B;
public:
	template<std::size_t S>
	static auto nttype_to_binary(const std::array<nucleotide_t<ATGC>,S> &source);
	template<std::floating_point T, std::size_t S>
	static auto nttype_to_binary_p(const std::array<nucleotide_p<ATGC,T>,S> &source);
	template<std::size_t S>
	static auto binary_to_nttype(const std::bitset<S> &source);
};

template<>
class convert<0x27> {
	static constexpr std::uint8_t ATGC = 0x27;
public:
	template<std::size_t S>
	static auto nttype_to_binary(const std::array<nucleotide_t<ATGC>,S> &source);
	template<std::floating_point T, std::size_t S>
	static auto nttype_to_binary_p(const std::array<nucleotide_p<ATGC,T>,S> &source);
	template<std::size_t S>
	static auto binary_to_nttype(const std::bitset<S> &source);
};

////////////////////////////////////////////////////////////////
//                                                            //
//                     class differential                     //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint8_t ATGC, std::size_t S>
auto differential::encode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
	std::array<nucleotide_t<ATGC>,S> code;
	nucleotide_t current_state = initial_state;

	for(std::size_t i=0; i<S; ++i){
		current_state += source[i];
		code[i] = current_state;
	}
	return code;
}

template<std::uint8_t ATGC, std::size_t S>
auto differential::decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state){
	std::array<nucleotide_t<ATGC>,S> decode;
	nucleotide_t prev = initial_state;

	for(std::size_t i=0; i<S; ++i){
		decode[i] = code[i]-prev;
		prev = code[i];
	}
	return decode;
}

template<std::uint8_t ATGC, std::floating_point T ,std::size_t S>
auto differential::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state){
	std::array<nucleotide_p<ATGC,T>,S> decode;
	auto previous = initial_state;

	for(std::size_t i=0; i<S; ++i){
		const auto &current = code[i];
		nucleotide_p<ATGC,T> temp = {0,0,0,0};
		for(auto j=0ui8; j<4; ++j) for(auto k=0ui8; k<4; ++k) temp[j] += previous[k] * current[k+j];
		decode[i] = temp;
		previous = current;
	}
	return decode;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                       class convert                        //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S>
auto convert<0x1B>::nttype_to_binary(const std::array<nucleotide_t<ATGC>,S> &source){
	std::bitset<S*2> result;
	for(std::size_t i=0; const auto &j: source){
		result[i++]=j.msb();
		result[i++]=j.lsb();
	}
	return result;
}

template<std::floating_point T, std::size_t S>
auto convert<0x1B>::nttype_to_binary_p(const std::array<nucleotide_p<ATGC,T>,S> &source){
	std::array<T,S*2> LLR;
	for(std::size_t i=0; const auto &j: source){
		auto PX0 = j[0]+j[2], PX1 = j[1]+j[3], P0X = j[0]+j[1], P1X = j[2]+j[3];
		if(P0X==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else if(P1X==0) LLR[i] = std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(P0X)-std::log(P1X);
		++i;
		if(PX0==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else if(PX1==0) LLR[i] = std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(PX0)-std::log(PX1);
		++i;
	}
	return LLR;
}

template<std::size_t S>
auto convert<0x1B>::binary_to_nttype(const std::bitset<S> &source){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> result;
	for(std::size_t i=0; auto &j: result){
		j = static_cast<int>(source.test(i++))<<1;
		j += static_cast<int>(source.test(i++));
	}
	return result;
}

template<std::size_t S>
auto convert<0x27>::nttype_to_binary(const std::array<nucleotide_t<ATGC>,S> &source){
	std::bitset<S*2> result;
	for(std::size_t i=0; const auto &j: source){
		auto msb = j.msb();
		result[i++] = msb;
		result[i++] = j.lsb()^msb;
	}
	return result;
}

template<std::floating_point T, std::size_t S>
auto convert<0x27>::nttype_to_binary_p(const std::array<nucleotide_p<ATGC,T>,S> &source){
	std::array<T,S*2> LLR;
	for(std::size_t i=0; const auto &j: source){
		auto PX0 = j[0]+j[3], PX1 = j[1]+j[2], P0X = j[0]+j[1], P1X = j[2]+j[3];
		if(P0X==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else if(P1X==0) LLR[i] = std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(P0X)-std::log(P1X);
		++i;
		if(PX0==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else if(PX1==0) LLR[i] = std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(PX0)-std::log(PX1);
		++i;
	}
	return LLR;
}

template<std::size_t S>
auto convert<0x27>::binary_to_nttype(const std::bitset<S> &source){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> result;
	for(std::size_t i=0; auto &j: result){
		auto msb = source.test(i++);
		j = (static_cast<int>(msb)<<1)+static_cast<int>(source.test(i++)^msb);
	}
	return result;
}

}

#endif