#ifndef INCLUDE_GUARD_dnas_codeDNASadapter
#define INCLUDE_GUARD_dnas_codeDNASadapter

#include <array>
#include <bitset>
#include <cmath>
#include <limits>
#include "DNASnttype.hpp"

namespace code::DNAS {

template<std::uint8_t ATGC> struct differential;

template<>
class differential<0x1B> {
	static constexpr std::uint8_t ATGC = 0x1B;
public:
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state = 0);
	template<std::floating_point T, std::size_t S>
	static auto decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state = {1,0,0,0});
};

template<>
class differential<0x27> {
	static constexpr std::uint8_t ATGC = 0x27;
public:
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state = 0);
	template<std::floating_point T, std::size_t S>
	static auto decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state = {1,0,0,0});
};

template<std::uint8_t ATGC>
struct convert {
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

template<std::size_t S>
auto differential<0x1B>::encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code;
	nucleotide_t current_state = initial_state;

	for(std::size_t i=0; i<S; i+=2){
		current_state += (static_cast<int>(source.test(i))<<1)+static_cast<int>(source.test(i+1));
		code[i>>1] = current_state;
	}
	return code;
}

template<std::size_t S>
auto differential<0x1B>::decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state){
	std::bitset<S*2> decode;
	nucleotide_t prev = initial_state;

	for(std::size_t i=0; i<S; ++i){
		auto j=i<<1;
		auto current = code[i]-prev;
		decode[j] = current.msb();
		decode[j+1] = current.lsb();
		prev = code[i];
	}
	return decode;
}

template<std::floating_point T ,std::size_t S>
auto differential<0x1B>::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state){
	std::array<T,S*2> LLR;
	auto previous = initial_state;

	for(std::size_t i=0; const auto &current: code){
		T PX0 = 0, PX1 = 0, P0X = 0, P1X = 0;
		for(auto j=0ui8; j<4; ++j) P0X += previous[j] * (current[j] + current[j+1]); //上位ビットが0: 遷移語が0 or 1になる組み合わせ
		for(auto j=0ui8; j<4; ++j) P1X += previous[j] * (current[j+2] + current[j+3]); //上位ビットが1: 遷移語が2 or 3になる組み合わせ
		for(auto j=0ui8; j<4; ++j) PX0 += previous[j] * (current[j] + current[j+2]); //下位ビットが0: 遷移語が0 or 2になる組み合わせ
		for(auto j=0ui8; j<4; ++j) PX1 += previous[j] * (current[j+1] + current[j+3]); //下位ビットが1: 遷移語が1 or 3になる組み合わせ
		previous = current;
		if(P0X==0) LLR[i] = std::numeric_limits<T>::infinity();
		else if(P1X==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(P0X)-std::log(P1X);
		++i;
		if(PX0==0) LLR[i] = std::numeric_limits<T>::infinity();
		else if(PX1==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(PX0)-std::log(PX1);
		++i;
	}
	return LLR;
}

template<std::size_t S>
auto differential<0x27>::encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code;
	nucleotide_t current_state = initial_state;

	for(std::size_t i=0; i<S; i+=2){
		auto msb = source.test(i);
		current_state += (static_cast<int>(msb)<<1)+static_cast<int>(source.test(i+1)^msb);
		code[i>>1] = current_state;
	}
	return code;
}

template<std::size_t S>
auto differential<0x27>::decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state){
	std::bitset<S*2> decode;
	nucleotide_t prev = initial_state;

	for(std::size_t i=0; i<S; ++i){
		auto j=i<<1;
		auto current = code[i]-prev;
		auto msb = current.msb();
		decode[j] = msb;
		decode[j+1] = current.lsb()^msb;
		prev = code[i];
	}
	return decode;
}

template<std::floating_point T ,std::size_t S>
auto differential<0x27>::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state){
	std::array<T,S*2> LLR;
	auto previous = initial_state;

	for(std::size_t i=0; const auto &current: code){
		T PX0 = 0, PX1 = 0, P0X = 0, P1X = 0;
		for(auto j=0ui8; j<4; ++j) P0X += previous[j] * (current[j] + current[j+1]); //上位ビットが0: 遷移語が0 or 1になる組み合わせ
		for(auto j=0ui8; j<4; ++j) P1X += previous[j] * (current[j+2] + current[j+3]); //上位ビットが1: 遷移語が2 or 3になる組み合わせ
		for(auto j=0ui8; j<4; ++j) PX0 += previous[j] * (current[j] + current[j+3]); //下位ビットが0: 遷移語が0 or 3になる組み合わせ
		for(auto j=0ui8; j<4; ++j) PX1 += previous[j] * (current[j+1] + current[j+2]); //下位ビットが1: 遷移語が1 or 2になる組み合わせ
		previous = current;
		if(P0X==0) LLR[i] = std::numeric_limits<T>::infinity();
		else if(P1X==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(P0X)-std::log(P1X);
		++i;
		if(PX0==0) LLR[i] = std::numeric_limits<T>::infinity();
		else if(PX1==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(PX0)-std::log(PX1);
		++i;
	}
	return LLR;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                       class convert                        //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint8_t ATGC>
template<std::size_t S>
auto convert<ATGC>::nttype_to_binary(const std::array<nucleotide_t<ATGC>,S> &source){
	std::bitset<S*2> code;
	for(std::size_t i=0; const auto &j: source){
		code[i++]=j.msb();
		code[i++]=j.lsb();
	}
	return code;
}

template<std::uint8_t ATGC>
template<std::floating_point T, std::size_t S>
auto convert<ATGC>::nttype_to_binary_p(const std::array<nucleotide_p<ATGC,T>,S> &source){
	std::array<T,S*2> LLR;
	for(std::size_t i=0; const auto &j: source){
		auto PX0 = j[0]+j[2], PX1 = j[1]+j[3], P0X = j[0]+j[1], P1X = j[2]+j[3];
		if(P0X==0) LLR[i] = std::numeric_limits<T>::infinity();
		else if(P1X==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(P0X)-std::log(P1X);
		++i;
		if(PX0==0) LLR[i] = std::numeric_limits<T>::infinity();
		else if(PX1==0) LLR[i] = -std::numeric_limits<T>::infinity();
		else LLR[i] = std::log(PX0)-std::log(PX1);
		++i;
	}
	return LLR;
}

template<std::uint8_t ATGC>
template<std::size_t S>
auto convert<ATGC>::binary_to_nttype(const std::bitset<S> &source){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code;
	for(std::size_t i=0; auto &j: code){
		j = static_cast<int>(source.test(i++))<<1;
		j += static_cast<int>(source.test(i++));
	}
	return code;
}

}

#endif