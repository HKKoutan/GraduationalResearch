#ifndef INCLUDE_GUARD_dnas_codeDNASRLL
#define INCLUDE_GUARD_dnas_codeDNASRLL

#include <array>
#include <bitset>
#include <cmath>
#include <limits>
#include "DNASnttype.hpp"

namespace code::DNAS {

template<std::uint8_t ATGC> struct VLRLL;

template<>
struct VLRLL<0x1B> {
	static constexpr std::uint8_t ATGC = 0x1B;
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state = 0);
};

template<>
struct VLRLL<0x27> {
	static constexpr std::uint8_t ATGC = 0x27;
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state = 0);
};

template<std::uint8_t ATGC> struct modifiedVLRLL;

template<>
struct modifiedVLRLL<0x1B> {
	static constexpr std::uint8_t ATGC = 0x1B;
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state, std::size_t initial_runlength);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state);
	template<std::floating_point T, std::size_t S>
	static auto decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state);
};

template<>
struct modifiedVLRLL<0x27> {
	static constexpr std::uint8_t ATGC = 0x27;
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state, std::size_t initial_runlength);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state);
	template<std::floating_point T, std::size_t S>
	static auto decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state);
};

////////////////////////////////////////////////////////////////
//                                                            //
//                        class VLRLL                         //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S>
auto VLRLL<0x1B>::encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code;
	std::size_t processing = 0;
	nucleotide_t current_state = initial_state;

	std::size_t i=0 ,j=0;
	while(j<S/2){
		++processing;
		switch(processing){
		// case 1:
		// case 3:
		// 	break;

		case 2:
		case 4:
			if(source.test(i+processing-2)){
				current_state += (source.test(i+processing-1)?0:3);
				if(!source.test(i+processing-1)){
					i+=processing;
					processing=0;
				}
			}else{
				current_state += (source.test(i+processing-1)?2:1);
				i+=processing;
				processing=0;
			}
			code[j++] = current_state;
			break;

		case 5:
			if(source.test(i+4)){
				current_state += 3;
				code[j++] = current_state;
				i+=processing;
				processing=0;
			}
			break;

		case 6:
			current_state = current_state+(source.test(i+processing-1)?2:1);
			code[j++] = current_state;
			i+=processing;
			processing=0;
			break;
		}
	}

	std::bitset<S> used;
	used.set();
	used >>= S-(i+processing);

	return std::make_tuple(code,used,processing>>1);
}

template<std::size_t S>
auto VLRLL<0x1B>::decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
	std::bitset<S*2> decode;
	nucleotide_t previous = initial_state;
	std::size_t zeros = 0;

	for(std::size_t i=0, j=0; i<S; ++i){
		switch(source[i]-previous){
		case 0:
			decode.set(j++);
			decode.set(j++);
			++zeros;
			break;

		case 1:
			decode.reset(j++);
			decode.reset(j++);
			zeros=0;
			break;

		case 2:
			decode.reset(j++);
			decode.set(j++);
			zeros=0;
			break;

		case 3:
			decode.set(j++);
			if(zeros<2) decode.reset(j++);
			zeros=0;
			break;
		}
		previous = source[i];
	}
	return decode;
}

template<std::size_t S>
auto VLRLL<0x27>::encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code;
	std::size_t processing = 0;
	nucleotide_t current_state = initial_state;

	std::size_t i=0 ,j=0;
	while(j<S/2){
		++processing;
		switch(processing){
		// case 1:
		// case 3:
		// 	break;

		case 2:
		case 4:
			if(source.test(i+processing-2)){
				current_state += (source.test(i+processing-1)?0:3);
				if(!source.test(i+processing-1)){
					i+=processing;
					processing=0;
				}
			}else{
				current_state += (source.test(i+processing-1)?1:2);
				i+=processing;
				processing=0;
			}
			code[j++] = current_state;
			break;

		case 5:
			if(source.test(i+4)){
				current_state += 3;
				code[j++] = current_state;
				i+=processing;
				processing=0;
			}
			break;

		case 6:
			current_state = current_state+(source.test(i+processing-1)?1:2);
			code[j++] = current_state;
			i+=processing;
			processing=0;
			break;
		}
	}

	std::bitset<S> used;
	used.set();
	used >>= S-(i+processing);

	return std::make_tuple(code,used,processing>>1);
}

template<std::size_t S>
auto VLRLL<0x27>::decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
	std::bitset<S*2> decode;
	nucleotide_t previous = initial_state;
	std::size_t zeros = 0;

	for(std::size_t i=0, j=0; i<S; ++i){
		switch(source[i]-previous){
		case 0:
			decode.set(j++);
			decode.set(j++);
			++zeros;
			break;

		case 1:
			decode.reset(j++);
			decode.set(j++);
			zeros=0;
			break;

		case 2:
			decode.reset(j++);
			decode.reset(j++);
			zeros=0;
			break;

		case 3:
			decode.set(j++);
			if(zeros<2) decode.reset(j++);
			zeros=0;
			break;
		}
		previous = source[i];
	}
	return decode;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                     class modifiedVLRLL                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S>
auto modifiedVLRLL<0x1B>::encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state, std::size_t initial_runlength){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code; 
	std::size_t processing = initial_runlength<<1;
	nucleotide_t current_state = initial_state;

	std::size_t i=(~processing)+1, j=0;
	while(j<S/2){
		++processing;
		switch(processing){
		// case 0:
		// case 1:
		// case 3:
		// case 5:
		// 	break;

		case 2:
		case 4:
			if(source.test(i+processing-2)){
				current_state += (source.test(i+processing-1)?0:3);
				if(!source.test(i+processing-1)){
					i+=processing;
					processing=0;
				}
			}else{
				current_state += (source.test(i+processing-1)?2:1);
				i+=processing;
				processing=0;
			}
			code[j++] = current_state;
			break;

		case 6:
			if(source.test(i+processing-2)){
				current_state += 3;
			}else{
				current_state += (source.test(i+processing-1)?2:1);
			}
			i+=processing;
			processing=0;
			code[j++] = current_state;
			break;
		}
	}
	return code;
}

template<std::size_t S>
auto modifiedVLRLL<0x1B>::decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state){
	std::bitset<S*2> decode;
	nucleotide_t previous = initial_state;

	for(std::size_t i=0, j=0; i<S; ++i){
		switch(code[i]-previous){
		case 0:
			decode.set(j++);
			decode.set(j++);
			break;

		case 1:
			decode.reset(j++);
			decode.reset(j++);
			break;

		case 2:
			decode.reset(j++);
			decode.set(j++);
			break;

		case 3:
			decode.set(j++);
			decode.reset(j++);
			break;
		}
		previous = code[i];
	}
	return decode;
}

template<std::floating_point T, std::size_t S>
auto modifiedVLRLL<0x1B>::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state){
	constexpr auto p3_to_11 = static_cast<T>(1.0/17.0);
	constexpr auto p3_to_10 = static_cast<T>(16.0/17.0);
	std::array<T,S*2> LLR;
	auto previous = initial_state;

	for(std::size_t i=0; const auto &current: code){
		T PX0 = 0, PX1 = 0, P0X = 0, P1X = 0;
		for(auto j=0ui8; j<4; ++j) P0X += previous[j] * (current[j+1] + current[j+2]); //上位ビットが0: 遷移語が1 or 2になる組み合わせ
		for(auto j=0ui8; j<4; ++j) P1X += previous[j] * (current[j] + current[j+3]); //上位ビットが1: 遷移語が0 or 3になる組み合わせ
		for(auto j=0ui8; j<4; ++j) PX0 += previous[j] * (current[j+1] + p3_to_10*current[j+3]); //下位ビットが0: 遷移語が1 or 3(*p3_to_10)になる組み合わせ
		for(auto j=0ui8; j<4; ++j) PX1 += previous[j] * (current[j] + current[j+2] + p3_to_11*current[j+3]); //下位ビットが1: 遷移語が0 or 2 or 3(*p3_to_11)になる組み合わせ
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
auto modifiedVLRLL<0x27>::encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state, std::size_t initial_runlength){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code; 
	std::size_t processing = initial_runlength<<1;
	nucleotide_t current_state = initial_state;

	std::size_t i=(~processing)+1, j=0;
	while(j<S/2){
		++processing;
		switch(processing){
		// case 0:
		// case 1:
		// case 3:
		// case 5:
		// 	break;

		case 2:
		case 4:
			if(source.test(i+processing-2)){
				current_state += (source.test(i+processing-1)?0:3);
				if(!source.test(i+processing-1)){
					i+=processing;
					processing=0;
				}
			}else{
				current_state += (source.test(i+processing-1)?1:2);
				i+=processing;
				processing=0;
			}
			code[j++] = current_state;
			break;

		case 6:
			if(source.test(i+processing-2)){
				current_state += 3;
			}else{
				current_state += (source.test(i+processing-1)?1:2);
			}
			i+=processing;
			processing=0;
			code[j++] = current_state;
			break;
		}
	}
	return code;
}

template<std::size_t S>
auto modifiedVLRLL<0x27>::decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
	std::bitset<S*2> decode;
	nucleotide_t previous = initial_state;

	for(std::size_t i=0, j=0; i<S; ++i){
		switch(source[i]-previous){
		case 0:
			decode.set(j++);
			decode.set(j++);
			break;

		case 1:
			decode.reset(j++);
			decode.set(j++);
			break;

		case 2:
			decode.reset(j++);
			decode.reset(j++);
			break;

		case 3:
			decode.set(j++);
			decode.reset(j++);
			break;
		}
		previous = source[i];
	}
	return decode;
}

template<std::floating_point T, std::size_t S>
auto modifiedVLRLL<0x27>::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &source, nucleotide_p<ATGC,T> initial_state){
	constexpr auto p3_to_11 = static_cast<T>(1.0/17.0);
	constexpr auto p3_to_10 = static_cast<T>(16.0/17.0);
	std::array<T,S*2> LLR;
	auto previous = initial_state;

	for(std::size_t i=0; const auto &current: source){
		T PX0 = 0, PX1 = 0, P0X = 0, P1X = 0;
		for(auto j=0ui8; j<4; ++j) P0X += previous[j] * (current[j+1] + current[j+2]); //上位ビットが0: 遷移語が1 or 2になる組み合わせ
		for(auto j=0ui8; j<4; ++j) P1X += previous[j] * (current[j] + current[j+3]); //上位ビットが1: 遷移語が0 or 3になる組み合わせ
		for(auto j=0ui8; j<4; ++j) PX0 += previous[j] * (current[j+2] + p3_to_10*current[j+3]); //下位ビットが0: 遷移語が2 or 3(*p3_to_10)になる組み合わせ
		for(auto j=0ui8; j<4; ++j) PX1 += previous[j] * (current[j] + current[j+1] + p3_to_11*current[j+3]); //下位ビットが1: 遷移語が0 or 1 or 3(*p3_to_11)になる組み合わせ
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

}

#endif