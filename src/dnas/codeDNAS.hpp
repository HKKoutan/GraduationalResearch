#ifndef __code_CODEDNAS__
#define __code_CODEDNAS__

#include <array>
#include <bitset>
// #include <exception>
#include "DNASnttype.hpp"

namespace code::DNAS {

namespace VLRLL {
	template<std::size_t S>
	auto encode(const std::bitset<S> &source, nucleotide_t initial_state = 0);
	template<std::size_t S>
	auto decode(const std::array<nucleotide_t,S> &source, nucleotide_t initial_state = 0);
}
namespace modified_VLRLL {
	template<std::size_t S>
	auto encode(const std::bitset<S> &source, nucleotide_t initial_state = 0);
	template<std::size_t S>
	auto decode(const std::array<nucleotide_t,S> &source, nucleotide_t initial_state);
}
namespace interim_mapping {
	template<std::size_t S>
	auto map(const std::array<nucleotide_t,S> &source);
	template<std::size_t S>
	auto demap(const std::bitset<S> &source);
}
namespace flip_balancing {
	template<std::size_t R, std::uint8_t ATGC=0x1B>
	auto balance(const std::array<nucleotide_t,R> &cr, std::size_t qty_AT=0, std::size_t qty_GC=0);
	template<std::size_t R, std::uint8_t ATGC=0x1B>
	auto restore(const std::array<nucleotide_t,R> &crbar, const std::bitset<R> &flipinfo);
}
template<std::size_t L, std::uint8_t ATGC=0x1B>
auto count_AT(const std::array<nucleotide_t,L> &c, std::size_t qty_AT_init=0);


namespace VLRLL {

template<std::size_t S>
auto encode(const std::bitset<S> &source, nucleotide_t initial_state){
	std::array<nucleotide_t,S/2> code;
	std::size_t processing = 0u;
	nucleotide_t current_state = initial_state;

	std::size_t i=0u ,j=0u;
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
	//端数処理(0パディング)
	if(processing%2 != 0){
		current_state = current_state+(source.test(i+processing-1)?3:1);
		code[j++] = current_state;
	}

	std::bitset<S> used;
	used.set();
	used >>= S-(i+processing);

	return std::make_pair(code,used);
}

template<std::size_t S>
auto decode(const std::array<nucleotide_t,S> &source, nucleotide_t initial_state){
	std::bitset<S*2> decode;
	nucleotide_t previous = initial_state;
	std::size_t zeros = 0u;

	for(std::size_t i=0u, iend=S, j=0u; i<iend; ++i){
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

}

namespace modified_VLRLL {

template<std::size_t S>
auto encode(const std::bitset<S> &source, nucleotide_t initial_state){
	static_assert(S%2==0);
	std::array<nucleotide_t,S/2> code; 
	std::size_t processing = 0u;
	nucleotide_t current_state = initial_state;

	std::size_t i=0u, j=0u;
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
auto decode(const std::array<nucleotide_t,S> &source, nucleotide_t initial_state){
	std::bitset<S*2> decode;
	nucleotide_t previous = initial_state;

	for(std::size_t i=0u, iend=S/2, j=0u; i<iend; ++i){
		switch(source[i]-previous){
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
		previous = source[i];
	}
	return decode;
}

}

namespace interim_mapping {

template<std::size_t S>
auto map(const std::array<nucleotide_t,S> &source){
	std::bitset<S*2> code;
	for(std::size_t i=0u; const auto &j: source){
		code[i]=j.msb();
		code[i+1]=j.lsb();
		i+=2;
	}
	return code;
}

template<std::size_t S>
auto demap(const std::bitset<S> &source){
	static_assert(S%2==0);
	std::array<nucleotide_t,S/2> code;
	for(std::size_t i=0u; auto &j: code){
		j = (static_cast<int>(source.test(i))<<1)+static_cast<int>(source.test(i+1));
		i+=2;
	}
	return code;
}

}

namespace flip_balancing {

template<std::size_t S, std::uint8_t ATGC=0x1B>
auto balance(const std::array<nucleotide_t,S> &source, std::size_t qty_AT, std::size_t qty_GC){
	static_assert(ATGC==0x1B);//ATGC=0x1Bの場合のみ対応
	std::array<nucleotide_t,S> balanced;
	std::bitset<S> flipinfo; 
	nucleotide_t diff = 0;

	for(uint64_t i=0; i<R; ++i){
		nucleotide_t prev_state = source[i]+diff;
		if(qty_AT>qty_GC && !prev_state.msb()){
			balanced[i] = prev_state^3;
			flipinfo.set(i);
		}else if(qty_AT<qty_GC && prev_state.msb()){
			balanced[i] = prev_state^3;
			flipinfo.set(i);
		}else{
			balanced[i] = prev_state;
			flipinfo.reset(i);
		}
		diff = balanced[i]-source[i];
		qty_AT += !balanced[i].msb();
		qty_GC += balanced[i].msb();
	}
	return std::make_pair(balanced,flipinfo);
}

template<std::size_t S, std::uint8_t ATGC=0x1B>
auto restore(const std::array<nucleotide_t,S> &source, const std::bitset<S> &flipinfo){
	static_assert(ATGC==0x1B);//ATGC=0x1Bの場合のみ対応
	std::array<nucleotide_t,S> restored;
	nucleotide_t diff = 0;

	for(uint64_t i=0; i<S; i++){
		restored[i] = (flipinfo.test(i)?(source[i]^3):source[i])-diff;
		diff = source[i]-restored[i];
	}
	return restored;
}

}

template<std::size_t L, std::uint8_t ATGC=0x1B>
auto count_AT(const std::array<nucleotide_t,L> &c, std::size_t qty_AT_init=0){
	auto qty_AT=qty_AT_init;
	constexpr std::uint8_t A = (ATGC>>6)&3;
	constexpr std::uint8_t T = (ATGC>>4)&3;
	constexpr std::uint8_t G = (ATGC>>2)&3;
	constexpr std::uint8_t C = ATGC&3;
	static_assert((A!=T)&&(A!=G)&&(A!=C)&&(T!=G)&&(T!=C)&&(G!=C));
	for(const auto &i: c){
		switch(i){
		case A:
		case T:
			++qty_AT;
			break;

		// case G:
		// case C:
		// 	break;
		}
	}
	return qty_AT;
}

}

#endif