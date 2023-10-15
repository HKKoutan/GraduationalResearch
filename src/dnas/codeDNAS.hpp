#ifndef __code_CODEDNAS__
#define __code_CODEDNAS__

#include <iostream>
#include <cstdint>
#include <array>
#include <bitset>
#include <vector>
#include <exception>
#include "DNASnttype.hpp"

namespace code::DNAS {

template<std::size_t S>
auto VLRLL_encode(const std::bitset<S> &source, std::array<nucleotide_t,S/2> &code, nucleotide_t initial_state = 0){
	std::size_t processing;
	nucleotide_t current_state = initial_state;

	size_t i,j=0u;
	while(j<S/2){
		++processing;
		switch(processing){
		case 1:
		case 3:
			break;
		
		case 2:
		case 4:
			if(source.test(i+processing-2)){
				current_state += (source.test(i+processing-1)?0:3);
				if(!source.test(i+processing-1)) processing=0;
			}else{
				current_state += (source.test(i+processing-1)?2:1);
				processing=0;
			}
			code[j++] = current_state;
			break;

		case 5:
			if(source.test(i+processing-1)){
				current_state += 3;
				code[j++] = current_state;
				processing=0;
			}
			break;
		
		case 6:
			current_state = current_state+(source.test(i+processing-1)?2:1);
			code[j++] = current_state;
			processing=0;
			break;
		}
	}
	//端数処理(0パディング)
	if(processing%2 != 0){
		current_state = current_state+(source.test(i+processing-1)?3:1);
		code[j++] = current_state;
	}

	return i;
}

template<std::size_t S>
void modified_VLRLL_encode(const std::bitset<S> &source, std::array<nucleotide_t,S/2> &code, nucleotide_t initial_state = 0){
	static_assert(S%2==0);
	std::size_t processing;
	nucleotide_t current_state = initial_state;

	size_t i,j=0u;
	while(j<S/2){
		++processing;
		switch(processing){
		case 0:
		case 1:
		case 3:
		case 5:
			break;
		
		case 2:
		case 4:
			if(source.test(i+processing-2)){
				current_state += (source.test(i+processing-1)?0:3);
				if(!source.test(i+processing-1)) processing=0;
			}else{
				current_state += (source.test(i+processing-1)?2:1);
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
			code[j++] = current_state;
			processing=0;
			break;
		}
	}
}

template<std::size_t C>
void interim_map(const std::array<nucleotide_t,C/2> &source, std::bitset<C> &code){
	static_assert(C%2==0);
	std::size_t i=0u;
	for(const auto &j: source){
		code[i]=j.lsb();
		code[i+1]=j.msb();
		i+=2;
	}
}

template<std::size_t S>
void interim_demap(const std::bitset<S> &source, std::array<nucleotide_t,S/2> &code){
	static_assert(S%2==0);
	std::size_t i=0u;
	for(const auto &j: source){
		j=(source.test(i)?2:0)+static_cast<int>(source.test(i+1));
		i+=2;
	}
}	

template<std::size_t S>
void VLRLL_decode(const std::array<nucleotide_t,S/2> &source, std::bitset<S> &decode, nucleotide_t initial_state = 0){
	nucleotide_t previous = initial_state;
	std::size_t zeros = 0u;

	for(std::size_t i=0u, iend=S/2, j=0u; i<iend; ++i){
		switch(source[i]-previous){
		case 0:
			decode.set(j++);
			decode.set(j++);
			zeros++;
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
}

template<std::size_t S>
void modified_VLRLL_decode(const std::array<nucleotide_t,S/2> &source, std::bitset<S> &decode, nucleotide_t initial_state = 0){
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
			if(zeros<2) decode.reset(j++);
			break;
		}
		previous = source[i];
	}
}

template<std::size_t R, std::uint8_t ATGC=0x1B>
void nt_addequalizing_encode(const std::array<nucleotide_t,R> &cr, std::array<nucleotide_t,R> &crbar, std::bitset<R> &flipinfo, std::size_t qty_AT=0, std::size_t qty_GC=0){
	static_assert(ATGC==0x1B);//ATGC=0x1Bの場合のみ対応
	nucleotide_t diff = 0;

	for(uint64_t i=0; i<R; ++i){
		nucleotide_t prev_state = cr[i]+diff;
		if(qty_AT>qty_GC && !prev_state.msb()){
			crbar[i] = prev_state^3;
			flipinfo.set(i);
		}else if(qty_AT<qty_GC && prev_state.msb()){
			crbar[i] = prev_state^3;
			flipinfo.set(i);
		}else{
			crbar[i] = prev_state;
			flipinfo.reset(i);
		}
		diff = crbar[i]-cr[i];
		qty_AT += !crbar[i].msb();
		qty_GC += crbar[i].msb();
	}
}

template<std::size_t R, std::uint8_t ATGC=0x1B>
void nt_addequalizing_decode(const std::array<nucleotide_t,R> &crbar, const std::bitset<R> &flipinfo, std::array<nucleotide_t,R> &cr){
	static_assert(ATGC==0x1B);//ATGC=0x1Bの場合のみ対応
	nucleotide_t diff = 0;

	for(uint64_t i=0; i<R; i++){
		if(flipinfo.test(i)) cr[i] = (~crbar[i])-diff;
		else cr[i] = crbar[i]-diff;
		diff = crbar[i]-cr[i];
	}
}

template<std::size_t L, std::uint8_t ATGC=0x1B>
void nt_qty_count(const std::array<nucleotide_t,L> &c, std::size_t &qty_AT, std::size_t &qty_GC){
	qty_AT=0;
	qty_GC=0;
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

		case G:
		case C:
			++qty_GC;
			break;
		}
	}
}

}

#endif