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
			current_state = (current_state + (source.test(i+processing-1)?2:1)) & 3;
			code[j++] = current_state;
			processing=0;
			break;
		}
	}
	//端数処理(0パディング)
	if(processing%2 != 0){
		current_state = (current_state + (source.test(i+processing-1)?3:1)) & 3;
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

		void VLRLL_decode(const std::vector<nucleotide_t> &source, std::vector<bool> &decode, const nucleotide_t initial_state = 0);
		void modified_VLRLL_decode(const std::vector<nucleotide_t> &source, std::vector<bool> &decode, const nucleotide_t initial_state = 0);

		void nt_addequalizing_encode(const std::vector<nucleotide_t> &cr, std::vector<nucleotide_t> &crbar, std::vector<bool> &info, std::uint32_t qty_AT=0, std::uint32_t qty_GC=0);
		void nt_addequalizing_decode(const std::vector<nucleotide_t> &crbar, const std::vector<bool> &i, std::vector<nucleotide_t> &cr);
		void nt_qty_count(const std::vector<nucleotide_t> &c, std::uint32_t &qty_AT, std::uint32_t &qty_GC);

}

#endif