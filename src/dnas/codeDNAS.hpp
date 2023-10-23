#ifndef __code_CODEDNAS__
#define __code_CODEDNAS__

#include <array>
#include <bitset>
// #include <exception>
#include "DNASnttype.hpp"

namespace code::DNAS {

namespace differential {
	template<std::uint8_t ATGC, std::size_t S>
	auto encode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::uint8_t ATGC, std::size_t S>
	auto decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state = 0);
}
namespace VLRLL {
	template<std::uint8_t ATGC=0x1B, std::size_t S>
	auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::uint8_t ATGC, std::size_t S>
	auto decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state = 0);
}
namespace modified_VLRLL {
	template<std::uint8_t ATGC, std::size_t S>
	auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::uint8_t ATGC, std::size_t S>
	auto decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state);
}
class flip_balancing {//ATGC=0x1B
public:
	template<std::uint8_t ATGC, std::size_t R>
	auto balance(const std::array<nucleotide_t<ATGC>,R> &cr, std::size_t qty_AT=0, std::size_t qty_GC=0);
	template<std::uint8_t ATGC, std::size_t R>
	auto restore(const std::array<nucleotide_t<ATGC>,R> &crbar, const std::bitset<R> &flipinfo);
};
template<std::size_t BS=0>//ATGC=0x37
class division_balancing {
public:
	template<std::uint8_t ATGC, std::size_t S>
	auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
};
template<std::uint8_t ATGC, std::size_t S>
auto quarternary_to_binary(const std::array<nucleotide_t<ATGC>,S> &source);
template<std::uint8_t ATGC=0x1B, std::size_t S>
auto binary_to_quarternary(const std::bitset<S> &source);
template<std::uint8_t ATGC, std::size_t L>
auto count_AT(const std::array<nucleotide_t<ATGC>,L> &c, std::uint64_t qty_AT_init=0);


namespace differential {

template<std::uint8_t ATGC, std::size_t S>
auto encode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
	std::array<nucleotide_t<ATGC>,S> code;
	nucleotide_t current_state = initial_state;

	for(std::size_t i=0u, iend=S; i<iend; ++i){
		current_state += source[i];
		code[i] = current_state;
	}
	return code;
}

template<std::uint8_t ATGC, std::size_t S>
auto decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
	std::array<nucleotide_t<ATGC>,S> decode;
	nucleotide_t prev = initial_state;

	for(std::size_t i=0u, iend=S; i<iend; ++i){
		decode[i] = source[i]-prev;
		prev = source[i];
	}
	return decode;
}

}

namespace VLRLL {

template<std::uint8_t ATGC, std::size_t S>
auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code;
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

template<std::uint8_t ATGC, std::size_t S>
auto decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
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

template<std::uint8_t ATGC, std::size_t S>
auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code; 
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

template<std::uint8_t ATGC, std::size_t S>
auto decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
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

template<std::uint8_t ATGC, std::size_t S>
auto flip_balancing::balance(const std::array<nucleotide_t<ATGC>,S> &source, std::size_t qty_AT, std::size_t qty_GC){
	static_assert(ATGC==0x1B);
	std::array<nucleotide_t<ATGC>,S> balanced;
	std::bitset<S> flipinfo; 
	nucleotide_t diff = 0;

	for(uint64_t i=0; i<S; ++i){
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

template<std::uint8_t ATGC, std::size_t S>
auto flip_balancing::restore(const std::array<nucleotide_t<ATGC>,S> &source, const std::bitset<S> &flipinfo){
	static_assert(ATGC==0x1B);
	std::array<nucleotide_t<ATGC>,S> restored;
	nucleotide_t diff = 0;

	for(uint64_t i=0; i<S; i++){
		restored[i] = (flipinfo.test(i)?(source[i]^3):source[i])-diff;
		diff = source[i]-restored[i];
	}
	return restored;
}

template<std::size_t BS>
template<std::uint8_t ATGC, std::size_t S>
auto division_balancing<BS>::balance(const std::array<nucleotide_t<ATGC>,S> &source){
	constexpr std::size_t block_size = BS==0?S:BS;
	constexpr std::size_t div_size = block_size>>1;
	static_assert(block_size%2==0&&S%block_size==0);
	auto balanced = source;

	for(std::size_t i=0u, iend=S/block_size; i<iend; ++i){
		std::size_t block_head=i*block_size, block_tail=block_head+block_size;
		std::uint64_t qtyATblock=0, qtyATdiv, qtyAThalf;

		for(std::size_t j=block_head, jend=block_head+div_size; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyATdiv = qtyATblock;
		for(std::size_t j=block_head+div_size, jend=block_tail; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyAThalf = qtyATblock>>1;

		std::size_t div_head=block_head, div_tail=div_head+div_size;
		if(qtyATblock&1){
			qtyATdiv += source[div_tail++].is_AT();
			++qtyAThalf;
		}
		while(qtyATdiv!=qtyAThalf&&div_tail<block_tail){
			qtyATdiv -= source[div_head++].is_AT();
			qtyATdiv += source[div_tail++].is_AT();
		}
		for(std::size_t j=div_head, jend=div_tail; j<jend; ++j) balanced[j]+=(ATGC==0x1B?2:ATGC==0x27?1:0);
	}
	return balanced;
}

template<std::uint8_t ATGC, std::size_t S>
auto quarternary_to_binary(const std::array<nucleotide_t<ATGC>,S> &source){
	std::bitset<S*2> code;
	for(std::size_t i=0u; const auto &j: source){
		code[i++]=j.msb();
		code[i++]=j.lsb();
	}
	return code;
}

template<std::uint8_t ATGC, std::size_t S>
auto binary_to_quarternary(const std::bitset<S> &source){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code;
	for(std::size_t i=0u; auto &j: code){
		j = static_cast<int>(source.test(i++))<<1;
		j += static_cast<int>(source.test(i++));
	}
	return code;
}

template<std::uint8_t ATGC, std::size_t L>
auto count_AT(const std::array<nucleotide_t<ATGC>,L> &c, std::uint64_t qty_AT_init){
	auto qty_AT=qty_AT_init;
	for(const auto &i: c) qty_AT+=i.is_AT();
	return qty_AT;
}

}

#endif