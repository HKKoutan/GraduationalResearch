#ifndef INCLUDE_GUARD_dnas_codeDNAS
#define INCLUDE_GUARD_dnas_codeDNAS

#include <array>
#include <bitset>
#include <cmath>
#include <limits>
#include "DNASnttype.hpp"

namespace code::DNAS {

template<std::uint8_t ATGC> struct differential;

template<>
struct differential<0x1B> {
	static constexpr std::uint8_t ATGC = 0x1B;
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state = 0);
	// template<std::floating_point T, std::size_t S>
	// static auto decode(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_t<ATGC> initial_state = 0);
};

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

struct FlipBalancing {//ATGC=0x1B
	template<std::uint8_t ATGC, std::size_t R>
	static auto balance(const std::array<nucleotide_t<ATGC>,R> &cr, std::size_t qty_AT=0, std::size_t qty_GC=0);
	template<std::uint8_t ATGC, std::size_t R>
	static auto restore(const std::array<nucleotide_t<ATGC>,R> &crbar, const std::bitset<R> &flipinfo);
};

template<std::uint8_t ATGC, std::size_t BS, std::uint8_t FLAG> struct DivisionBalancing;//BS:ブロック長, FLAG:識別子[0:default 1:runlength 2:minchange 3:runlength]

template<std::size_t BS>
struct DivisionBalancing<0x1B,BS,0> {
	static constexpr std::uint8_t ATGC = 0x1B;
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
};

template<std::size_t BS>
struct DivisionBalancing<0x27,BS,1> {
	static constexpr std::uint8_t ATGC = 0x27;
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
};

template<>
struct DivisionBalancing<0x27,6,0x03> {
	static constexpr std::uint8_t ATGC = 0x27;
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
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

template<std::uint8_t ATGC, std::size_t L>
auto countAT(const std::array<nucleotide_t<ATGC>,L> &c);
template<std::uint8_t ATGC, std::size_t L>
auto countRunlength(const std::array<nucleotide_t<ATGC>,L> &c);

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
auto modifiedVLRLL<0x1B>::decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
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

template<std::floating_point T, std::size_t S>
auto modifiedVLRLL<0x1B>::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &source, nucleotide_p<ATGC,T> initial_state){
	constexpr auto p3_to_11 = static_cast<T>(1.0/17.0);
	constexpr auto p3_to_10 = static_cast<T>(16.0/17.0);
	std::array<T,S*2> LLR;
	auto previous = initial_state;

	for(std::size_t i=0; const auto &current: source){
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

////////////////////////////////////////////////////////////////
//                                                            //
//                    class flip_balancing                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint8_t ATGC, std::size_t S>
auto FlipBalancing::balance(const std::array<nucleotide_t<ATGC>,S> &source, std::size_t qty_AT, std::size_t qty_GC){
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
auto FlipBalancing::restore(const std::array<nucleotide_t<ATGC>,S> &source, const std::bitset<S> &flipinfo){
	static_assert(ATGC==0x1B);
	std::array<nucleotide_t<ATGC>,S> restored;
	nucleotide_t diff = 0;

	for(uint64_t i=0; i<S; i++){
		restored[i] = (flipinfo.test(i)?(source[i]^3):source[i])-diff;
		diff = source[i]-restored[i];
	}
	return restored;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                  class DivisionBalancing                   //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t BS>
template<std::size_t S>
auto DivisionBalancing<0x1B,BS,0>::balance(const std::array<nucleotide_t<ATGC>,S> &source){
	constexpr std::size_t block_size = BS==0?S:BS;
	constexpr std::size_t div_size = block_size>>1;
	static_assert(block_size%2==0&&S%block_size==0);
	auto balanced = source;

	for(std::size_t i=0u, iend=S/block_size; i<iend; ++i){
		section block(i*block_size, block_size);
		std::uint64_t qtyATblock=0, qtyATdiv, qtyAThalf;

		for(std::size_t j=block.head, jend=block.head+div_size; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyATdiv = qtyATblock;
		for(std::size_t j=block.head+div_size, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyAThalf = qtyATblock>>1;

		section div(block.head, div_size);
		if(qtyATblock&1){
			qtyATdiv += source[div.tail++].is_AT();
			++qtyAThalf;
		}
		while(qtyATdiv!=qtyAThalf && div.tail<block.tail){
			qtyATdiv -= source[div.head++].is_AT();
			qtyATdiv += source[div.tail++].is_AT();
		}
		for(std::size_t j=div.head; j<div.tail; ++j) balanced[j]+=2;
	}
	return balanced;
}

template<std::size_t BS>
template<std::size_t S>
auto DivisionBalancing<0x27,BS,1>::balance(const std::array<nucleotide_t<ATGC>,S> &source){
	constexpr std::size_t block_size = BS==0?S:BS;
	constexpr std::size_t div_size = block_size>>1;
	static_assert(block_size%2==0&&S%block_size==0);
	auto balanced = source;
	nucleotide_t<ATGC> diff = 0;

	for(std::size_t i=0; i<S/block_size; ++i){
		section block(i*block_size, block_size);
		std::uint64_t qtyATblock=0, qtyATdiv, qtyAThalf;

		for(std::size_t j=block.head, jend=block.head+div_size; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyATdiv = qtyATblock;
		for(std::size_t j=block.head+div_size, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyAThalf = qtyATblock>>1;

		section div(block.head, div_size);
		std::pair<nucleotide_t<ATGC>,nucleotide_t<ATGC>> change = {1,3};
		//奇数の場合の処理
		if(qtyATblock&1){
			qtyATdiv += source[div.tail++].is_AT();
			++qtyAThalf;
		}
		//分割区間を探す
		while(qtyATdiv!=qtyAThalf && div.tail<block.tail){
			qtyATdiv -= source[div.head++].is_AT();
			qtyATdiv += source[div.tail++].is_AT();
		}
		//分割位置の連長を調べる
		if(source[div.head-1]==source[div.head]+1){
			section run(block.head-1, 1);
			while(source[run.head-1]==source[run.head]) --run.head;
			while(source[run.tail]==source[run.tail+1]) ++run.tail;
			if(run.size()>3) change.first = 3;
		}
		if(source[div.tail-1]==source[div.tail]+3){
			section run(block.head-1, 1);
			while(source[run.head-1]==source[run.head]) --run.head;
			while(source[run.tail]==source[run.tail+1]) ++run.tail;
			if(run.size()>3) change.second = 1;
		}
		//適用
		for(std::size_t j=block.head; j<div.head; ++j) balanced[j]+=diff;
		for(std::size_t j=div.head; j<div.tail; ++j) balanced[j]+=diff+change.first;
		diff += change.first+change.second;
		for(std::size_t j=div.tail; j<block.tail; ++j) balanced[j]+=diff;
	}
	return balanced;
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

////////////////////////////////////////////////////////////////
//                                                            //
//                       miscellaneous                        //
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

}

#endif