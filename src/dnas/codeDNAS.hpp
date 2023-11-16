#ifndef INCLUDE_GUARD_dnas_codeDNAS
#define INCLUDE_GUARD_dnas_codeDNAS

#include <array>
#include <bitset>
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

template<std::uint8_t ATGC> struct modifiedVLRLL;

template<>
struct modifiedVLRLL<0x1B> {
	static constexpr std::uint8_t ATGC = 0x1B;
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state);
	// template<std::floating_point T, std::size_t S>
	// static auto decode(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_t<ATGC> initial_state);
};

struct FlipBalancing {//ATGC=0x1B
	template<std::uint8_t ATGC, std::size_t R>
	static auto balance(const std::array<nucleotide_t<ATGC>,R> &cr, std::size_t qty_AT=0, std::size_t qty_GC=0);
	template<std::uint8_t ATGC, std::size_t R>
	static auto restore(const std::array<nucleotide_t<ATGC>,R> &crbar, const std::bitset<R> &flipinfo);
};

template<std::uint8_t ATGC, std::size_t BS, std::uint8_t FLAG> struct DivisionBalancing;//BS:ブロック長, FLAG:操作情報を記録するかのビットフラグ LPBxxxxx L:位置 P:偶奇 B:境界

template<std::size_t BS>
struct DivisionBalancing<0x1B,BS,0> {
	static constexpr std::uint8_t ATGC = 0x1B;
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
};

template<std::size_t BS>
struct DivisionBalancing<0x27,BS,0> {
	static constexpr std::uint8_t ATGC = 0x27;
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
};

template<>
struct DivisionBalancing<0x27,8,0x03> {
	static constexpr std::uint8_t ATGC = 0x27;
	template<std::size_t S>
	static auto encode(const std::array<nucleotide_t<ATGC>,S> &source);
};

template<std::uint8_t ATGC>
struct convert {
	template<std::size_t S>
	static auto quarternary_to_binary(const std::array<nucleotide_t<ATGC>,S> &source);
	template<std::size_t S>
	static auto binary_to_quarternary(const std::bitset<S> &source);
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
auto VLRLL<0x1B>::decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
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

////////////////////////////////////////////////////////////////
//                                                            //
//                     class modifiedVLRLL                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S>
auto modifiedVLRLL<0x1B>::encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state){
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

template<std::size_t S>
auto modifiedVLRLL<0x1B>::decode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state){
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
		while(qtyATdiv!=qtyAThalf && div_tail<block_tail){
			qtyATdiv -= source[div_head++].is_AT();
			qtyATdiv += source[div_tail++].is_AT();
		}
		for(std::size_t j=div_head; j<div_tail; ++j) balanced[j]+=2;
	}
	return balanced;
}

template<std::size_t BS>
template<std::size_t S>
auto DivisionBalancing<0x27,BS,0>::balance(const std::array<nucleotide_t<ATGC>,S> &source){
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
		while(qtyATdiv!=qtyAThalf && div_tail<block_tail){
			qtyATdiv -= source[div_head++].is_AT();
			qtyATdiv += source[div_tail++].is_AT();
		}
		for(std::size_t j=div_head; j<div_tail; ++j) balanced[j]+=1;
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
auto convert<ATGC>::quarternary_to_binary(const std::array<nucleotide_t<ATGC>,S> &source){
	std::bitset<S*2> code;
	for(std::size_t i=0u; const auto &j: source){
		code[i++]=j.msb();
		code[i++]=j.lsb();
	}
	return code;
}

template<std::uint8_t ATGC>
template<std::size_t S>
auto convert<ATGC>::binary_to_quarternary(const std::bitset<S> &source){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code;
	for(std::size_t i=0u; auto &j: code){
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