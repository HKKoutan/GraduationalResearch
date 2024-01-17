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
class VLRLL<0x1B> {
	static constexpr std::uint8_t ATGC = 0x1B;
public:
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state = 0);
};

template<>
class VLRLL<0x27> {
	static constexpr std::uint8_t ATGC = 0x27;
public:
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, nucleotide_t<ATGC> initial_state = 0);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state = 0);
};

template<std::uint8_t ATGC> struct puncturedRLL;

// template<>
// class puncturedRLL<0x1B> {
// 	static constexpr std::uint8_t ATGC = 0x1B;
// public:
// 	template<std::size_t S>
// 	static auto encode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state, std::size_t initial_runlength);
// 	template<std::size_t S>
// 	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state);
// 	template<std::floating_point T, std::size_t S>
// 	static auto decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state);
// };

template<>
class puncturedRLL<0x1B> {
	static constexpr std::uint8_t ATGC = 0x1B;
public:
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, std::size_t initial_zeros);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code);
	template<std::floating_point T, std::size_t S>
	static auto decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code);
};

// template<>
// class puncturedRLL<0x27> {
// 	static constexpr std::uint8_t ATGC = 0x27;
// public:
// 	template<std::size_t S>
// 	static auto encode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state, std::size_t initial_runlength);
// 	template<std::size_t S>
// 	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state);
// 	template<std::floating_point T, std::size_t S>
// 	static auto decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state);
// };

template<>
class puncturedRLL<0x27> {
	static constexpr std::uint8_t ATGC = 0x27;
public:
	template<std::size_t S>
	static auto encode(const std::bitset<S> &source, std::size_t initial_zeros);
	template<std::size_t S>
	static auto decode(const std::array<nucleotide_t<ATGC>,S> &code);
	template<std::floating_point T, std::size_t S>
	static auto decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code);
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
//                      class puncturedRLL                    //
//                                                            //
////////////////////////////////////////////////////////////////

// template<std::size_t S>
// auto puncturedRLL<0x1B>::encode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state, std::size_t initial_runlength){
// 	std::array<nucleotide_t<ATGC>,S> code; 
// 	std::size_t runlength = initial_runlength;
// 	nucleotide_t current_state = initial_state;

// 	for(std::size_t i=0; i<S; ++i){
// 		auto si = source[i]+1;
// 		if(si==0){
// 			if(runlength==2){
// 				current_state += 3;
// 				runlength = 0;
// 			}else{
// 				current_state += si;
// 				runlength++;
// 			}
// 		}else{
// 			current_state += si;
// 			runlength = 0;
// 		}
// 		code[i] = current_state;
// 	}

// 	return code;
// }

// template<std::size_t S>
// auto puncturedRLL<0x1B>::decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state){
// 	std::array<nucleotide_t<ATGC>,S> decode;
// 	nucleotide_t previous = initial_state;

// 	for(std::size_t i=0; i<S; ++i){
// 		decode[i] = code[i]-previous+3;
// 		previous = code[i];
// 	}

// 	return decode;
// }

// template<std::floating_point T, std::size_t S>
// auto puncturedRLL<0x1B>::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state){
// 	constexpr auto p3_to_11 = static_cast<T>(1.0/17.0);
// 	constexpr auto p3_to_10 = static_cast<T>(16.0/17.0);
// 	std::array<nucleotide_p<ATGC,T>,S> decode;
// 	auto previous = initial_state;

// 	for(std::size_t i=0; i<S; ++i){
// 		auto &current = code[i];
// 		auto &di = decode[i];
// 		for(auto j=0ui8; j<4; ++j) di[0] += previous[j] * (current[j+1] + p3_to_11*current[j]);
// 		for(auto j=0ui8; j<4; ++j) di[1] += previous[j] * current[j+2];
// 		for(auto j=0ui8; j<4; ++j) di[2] += previous[j] * current[j+3];
// 		for(auto j=0ui8; j<4; ++j) di[3] += previous[j] * (p3_to_10*current[j]); //下位ビットが1: 遷移語が0 or 2 or 3(*p3_to_11)になる組み合わせ
// 		previous = current;
// 	}

// 	return decode;
// }

template<std::size_t S>
auto puncturedRLL<0x1B>::encode(const std::bitset<S> &source, std::size_t initial_zeros){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code; 
	std::size_t processing = initial_zeros<<1;

	std::size_t i=(~processing)+1/*processingの2の補数*/ ,j=0;
	while(j<S/2){
		processing+=2;
		switch(processing){
		case 2:
		case 4:
			if(source.test(i+processing-2)){
				if(source.test(i+processing-1)){
					code[j++] = 0;
				}else{
					code[j++] = 3;
					i+=processing;
					processing=0;
				}
			}else{
				code[j++] = (source.test(i+processing-1)?2:1);
				i+=processing;
				processing=0;
			}
			break;

		case 6:
			if(source.test(i+4)){
				code[j++] += 3;
				i+=processing;
				processing=0;
			}else{
				code[j++] = (source.test(i+5)?2:1);
				i+=processing;
				processing=0;
			}
			break;
		}
	}

	return code;
}

template<std::size_t S>
auto puncturedRLL<0x1B>::decode(const std::array<nucleotide_t<ATGC>,S> &code){
	std::bitset<2*S> result;
	for(std::size_t i=0; auto j: code){
		j+=3;
		result[i++]=j.msb();
		result[i++]=j.lsb();
	}
	return result;
}

template<std::floating_point T, std::size_t S>
auto puncturedRLL<0x1B>::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code){
	constexpr auto p3_to_11 = static_cast<T>(1.0/17.0);
	constexpr auto p3_to_10 = static_cast<T>(16.0/17.0);
	std::array<T,S*2> LLR;
	for(std::size_t i=0; const auto &j: code){
		auto PX0 = j[1]+p3_to_10*j[3], PX1 = j[0]+j[2]+p3_to_11*j[3], P0X = j[1]+j[2], P1X = j[0]+j[3];
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

// template<std::size_t S>
// auto puncturedRLL<0x27>::encode(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> initial_state, std::size_t initial_runlength){
// 	std::array<nucleotide_t<ATGC>,S> code; 
// 	std::size_t runlength = initial_runlength;
// 	nucleotide_t current_state = initial_state;

// 	for(std::size_t i=0; i<S; ++i){
// 		auto si = source[i]+2;
// 		if(si==0){
// 			if(runlength==2){
// 				current_state += 1;
// 				runlength = 0;
// 			}else{
// 				current_state += si;
// 				runlength++;
// 			}
// 		}else{
// 			current_state += si;
// 			runlength = 0;
// 		}
// 		code[i] = current_state;
// 	}

// 	return code;
// }

// template<std::size_t S>
// auto puncturedRLL<0x27>::decode(const std::array<nucleotide_t<ATGC>,S> &code, nucleotide_t<ATGC> initial_state){
// 	std::array<nucleotide_t<ATGC>,S> decode;
// 	nucleotide_t previous = initial_state;

// 	for(std::size_t i=0; i<S; ++i){
// 		decode[i] = code[i]-previous+2;
// 		previous = code[i];
// 	}

// 	return decode;
// }

// template<std::floating_point T, std::size_t S>
// auto puncturedRLL<0x27>::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code, nucleotide_p<ATGC,T> initial_state){
// 	constexpr auto p3_to_11 = static_cast<T>(1.0/17.0);
// 	constexpr auto p3_to_10 = static_cast<T>(16.0/17.0);
// 	std::array<nucleotide_p<ATGC,T>,S> decode;
// 	auto previous = initial_state;

// 	for(std::size_t i=0; i<S; ++i){
// 		auto &current = code[i];
// 		auto &di = decode[i];
// 		for(auto j=0ui8; j<4; ++j) di[0] += previous[j] * (current[j+2] + p3_to_11*current[j+3]);
// 		for(auto j=0ui8; j<4; ++j) di[1] += previous[j] * (p3_to_10*current[j+3]);
// 		for(auto j=0ui8; j<4; ++j) di[2] += previous[j] * current[j];
// 		for(auto j=0ui8; j<4; ++j) di[3] += previous[j] * current[j+1]; //下位ビットが1: 遷移語が0 or 2 or 3(*p3_to_11)になる組み合わせ
// 		previous = current;
// 	}

// 	return decode;
// }

template<std::size_t S>
auto puncturedRLL<0x27>::encode(const std::bitset<S> &source, std::size_t initial_zeros){
	static_assert(S%2==0);
	std::array<nucleotide_t<ATGC>,S/2> code; 
	std::size_t processing = initial_zeros<<1;

	std::size_t i=(~processing)+1/*processingの2の補数*/ ,j=0;
	while(j<S/2){
		processing+=2;
		switch(processing){
		case 2:
		case 4:
			if(source.test(i+processing-2)){
				if(source.test(i+processing-1)){
					code[j++] = 0;
				}else{
					code[j++] = 3;
					i+=processing;
					processing=0;
				}
			}else{
				code[j++] = (source.test(i+processing-1)?1:2);
				i+=processing;
				processing=0;
			}
			break;

		case 6:
			if(source.test(i+4)){
				code[j++] += 3;
				i+=processing;
				processing=0;
			}else{
				code[j++] = (source.test(i+5)?1:2);
				i+=processing;
				processing=0;
			}
			break;
		}
	}

	return code;
}

template<std::size_t S>
auto puncturedRLL<0x27>::decode(const std::array<nucleotide_t<ATGC>,S> &code){
	std::bitset<2*S> result;
	for(std::size_t i=0; auto j: code){
		j = 2-j;
		auto msb = j.msb();
		result[i++] = msb;
		result[i++] = j.lsb()^msb;
	}
	return result;
}

template<std::floating_point T, std::size_t S>
auto puncturedRLL<0x27>::decode_p(const std::array<nucleotide_p<ATGC,T>,S> &code){
	constexpr auto p3_to_11 = static_cast<T>(1.0/17.0);
	constexpr auto p3_to_10 = static_cast<T>(16.0/17.0);
	std::array<T,S*2> LLR;
	for(std::size_t i=0; const auto &j: code){
		auto PX0 = j[2]+p3_to_10*j[3], PX1 = j[0]+j[1]+p3_to_11*j[3], P0X = j[1]+j[2], P1X = j[0]+j[3];
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

}

#endif