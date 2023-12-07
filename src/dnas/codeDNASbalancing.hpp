#ifndef INCLUDE_GUARD_dnas_codeDNASbalancing
#define INCLUDE_GUARD_dnas_codeDNASbalancing

#include <array>
#include <bitset>
#include <cmath>
// #include <limits>
#include <stdexcept>
#include "DNASnttype.hpp"

namespace code::DNAS {

struct FlipBalancing {//ATGC=0x1B
	template<std::uint8_t ATGC, std::size_t R>
	static auto balance(const std::array<nucleotide_t<ATGC>,R> &cr, std::size_t qty_AT=0, std::size_t qty_GC=0);
	template<std::uint8_t ATGC, std::size_t R>
	static auto restore(const std::array<nucleotide_t<ATGC>,R> &crbar, const std::bitset<R> &flipinfo);
};

template<std::uint8_t ATGC, std::size_t BS, std::uint8_t FLAG> struct DivisionBalancing;//BS:ブロック長, FLAG:識別子[1:runlength 2:lesschange 4:pitch]

namespace {
	template<std::size_t BS>
	constexpr auto DivisionBalancingDistribution();
	template<std::size_t BS>
	constexpr auto DivisionBalancingDistribution(double tolerance);
}

template<std::size_t BS>
class DivisionBalancing<0x1B,BS,0> {
	static constexpr std::uint8_t ATGC = 0x1B;
	static constexpr std::array<float,BS> dist = DivisionBalancingDistribution<BS>();
public:
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
	template<std::floating_point T, std::size_t S>
	static auto restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source);
};

template<std::size_t BS>
class DivisionBalancing<0x27,BS,1> {
	static constexpr std::uint8_t ATGC = 0x27;
	static constexpr std::array<float,BS> dist = DivisionBalancingDistribution<BS>();
public:
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> diff_init = 0);
	template<std::floating_point T, std::size_t S>
	static auto restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source);
};

template<std::size_t BS>
class DivisionBalancing<0x1B,BS,2> {
	static constexpr std::uint8_t ATGC = 0x1B;
	const double tolerance;
	const std::array<float,BS> dist;
public:
	DivisionBalancing(double tolerance):tolerance(tolerance),dist(DivisionBalancingDistribution<BS>(tolerance)){}
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
	template<std::floating_point T, std::size_t S>
	static auto restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source);
};

template<std::size_t BS>
class DivisionBalancing<0x1B,BS,6> {
	static constexpr std::uint8_t ATGC = 0x1B;
	const double tolerance;
	const std::array<float,BS> dist;
public:
	// DivisionBalancing(double tolerance):tolerance(tolerance),dist(DivisionBalancingDistribution<BS>(tolerance)){}
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
	template<std::floating_point T, std::size_t S>
	static auto restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source);
};

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
	constexpr std::size_t blocksize = BS==0?S:BS;
	constexpr std::size_t divsize = blocksize>>1;
	static_assert(blocksize%2==0&&S%blocksize==0);
	auto balanced = source;

	for(std::size_t i=0, iend=S/blocksize; i<iend; ++i){
		section block(i*blocksize, blocksize);
		std::uint64_t qtyATblock=0, qtyATdiv, stride;

		for(std::size_t j=block.head, jend=block.head+divsize; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyATdiv = qtyATblock<<1;
		for(std::size_t j=block.head+divsize, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();

		section div(block.head<<1, blocksize);
		stride = qtyATdiv>qtyATblock?qtyATdiv-qtyATblock:qtyATblock-qtyATdiv;
		while(stride>0){
			for(std::size_t k=0; k<stride; ++k){
				qtyATdiv += source[(div.tail++)>>1].is_AT();
				qtyATdiv -= source[(div.head++)>>1].is_AT();
			}
			stride = qtyATdiv>qtyATblock?qtyATdiv-qtyATblock:qtyATblock-qtyATdiv;
		}
		for(std::size_t j=div.head>>1, jend=(div.tail+1)>>1; j<jend; ++j) balanced[j]+=2;
	}
	return balanced;
}

template<std::size_t BS>
template<std::floating_point T, std::size_t S>
auto DivisionBalancing<0x1B,BS,0>::restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source){
	constexpr std::size_t blocksize = BS==0?S:BS;
	std::array<nucleotide_p<ATGC,T>,S> result = source;
	for(std::size_t i=0, iend=S/blocksize; i<iend; ++i){
		for(std::size_t j=0; j<blocksize; ++j){
			const auto ij = i*blocksize+j;
			const auto &sj = source[ij];
			auto &rj = result[ij];
			for(auto k=0ui8; k<4; ++k) rj[k] += (sj[k+2]-sj[k])*dist[j]*0.08f;
		}
	}
	return result;
}

template<std::size_t BS>
template<std::size_t S>
auto DivisionBalancing<0x27,BS,1>::balance(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> diff_init){
	constexpr std::size_t block_size = BS==0?S:BS;
	constexpr std::size_t div_size = block_size>>1;
	static_assert(block_size%2==0&&S%block_size==0);
	auto balanced = source;
	auto diff = diff_init;

	for(std::size_t i=0; i<S/block_size; ++i){
		section block(i*block_size, block_size);
		std::uint64_t qtyATblock=0, qtyATdiv, qtyAThalf;

		for(std::size_t j=block.head, jend=block.head+div_size; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyATdiv = qtyATblock;
		for(std::size_t j=block.head+div_size, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyAThalf = qtyATblock>>1;

		section div(block.head, div_size);
		std::pair<nucleotide_t<ATGC>,nucleotide_t<ATGC>> change = {1,(diff==0?3:1)};
		//奇数の場合の処理
		if(qtyATblock&1){
			qtyATdiv += source[div.tail++].is_AT();
			++qtyAThalf;
		}
		//分割区間を探す
		while(div.tail<block.tail && qtyATdiv!=qtyAThalf){
			qtyATdiv += source[div.tail++].is_AT();
			qtyATdiv -= source[div.head++].is_AT();
		}
		//分割位置の連長を調べる
		if(div.head!=0&&source[div.head-1]==source[div.head]+1){
			section run(div.head-1, 2);
			while(run.head!=0&&source[run.head-1]==source[run.head]) --run.head;
			while(run.tail!=source.size()&&source[run.tail-1]==source[run.tail]) ++run.tail;
			if(run.size()>3){
				change.first += 2;
				change.second += 2;
			}
		}
		if(div.tail!=source.size()&&source[div.tail-1]==source[div.tail]+change.second){
			section run(div.tail-1, 2);
			while(run.head!=0&&source[run.head-1]==source[run.head]) --run.head;
			while(run.tail!=source.size()&&source[run.tail-1]==source[run.tail]) ++run.tail;
			if(run.size()>3) change.second += 2;
		}
		//適用
		for(std::size_t j=block.head; j<div.head; ++j) balanced[j]+=diff;
		for(std::size_t j=div.head; j<div.tail; ++j) balanced[j]+=diff+change.first;
		diff += change.first+change.second;
		for(std::size_t j=div.tail; j<block.tail; ++j) balanced[j]+=diff;
	}
	return balanced;
}

template<std::size_t BS>
template<std::floating_point T, std::size_t S>
auto DivisionBalancing<0x27,BS,1>::restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source){
	constexpr std::size_t block_size = BS==0?S:BS;
	constexpr std::size_t div_size = block_size>>1;
	constexpr T divprob = static_cast<T>(1/static_cast<double>(div_size));
	constexpr T b1p = static_cast<T>(7.0/8.0);
	constexpr T b3p = 1-b1p;
	auto result = source;
	for(std::size_t i=0, iend=S/block_size; i<iend; ++i){
		section block(i*block_size, block_size);
		T secprob = 0;
		for(std::size_t j=block.head, jend=block.head+div_size; j<jend; ++j){
			const auto &sj = source[j];
			auto &rj = result[j];
			secprob += divprob;
			for(auto k=0ui8; k<4; ++k) rj[k] += (sj[k+1]*b1p+sj[k+3]*b3p-sj[k])*secprob;
		}
		secprob += divprob/2;
		for(std::size_t j=block.head+div_size, jend=block.tail; j<jend; ++j){
			const auto &sj = source[j];
			auto &rj = result[j];
			secprob -= divprob;
			for(auto k=0ui8; k<4; ++k) rj[k] += (sj[k+1]*b1p+sj[k+3]*b3p-sj[k])*secprob;
		}
	}
	return result;
}

// template<std::size_t BS>
// template<std::size_t S>
// auto DivisionBalancing<0x1B,BS,2>::balance(const std::array<nucleotide_t<ATGC>,S> &source){
// 	constexpr std::size_t block_size = BS==0?S:BS;
// 	constexpr std::size_t div_size = block_size>>1;
// 	static_assert(block_size%2==0&&S%block_size==0);
// 	if(tolerance>0.5||tolerance<0) throw std::invalid_argument("invalid tolerance value");
// 	const std::uint64_t qtyhalflow = static_cast<std::uint64_t>(std::ceil(static_cast<double>(block_size)*(0.5-tolerance)));
// 	const std::uint64_t qtyhalfhigh = static_cast<std::uint64_t>(std::floor(static_cast<double>(block_size)*(0.5+tolerance)));
// 	auto balanced = source;

// 	for(std::size_t i=0, iend=S/block_size; i<iend; ++i){
// 		section block(i*block_size, block_size);
// 		std::uint64_t qtyATblock=0, qtyATdiv, qtyAThalf; //qtyATlow, qtyAThigh;

// 		for(std::size_t j=block.head, jend=block.head+div_size; j<jend; ++j) qtyATblock += source[j].is_AT();
// 		qtyATdiv = qtyATblock;
// 		for(std::size_t j=block.head+div_size, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();
// 		if(qtyATblock>qtyhalfhigh||qtyATblock<qtyhalflow){//許容範囲に収まってたらスキップ
// 			qtyAThalf = qtyATblock>>1;

// 			section div(block.head, div_size);
// 			if(qtyATblock&1){
// 				qtyATdiv += source[div.tail++].is_AT();
// 				++qtyAThalf;
// 			}
// 			while(div.tail<block.tail && qtyATdiv!=qtyAThalf){
// 				qtyATdiv += source[div.tail++].is_AT();
// 				qtyATdiv -= source[div.head++].is_AT();
// 			}
// 			for(std::size_t j=div.head; j<div.tail; ++j) balanced[j]+=2;
// 		}
// 	}
// 	return balanced;
// }

// template<std::size_t BS>
// template<std::floating_point T, std::size_t S>
// auto DivisionBalancing<0x1B,BS,2>::restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source){
// 	constexpr std::size_t block_size = BS==0?S:BS;
// 	constexpr std::size_t div_size = block_size>>1;
// 	const T divprob = static_cast<T>(modprob/static_cast<double>(div_size));
// 	auto result = source;
// 	for(std::size_t i=0, iend=S/block_size; i<iend; ++i){
// 		section block(i*block_size, block_size);
// 		T secprob = 0;
// 		for(std::size_t j=block.head, jend=block.head+div_size; j<jend; ++j){
// 			const auto &sj = source[j];
// 			auto &rj = result[j];
// 			secprob += divprob;
// 			for(auto k=0ui8; k<4; ++k) rj[k] += (sj[k+2]-sj[k])*secprob;
// 		}
// 		secprob += divprob/2;
// 		for(std::size_t j=block.head+div_size, jend=block.tail; j<jend; ++j){
// 			const auto &sj = source[j];
// 			auto &rj = result[j];
// 			secprob -= divprob;
// 			for(auto k=0ui8; k<4; ++k) rj[k] += (sj[k+2]-sj[k])*secprob;
// 		}
// 	}
// 	return result;
// }

// template<std::size_t BS>
// template<std::size_t S>
// auto DivisionBalancing<0x1B,BS,6>::balance(const std::array<nucleotide_t<ATGC>,S> &source){
// 	constexpr std::size_t block_size = BS==0?S:BS;
// 	constexpr std::size_t div_size = block_size>>1;
// 	static_assert(block_size%2==0&&S%block_size==0);
// 	if(tolerance>0.5||tolerance<0) throw std::invalid_argument("invalid tolerance value");
// 	const std::uint64_t qtytolerance = static_cast<std::uint64_t>(std::floor(static_cast<double>(block_size)*tolerance));
// 	const std::uint64_t qtyhalflow = div_size-qtytolerance;
// 	const std::uint64_t qtyhalfhigh = div_size+qtytolerance;
// 	const std::size_t division_pitch = ((qtytolerance>>1)<<1)+1;//分割位置の候補の間隔
// 	auto balanced = source;

// 	for(std::size_t i=0, iend=S/block_size; i<iend; ++i){
// 		section block(i*block_size, block_size);
// 		const std::uint64_t qtyATtolerance = qtytolerance>>1;
// 		std::uint64_t qtyATblock=0, qtyATdiv, qtyAThalf, qtyATlow, qtyAThigh;

// 		for(std::size_t j=block.head, jend=block.head+div_size; j<jend; ++j) qtyATblock += source[j].is_AT();
// 		qtyATdiv = qtyATblock;
// 		for(std::size_t j=block.head+div_size, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();
// 		if(qtyATblock>qtyhalfhigh||qtyATblock<qtyhalflow){//許容範囲に収まってたらスキップ
// 			qtyAThalf = qtyATblock>>1;
// 			qtyATlow = qtyAThalf>qtyATtolerance?qtyAThalf-qtyATtolerance:0;
// 			qtyAThigh = qtyAThalf+qtyATtolerance<div_size?qtyAThalf+qtyATtolerance:div_size;

// 			section div(block.head, div_size);
// 			if(qtyATblock&1){
// 				qtyATdiv += source[div.tail++].is_AT();
// 				++qtyATlow;
// 				++qtyAThigh;
// 			}
// 			while(qtyATdiv<qtyATlow || qtyATdiv>qtyAThigh){
// 				for(std::size_t k=0; div.tail<block.tail&&k<division_pitch; ++k){
// 					qtyATdiv += source[div.tail++].is_AT();
// 					qtyATdiv -= source[div.head++].is_AT();
// 				}
// 				if(div.tail>=block.tail) throw std::exception();
// 			}
// 			for(std::size_t j=div.head; j<div.tail; ++j) balanced[j]+=2;
// 		}
// 	}
// 	return balanced;
// }

// template<std::size_t BS>
// template<std::floating_point T, std::size_t S>
// auto DivisionBalancing<0x1B,BS,6>::restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source){
// 	constexpr std::size_t block_size = BS==0?S:BS;
// 	constexpr std::size_t div_size = block_size>>1;
// 	const std::size_t division_pitch = ((static_cast<std::uint64_t>(std::floor(static_cast<double>(block_size)*tolerance))>>1)<<1)+1;//分割位置の候補の間隔
// 	const T divprob = static_cast<T>(modprob/static_cast<double>((div_size-1)/division_pitch+1));
// 	const T divprobhalf = divprob/2;
// 	auto result = source;
// 	for(std::size_t i=0, iend=S/block_size; i<iend; ++i){
// 		section block(i*block_size, block_size);
// 		T secprob = 0;
// 		for(std::size_t j=block.head, jend=block.head+div_size; j<jend; ++j){
// 			const auto &sj = source[j];
// 			auto &rj = result[j];
// 			if((j-block.head)%division_pitch==0) secprob += divprob;
// 			for(auto k=0ui8; k<4; ++k) rj[k] += (sj[k+2]-sj[k])*secprob;
// 		}
// 		for(std::size_t j=block.head+div_size, jend=block.tail; j<jend; ++j){
// 			const auto &sj = source[j];
// 			auto &rj = result[j];
// 			if((j-(block.head+div_size))%division_pitch==0) secprob -= divprobhalf;
// 			if((j-(block.head+div_size))%division_pitch==1||division_pitch==1) secprob -= divprobhalf;
// 			for(auto k=0ui8; k<4; ++k) rj[k] += (sj[k+2]-sj[k])*secprob;
// 		}
// 	}
// 	return result;
// }

////////////////////////////////////////////////////////////////
//                                                            //
//               DivisionBalancingDistribution                //
//                                                            //
////////////////////////////////////////////////////////////////
namespace {

template<std::size_t BS>
constexpr auto DivisionBalancingDistribution(){
	return std::array<float,BS>{};
}

template<>
constexpr auto DivisionBalancingDistribution<16>(){
	return std::array{
		19734.f/65536.f,
		31020.f/65536.f,
		39618.f/65536.f,
		46648.f/65536.f,
		52558.f/65536.f,
		57572.f/65536.f,
		61818.f/65536.f,
		65536.f/65536.f,
		52666.f/65536.f,
		39796.f/65536.f,
		30358.f/65536.f,
		22768.f/65536.f,
		16438.f/65536.f,
		11088.f/65536.f,
		6578.f/65536.f,
		2860.f/65536.f
	};
}

template<std::size_t BS>
constexpr auto DivisionBalancingDistribution(double tolerance){
	return std::array<float,BS>{};
}

template<>
constexpr auto DivisionBalancingDistribution<16>(double tolerance){
	if(tolerance == 0) return DivisionBalancingDistribution<16>();
	if(tolerance == 0.125){
		return std::array{
			11866.f/65536.f,
			12686.f/65536.f,
			13206.f/65536.f,
			13526.f/65536.f,
			13702.f/65536.f,
			13770.f/65536.f,
			13770.f/65536.f,
			13770.f/65536.f,
			2968.f/65536.f,
			1764.f/65536.f,
			1004.f/65536.f,
			524.f/65536.f,
			228.f/65536.f,
			64.f/65536.f,
			0.f/65536.f,
			0.f/65536.f
		};
	}
}

}

}

#endif