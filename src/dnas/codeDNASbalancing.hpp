#ifndef INCLUDE_GUARD_dnas_codeDNASbalancing
#define INCLUDE_GUARD_dnas_codeDNASbalancing

#include <array>
#include <bitset>
#include <cmath>
#include <cassert>
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
	template<std::size_t BS, std::uint8_t FLAG>
	constexpr auto DivisionBalancingDistribution();
	template<std::size_t BS, std::uint8_t FLAG>
	constexpr auto DivisionBalancingDistribution(double tolerance);
}

template<std::size_t BS>
class DivisionBalancing<0x1B,BS,0> {
	static constexpr std::uint8_t ATGC = 0x1B;
	static constexpr std::array<float,BS+1> dist = DivisionBalancingDistribution<BS,0>();
public:
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source);
	template<std::floating_point T, std::size_t S>
	static auto restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source, T prevdist=0);
};

template<std::size_t BS>
class DivisionBalancing<0x27,BS,1> {
	static constexpr std::uint8_t ATGC = 0x27;
	static constexpr std::array<float,BS+1> dist = DivisionBalancingDistribution<BS,0>();//別の分布が必要
public:
	template<std::size_t S>
	static auto balance(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> prevchange = 0);
	template<std::floating_point T, std::size_t S>
	static auto restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source, T prevdist=0);
};

template<std::size_t BS>
class DivisionBalancing<0x1B,BS,2> {
	static constexpr std::uint8_t ATGC = 0x1B;
	const double tolerance;
	const std::array<float,BS+1> dist;
public:
	DivisionBalancing(double tolerance):tolerance(tolerance),dist(DivisionBalancingDistribution<BS,2>(tolerance)){}
	template<std::size_t S>
	auto balance(const std::array<nucleotide_t<ATGC>,S> &source) const;
	template<std::floating_point T, std::size_t S>
	auto restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source, T prevdist=0) const;
};

template<std::size_t BS>
class DivisionBalancing<0x1B,BS,6> {
	static constexpr std::uint8_t ATGC = 0x1B;
	const double tolerance;
	const std::array<float,BS+1> dist;
public:
	DivisionBalancing(double tolerance):tolerance(tolerance),dist(DivisionBalancingDistribution<BS,6>(tolerance)){}
	template<std::size_t S>
	auto balance(const std::array<nucleotide_t<ATGC>,S> &source) const;
	template<std::floating_point T, std::size_t S>
	auto restore_p(const std::array<nucleotide_p<ATGC,T>,S> &source, T prevdist=0) const;
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
		std::uint64_t qtyATblock=0, qtyATdiv, qtyATdivwide;

		for(std::size_t j=block.head, jend=block.head+divsize; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyATdiv = qtyATblock<<1;
		for(std::size_t j=block.head+divsize, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();

		section div(block.head<<1, blocksize);
		qtyATdivwide = qtyATdiv;
		while(qtyATdivwide!=qtyATblock){
			for(std::size_t k=0, kend=qtyATdiv>qtyATblock?qtyATdiv-qtyATblock:qtyATblock-qtyATdiv; k<kend; ++k){
				qtyATdiv += source[(div.tail++)>>1].is_AT();
				qtyATdiv -= source[(div.head++)>>1].is_AT();
			}
			qtyATdivwide = qtyATdiv + (source[div.tail>>1].is_AT()+source[div.head>>1].is_AT()-1)*(div.head&1);
		}
		for(std::size_t j=div.head>>1, jend=(div.tail+1)>>1; j<jend; ++j) balanced[j]+=2;
	}
	return balanced;
}

template<std::size_t BS>
template<std::floating_point T, std::size_t S>
auto DivisionBalancing<0x1B,BS,0>::restore_p(const std::array<nucleotide_p<ATGC,T>,S> &differential, T prevdist){
	std::array<nucleotide_p<ATGC,T>,S> result = differential;
	if constexpr(BS!=0){
		for(std::size_t i=0; i<S; ++i){
			const auto &di = differential[i];
			auto &ri = result[i];
			auto j = i%BS;
			T prob = dist[j];
			if(j==0){
				prob -= prevdist;
				prevdist = dist[BS];
			}
			for(auto k=0ui8; k<4; ++k) ri[k] += (di[k+2]-di[k])*prob;
		}
	}else{
		prevdist = 0;
	}
	return std::make_pair(result, prevdist);
}

template<std::size_t BS>
template<std::size_t S>
auto DivisionBalancing<0x27,BS,1>::balance(const std::array<nucleotide_t<ATGC>,S> &source, nucleotide_t<ATGC> prevchange){
	constexpr std::size_t blocksize = BS==0?S:BS;
	constexpr std::size_t divsize = blocksize>>1;
	static_assert(blocksize%2==0&&S%blocksize==0);
	std::array<nucleotide_t<ATGC>,S> balanced;
	auto change = prevchange;

	for(std::size_t i=0; i<S/blocksize; ++i){
		section block(i*blocksize, blocksize);
		std::uint64_t qtyATblock=0, qtyATdiv, qtyATdivwide;

		for(std::size_t j=block.head, jend=block.head+divsize; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyATdiv = qtyATblock<<1;
		for(std::size_t j=block.head+divsize, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();

		section div(block.head<<1, blocksize);
		qtyATdivwide = qtyATdiv;
		while(qtyATdivwide!=qtyATblock){
			for(std::size_t k=0, kend=qtyATdiv>qtyATblock?qtyATdiv-qtyATblock:qtyATblock-qtyATdiv; k<kend; ++k){
				qtyATdiv += source[(div.tail++)>>1].is_AT();
				qtyATdiv -= source[(div.head++)>>1].is_AT();
			}
			qtyATdivwide = qtyATdiv + (source[div.tail>>1].is_AT()+source[div.head>>1].is_AT()-1)*(div.head&1);
		}
		//適用
		div.head = div.head>>1;
		div.tail = (div.tail+1)>>1;
		for(std::size_t j=block.head, jend=div.head; j<jend; ++j) balanced[j] = source[j]+change;
		if(div.head!=0&&source[div.head-1]==source[div.head]+1) change += 3;//分割位置の前後が連続しているかどうか調べる
		else change += 1;
		for(std::size_t j=div.head, jend=div.tail; j<jend; ++j) balanced[j] = source[j]+change;
		if(div.tail<source.size()&&source[div.tail-1]==source[div.tail]+3) change += 1;
		else change += 3;
		for(std::size_t j=div.tail; j<block.tail; ++j) balanced[j] = source[j]+change;
	}
	return std::make_pair(balanced, change);
}

template<std::size_t BS>
template<std::floating_point T, std::size_t S>
auto DivisionBalancing<0x27,BS,1>::restore_p(const std::array<nucleotide_p<ATGC,T>,S> &differential, T prevdist){
	constexpr T b1p = static_cast<T>(3.0/4.0);
	constexpr T b3p = 1-b1p;
	std::array<nucleotide_p<ATGC,T>,S> result = differential;
	if constexpr(BS!=0){
		for(std::size_t i=0; i<S; ++i){
			const auto &di = differential[i];
			auto &ri = result[i];
			auto j = i%BS;
			T prob = dist[j];
			if(j==0){
				prob -= prevdist;
				prevdist = dist[BS];
			}
			for(auto k=0ui8; k<4; ++k) ri[k] += (di[k+1]*b1p+di[k+3]*b3p-di[k]*(b1p+b3p))*prob;
		}
	}else{
		prevdist = 0;
	}
	return std::make_pair(result, prevdist);
}

template<std::size_t BS>
template<std::size_t S>
auto DivisionBalancing<0x1B,BS,2>::balance(const std::array<nucleotide_t<ATGC>,S> &source) const{
	constexpr std::size_t blocksize = BS==0?S:BS;
	constexpr std::size_t divsize = blocksize>>1;
	static_assert(blocksize%2==0&&S%blocksize==0);
	const std::uint64_t qtytolerance = static_cast<std::uint64_t>(static_cast<double>(blocksize)*tolerance);
	const std::uint64_t qtyhalflow = divsize-qtytolerance;
	const std::uint64_t qtyhalfhigh = divsize+qtytolerance;
	auto balanced = source;

	for(std::size_t i=0, iend=S/blocksize; i<iend; ++i){
		section block(i*blocksize, blocksize);
		std::uint64_t qtyATblock=0, qtyATdiv, qtyATdivwide, stride;

		for(std::size_t j=block.head, jend=block.head+divsize; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyATdiv = qtyATblock<<1;
		for(std::size_t j=block.head+divsize, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();

		if(qtyATblock>qtyhalfhigh||qtyATblock<qtyhalflow){//許容範囲に収まってたらスキップ
			section div(block.head<<1, blocksize);
			stride = qtyATdiv>qtyATblock?qtyATdiv-qtyATblock:qtyATblock-qtyATdiv;
			while(stride>qtytolerance){
				for(std::size_t k=0, kend=qtyATdiv>qtyATblock?qtyATdiv-qtyATblock:qtyATblock-qtyATdiv; k<kend; ++k){
					qtyATdiv += source[(div.tail++)>>1].is_AT();
					qtyATdiv -= source[(div.head++)>>1].is_AT();
				}
				qtyATdivwide = qtyATdiv + (source[div.tail>>1].is_AT()+source[div.head>>1].is_AT()-1)*(div.head&1);
				stride = qtyATdivwide>qtyATblock?qtyATdivwide-qtyATblock:qtyATblock-qtyATdivwide;
			}
			for(std::size_t j=div.head>>1, jend=(div.tail+1)>>1; j<jend; ++j) balanced[j]+=2;
		}
	}
	return balanced;
}

template<std::size_t BS>
template<std::floating_point T, std::size_t S>
auto DivisionBalancing<0x1B,BS,2>::restore_p(const std::array<nucleotide_p<ATGC,T>,S> &differential, T prevdist) const{
	std::array<nucleotide_p<ATGC,T>,S> result = differential;
	if constexpr(BS!=0){
		for(std::size_t i=0; i<S; ++i){
			const auto &di = differential[i];
			auto &ri = result[i];
			auto j = i%BS;
			T prob = dist[j];
			if(j==0){
				prob -= prevdist;
				prevdist = dist[BS];
			}
			for(auto k=0ui8; k<4; ++k) ri[k] += (di[k+2]-di[k])*prob;
		}
	}else{
		prevdist = 0;
	}
	return std::make_pair(result, prevdist);
}

template<std::size_t BS>
template<std::size_t S>
auto DivisionBalancing<0x1B,BS,6>::balance(const std::array<nucleotide_t<ATGC>,S> &source) const{
	constexpr std::size_t blocksize = BS==0?S:BS;
	constexpr std::size_t divsize = blocksize>>1;
	static_assert(blocksize%2==0&&S%blocksize==0);
	const std::uint64_t qtytolerance = static_cast<std::uint64_t>(static_cast<double>(blocksize)*tolerance);
	const std::uint64_t qtyhalflow = divsize-qtytolerance;
	const std::uint64_t qtyhalfhigh = divsize+qtytolerance;
	const std::uint64_t pitch = qtytolerance*2+1;
	auto balanced = source;

	for(std::size_t i=0, iend=S/blocksize; i<iend; ++i){
		section block(i*blocksize, blocksize);
		std::uint64_t qtyATblock=0, qtyATdiv, qtyATdivwide, stride;

		for(std::size_t j=block.head, jend=block.head+divsize; j<jend; ++j) qtyATblock += source[j].is_AT();
		qtyATdiv = qtyATblock<<1;
		for(std::size_t j=block.head+divsize, jend=block.tail; j<jend; ++j) qtyATblock += source[j].is_AT();

		if(qtyATblock>qtyhalfhigh||qtyATblock<qtyhalflow){//許容範囲に収まってたらスキップ
			section div(block.head<<1, blocksize);
			stride = qtyATdiv>qtyATblock?qtyATdiv-qtyATblock:qtyATblock-qtyATdiv;
			while(stride>qtytolerance){
				for(std::size_t k=0; k<pitch; ++k){
					qtyATdiv += source[(div.tail++)>>1].is_AT();
					qtyATdiv -= source[(div.head++)>>1].is_AT();
				}
				qtyATdivwide = qtyATdiv + (source[div.tail>>1].is_AT()+source[div.head>>1].is_AT()-1)*(div.head&1);
				stride = qtyATdivwide>qtyATblock?qtyATdivwide-qtyATblock:qtyATblock-qtyATdivwide;
			}
			for(std::size_t j=div.head>>1, jend=(div.tail+1)>>1; j<jend; ++j) balanced[j]+=2;
		}

		std::uint64_t qtyATblock2=0;
		for(std::size_t j=block.head, jend=block.tail; j<jend; ++j) qtyATblock2 += balanced[j].is_AT();
		assert(qtyATblock2<=qtyhalfhigh&&qtyATblock2>=qtyhalflow);
	}
	return balanced;
}

template<std::size_t BS>
template<std::floating_point T, std::size_t S>
auto DivisionBalancing<0x1B,BS,6>::restore_p(const std::array<nucleotide_p<ATGC,T>,S> &differential, T prevdist) const{
	std::array<nucleotide_p<ATGC,T>,S> result = differential;
	if constexpr(BS!=0){
		for(std::size_t i=0; i<S; ++i){
			const auto &di = differential[i];
			auto &ri = result[i];
			auto j = i%BS;
			T prob = dist[j];
			if(j==0){
				prob -= prevdist;
				prevdist = dist[BS];
			}
			for(auto k=0ui8; k<4; ++k) ri[k] += (di[k+2]-di[k])*prob;
		}
	}else{
		prevdist = 0;
	}
	return std::make_pair(result, prevdist);
}

////////////////////////////////////////////////////////////////
//                                                            //
//               DivisionBalancingDistribution                //
//                                                            //
////////////////////////////////////////////////////////////////
namespace {

template<std::size_t BS, std::uint8_t FLAG>
constexpr auto DivisionBalancingDistribution(){
	return std::array<float,BS+1>{};
}

template<>
constexpr auto DivisionBalancingDistribution<16,0>(){
	constexpr float sum = 65536.f;
	return std::array{
		25740.f/sum,
		12012.f/sum,
		8316.f/sum,
		6300.f/sum,
		4900.f/sum,
		3780.f/sum,
		2772.f/sum,
		1716.f/sum,
		12870.f/sum,
		18876.f/sum,
		10164.f/sum,
		7308.f/sum,
		5600.f/sum,
		4340.f/sum,
		3276.f/sum,
		2244.f/sum,
		858.f/sum
	};
}

template<>
constexpr auto DivisionBalancingDistribution<8,0>(){
	constexpr float sum = 256.f;
	return std::array{
		140.f/sum,
		60.f/sum,
		36.f/sum,
		20.f/sum,
		70.f/sum,
		100.f/sum,
		48.f/sum,
		28.f/sum,
		10.f/sum
	};
}

template<std::size_t BS, std::uint8_t FLAG>
constexpr auto DivisionBalancingDistribution(double tolerance){
	return std::array<float,BS+1>{};
}

template<>
constexpr auto DivisionBalancingDistribution<16,2>(double tolerance){
	constexpr float sum = 65536.f;
	if(tolerance == 0) return DivisionBalancingDistribution<16,0>();
	if(tolerance == 0.125){
		return std::array{
			12496.f/sum,
			630.f/sum,
			350.f/sum,
			182.f/sum,
			84.f/sum,
			28.f/sum,
			0.f/sum,
			0.f/sum,
			10802.f/sum,
			1834.f/sum,
			570.f/sum,
			310.f/sum,
			158.f/sum,
			72.f/sum,
			24.f/sum,
			0.f/sum,
			0.f/sum
		};
	}else return std::array<float,17>{};
}

template<>
constexpr auto DivisionBalancingDistribution<8,2>(double tolerance){
	constexpr float sum = 256.f;
	if(tolerance == 0) return DivisionBalancingDistribution<8,2>();
	if(tolerance == 0.125){
		return std::array{
			68.f/sum,
			4.f/sum,
			2.f/sum,
			0.f/sum,
			50.f/sum,
			18.f/sum,
			4.f/sum,
			2.f/sum,
			0.f/sum
		};
	}else return std::array<float,9>{};
}

template<>
constexpr auto DivisionBalancingDistribution<16,6>(double tolerance){
	constexpr float sum = 65536.f;
	if(tolerance == 0) return DivisionBalancingDistribution<16,0>();
	if(tolerance == 0.125){
		return std::array{
			10802.f/sum,
			0.f/sum,
			2514.f/sum,
			0.f/sum,
			0.f/sum,
			402.f/sum,
			0.f/sum,
			52.f/sum,
			10802.f/sum,
			0.f/sum,
			0.f/sum,
			2514.f/sum,
			0.f/sum,
			402.f/sum,
			0.f/sum,
			0.f/sum,
			52.f/sum
		};
	}else return std::array<float,17>{};
}

template<>
constexpr auto DivisionBalancingDistribution<8,6>(double tolerance){
	constexpr float sum = 256.f;
	if(tolerance == 0) return DivisionBalancingDistribution<8,2>();
	if(tolerance == 0.125){
		return std::array{
			50.f/sum,
			22.f/sum,
			0.f/sum,
			2.f/sum,
			50.f/sum,
			0.f/sum,
			22.f/sum,
			2.f/sum,
			0.f/sum
		};
	}else return std::array<float,9>{};
}


}

}

#endif