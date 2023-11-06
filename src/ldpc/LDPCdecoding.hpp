#ifndef INCLUDE_GUARD_ldpc_LDPCdecoding
#define INCLUDE_GUARD_ldpc_LDPCdecoding

#include <bitset>
#include <cstdint>
#include <bit>
#include <algorithm>
#include <cmath>
#include <limits>
#include <exception>
#include "LDPCCheckMatrix.hpp"
#include "LDPCfuncGallager.hpp"

namespace code::LDPC {

template<class T>
concept DecoderType = requires(std::remove_reference_t<T> &x){
	T::rowupdate();
};

template<CheckMatrix T>
class Iterative_decoding {
	using fptype = float;
	using signtype = std::uint32_t;
	static_assert(sizeof(fptype)==sizeof(signtype));
	static constexpr signtype signmask = 1u<<(sizeof(signtype)*8u-1u);
	static constexpr std::size_t S = T::sourcesize();
	static constexpr std::size_t C = T::codesize();
	static constexpr std::size_t Hsize = C-S;

	inline static T H;//検査行列

	//alphabeta, alphabetapの型と初期化関数をTの型に応じて変える
	template<CheckMatrix U>
	struct decodingData {
		static constexpr std::size_t S = U::sourcesize();
		static constexpr std::size_t C = U::codesize();
		using datatype = std::vector<std::pair<std::array<fptype,C>,std::array<fptype,C>>>;
		using pointertype = std::array<std::vector<std::pair<fptype*,const fptype*>>,C-S>;
		static void init();
	};
	template<std::size_t S, std::size_t C, std::size_t W>
	struct decodingData<CheckMatrix_regular<S,C,W>> {
		using datatype = std::array<std::pair<std::array<fptype,C>,std::array<fptype,C>>,W*S/C>;
		using pointertype = std::array<std::array<std::pair<fptype*,const fptype*>,W>,C-S>;
		static void init();
	};

	inline static thread_local decodingData<T>::datatype alphabeta;
	inline static thread_local decodingData<T>::pointertype alphabetap;

public:
	explicit Iterative_decoding(const T &H);
	void decode_init();//decodeで使用する変数の初期化
	template<DecoderType D, std::floating_point U>
	bool iterate(std::array<U,C> &LPR, const std::array<U,C> &LLR);
	template<std::floating_point U>
	auto estimate(const std::array<U,C> &LEVR) const;//推定符号語を求める
	//rowupdate
	struct SumProduct {
		inline static const func_Gallager_halftable fg;
		static void rowupdate();
	};
	struct MinSum {
		static void rowupdate();
	};
};

////////////////////////////////////////////////////////////////
//                                                            //
//                  class Iterative_decoding                  //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
template<CheckMatrix U>
void Iterative_decoding<T>::decodingData<U>::init(){
	std::array<std::size_t,C> Hheight{};
	for(const auto &Hi: H) for(auto j: Hi) ++Hheight[j];
	alphabeta.resize(std::ranges::max(Hheight));

	std::array<std::vector<std::uint64_t>,C> HT{};//Hの転置
	for(std::size_t i=0; i<Hsize; ++i) for(auto j: H[i]) HT[j].push_back(i);

	for(std::size_t i=0; i<Hsize; ++i){
		auto &Hi = H[i];
		auto &abpi = alphabetap[i];
		//Hとalphabetapの要素数を揃える
		abpi.resize(Hi.size());
		//alpha<-alphap beta<-betap
		for(std::size_t j=0u, jend=Hi.size(); j<jend; ++j){
			auto &hij = Hi[j];
			auto &Hj = HT[hij];
			std::size_t k=0;
			while(Hj[k]!=i) ++k;
			auto &[ai, bi] = alphabeta[k];
			abpi[j] = std::make_pair(&ai[hij],&bi[hij]);
		}
	}
}

template<CheckMatrix T>
template<std::size_t S, std::size_t C, std::size_t W>
void Iterative_decoding<T>::decodingData<CheckMatrix_regular<S,C,W>>::init(){
	std::array<std::vector<std::uint64_t>,C> HT{};//Hの転置
	for(std::size_t i=0; i<Hsize; ++i) for(auto j: H[i]) HT[j].push_back(i);

	for(std::size_t i=0; i<Hsize; ++i){
		auto &Hi = H[i];
		auto &abpi = alphabetap[i];
		//alpha<-alphap beta<-betap
		for(std::size_t j=0u, jend=Hi.size(); j<jend; ++j){
			auto &hij = Hi[j];
			auto &Hj = HT[hij];
			std::size_t k=0;
			while(Hj[k]!=i) ++k;
			auto &[ai, bi] = alphabeta[k];
			abpi[j] = std::make_pair(&ai[hij],&bi[hij]);
		}
	}
}

template<CheckMatrix T>
Iterative_decoding<T>::Iterative_decoding(const T &Hp){
	static bool init;
	if(!init){
		H = Hp;
		init = true;
	}
}

template<CheckMatrix T>
void Iterative_decoding<T>::decode_init(){
	static thread_local bool init;
	if(!init){
		decodingData<T>::init();
		init = true;
	}
	for(auto &[ai, bi]: alphabeta) for(auto &bij: bi) bij = 0;
}

template<CheckMatrix T>
template<DecoderType D, std::floating_point U>
bool Iterative_decoding<T>::iterate(std::array<U,C> &LPR, const std::array<U,C> &LLR){
	//apply LLR
	for(auto &[ai, bi]: alphabeta) for(std::size_t j=0; j<C; ++j) bi[j] += LLR[j];
	//row update
	D::rowupdate();
	//column update
	for(auto &lpj: LPR) lpj = 0;
	for(auto &[ai, bi]: alphabeta) for(std::size_t j=0; j<C; ++j) LPR[j] += ai[j];
	for(auto &[ai, bi]: alphabeta) for(std::size_t j=0; j<C; ++j) bi[j] = LPR[j]-ai[j];
	for(std::size_t j=0; j<C; ++j) LPR[j] += LLR[j];
	//parity check
	for(const auto &Hi : H){
		auto parity = false;
		for(const auto &j : Hi) parity ^= LPR[j]<0;
		if(parity) return false;
	}
	return true;
}

template<CheckMatrix T>
template<std::floating_point U>
auto Iterative_decoding<T>::estimate(const std::array<U,C> &LEVR) const{
	std::bitset<S> est;
	for(std::size_t j=0; j<S; ++j) est[j] = LEVR[j]<0;
	return est;
}

template<CheckMatrix T>
void Iterative_decoding<T>::SumProduct::rowupdate(){
	for(auto &[ai, bi]: alphabeta) for(auto &bij: bi) bij = fg(bij);
	for(auto &abpi: alphabetap){
		fptype abssum = 0;
		signtype signprod = 0;
		for(const auto [apij,bpij]: abpi){
			auto bij = *bpij;
			abssum += std::fabs(bij);
			signprod ^= std::bit_cast<signtype>(bij);
		}
		for(const auto [apij,bpij]: abpi){
			auto bij = *bpij;
			auto absval = fg(abssum-std::fabs(bij));
			auto sign = (std::bit_cast<signtype>(bij)^signprod)&signmask;
			*apij = std::bit_cast<fptype>(sign|std::bit_cast<signtype>(absval));
		}
	}
}

template<CheckMatrix T>
void Iterative_decoding<T>::MinSum::rowupdate(){
	for(auto &abpi: alphabetap){
		signtype signprod = 0;
		for(const auto [apij,bpij]: abpi){
			signprod ^= std::bit_cast<signtype>(*bpij);
		}
		for(std::size_t j=0, jend=abpi.size(); j<jend; ++j){
			auto [apij,bpij] = abpi[j];
			auto min = std::numeric_limits<fptype>::infinity();
			for(std::size_t k=0u, kend=abpi.size(); k<kend; ++k) if(j != k){
				auto temp = std::fabs(*abpi[k].second);
				if(temp<min) min = temp;
			}
			auto sign = (std::bit_cast<const signtype>(*bpij)^signprod)&signmask;
			*apij = std::bit_cast<fptype>(sign|std::bit_cast<signtype>(min));
		}
	}
}

}

#endif