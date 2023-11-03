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
	using masktype = std::uint32_t;
	static_assert(sizeof(fptype)==sizeof(masktype));
	static constexpr masktype signmask = 1u<<(sizeof(masktype)*8u-1u);
	static constexpr std::size_t S = T::sourcesize();
	static constexpr std::size_t C = T::codesize();
	static constexpr std::size_t Hsize = C-S;

	inline static T H;//検査行列
	//TODO: alphaとbetaを共有
	inline static thread_local std::vector<std::array<fptype,C>> alphabeta;
	inline static thread_local std::array<std::vector<fptype*>,C-S> alphabetap;
	inline static thread_local std::vector<std::array<masktype,C>> alphabetamask;//alphabetaの有効な要素に対応する要素が1埋めされた配列
	//メンバ初期化関数
	static auto alphabeta_height();
	static auto makeHT();
	static auto alphabetap_init(std::array<std::vector<std::uint64_t>,C> &HT);
	static auto alphabetamask_init(std::array<std::vector<std::uint64_t>,C> &HT);
public:
	explicit Iterative_decoding(const T &H);
	void decode_init();//decodeで使用する変数の初期化
	template<DecoderType D, std::floating_point U>
	bool iterate(std::array<U,C> &LPR, const std::array<U,C> &LLR);
	template<std::floating_point U>
	auto estimate(const std::array<U,C> &LEVR) const;//推定符号語を求める
	//rowupdate
	struct SumProduct {
		inline static const func_Gallager_table fg;
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
auto Iterative_decoding<T>::alphabeta_height(){
	std::array<std::size_t,C> Hheight;
	for(const auto &Hi: H) for(auto j: Hi) ++Hheight[j];
	return std::ranges::max(Hheight);
}

template<CheckMatrix T>
auto Iterative_decoding<T>::makeHT(){
	std::array<std::vector<std::uint64_t>,C> HT;
	for(std::size_t i=0; i<Hsize; ++i) for(auto j: H[i]) HT[j].push_back(i);
	return HT;
}

template<CheckMatrix T>
auto Iterative_decoding<T>::alphabetap_init(std::array<std::vector<std::uint64_t>,C> &HT){
	decltype(alphabetap) abp{};
	for(std::size_t i=0; i<Hsize; ++i){
		auto &Hi = H[i];
		auto &abpi = abp[i];
		//Hとalphabetapの要素数を揃える
		abpi.resize(Hi.size());
		//alphabeta<-alphabetap
		for(std::size_t j=0u, jend=Hi.size(); j<jend; ++j){
			auto &hij = Hi[j];
			auto &Hj = HT[hij];
			auto &abi = alphabeta[std::ranges::find(Hj, i)-Hj.begin()];
			abpi[j] = &abi[hij];
		}
	}
	return abp;
}

template<CheckMatrix T>
auto Iterative_decoding<T>::alphabetamask_init(std::array<std::vector<std::uint64_t>,C> &HT){
	decltype(alphabetamask) abm(alphabeta.size());
	for(std::size_t j=0; j<C; ++j) for(std::size_t i=0, iend=HT[j].size(); i<iend; ++i) abm[i][j]=0xffffffff;
	return abm;
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
		alphabeta.resize(alphabeta_height());
		auto HT = makeHT();
		alphabetap = alphabetap_init(HT);
		alphabetamask = alphabetamask_init(HT);
		init = true;
	}
	for(auto &bi: alphabeta) for(auto &bij: bi) bij = 0;
}

template<CheckMatrix T>
template<DecoderType D, std::floating_point U>
bool Iterative_decoding<T>::iterate(std::array<U,C> &LPR, const std::array<U,C> &LLR){
	//apply LLR
	for(auto &bi: alphabeta) for(std::size_t j=0; j<C; ++j) bi[j] += LLR[j];
	//row update
	D::rowupdate();
	for(std::size_t i=0, iend=alphabeta.size(); i<iend;++i){
		auto &ai = alphabeta[i];
		auto &mi = alphabetamask[i];
		for(std::size_t j=0; j<C; ++j) ai[j] = std::bit_cast<fptype>(std::bit_cast<masktype>(ai[j])&mi[j]);
	}
	//column update
	for(auto &lpj: LPR) lpj = 0;
	for(auto &ai: alphabeta) for(std::size_t j=0; j<C; ++j) LPR[j] += ai[j];
	for(auto &abi: alphabeta) for(std::size_t j=0; j<C; ++j) abi[j] = LPR[j]-abi[j];
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
	for(std::size_t j=0u, jend=S; j<jend; ++j) est[j] = LEVR[j]<0;
	return est;
}

template<CheckMatrix T>
void Iterative_decoding<T>::SumProduct::rowupdate(){
	for(auto &bi: alphabeta) for(auto &bij: bi) bij = fg(bij);
	for(auto &abpi: alphabetap){
		fptype abssum = 0;
		masktype signprod = 0;
		for(auto bpij: abpi){
			auto bij = *bpij;
			abssum += std::fabs(bij);
			signprod ^= std::bit_cast<masktype>(bij);
		}
		for(std::size_t j=0u, jend=abpi.size(); j<jend; ++j){
			auto bij = *abpi[j];
			auto absval = fg(abssum-std::fabs(bij));
			auto sign = (std::bit_cast<masktype>(bij)^signprod)&signmask;
			*abpi[j] = std::bit_cast<fptype>(sign|std::bit_cast<masktype>(absval));
		}
	}
}

template<CheckMatrix T>
void Iterative_decoding<T>::MinSum::rowupdate(){
	for(auto &abpi: alphabetap){
		masktype signprod = 0;
		for(auto bpij: abpi){
			signprod ^= std::bit_cast<masktype>(*bpij);
		}
		for(std::size_t j=0u, jend=abpi.size(); j<jend; ++j){
			auto min = std::numeric_limits<fptype>::infinity();
			for(std::size_t k=0u, kend=abpi.size(); k<kend; ++k) if(j != k){
				auto temp = std::fabs(*abpi[k]);
				if(temp<min) min = temp;
			}
			auto sign = (std::bit_cast<const masktype>(*abpi[j])^signprod)&signmask;
			*abpi[j] = std::bit_cast<fptype>(sign|std::bit_cast<masktype>(min));
		}
	}
}

}

#endif