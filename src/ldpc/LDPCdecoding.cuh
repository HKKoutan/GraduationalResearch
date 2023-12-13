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
#include "LDPCboxplus.hpp"

namespace code::LDPC {

template<CheckMatrix T>
class Iterative_decoding {
	using fptype = float;
	static constexpr std::size_t S = T::sourcesize();
	static constexpr std::size_t C = T::codesize();

	inline static T H;//検査行列

	inline static thread_local std::vector<std::pair<std::array<fptype,C>,std::array<fptype,C>>> alphabeta;
	inline static thread_local std::array<std::vector<std::pair<fptype*,const fptype*>>,C-S> alphabetap;
	static void init();
public:
	explicit Iterative_decoding(const T &H);
	void decode_init();//decodeで使用する変数の初期化
	template<std::floating_point U, boxplusclass P>
	bool iterate(std::array<U,C> &LPR, const std::array<U,C> &LLR, const P &bp);
};

template<std::size_t S, std::size_t C, std::size_t W>
class Iterative_decoding<CheckMatrix_regular<S,C,W>> {
	using fptype = float;

	inline static CheckMatrix_regular<S,C,W> H;//検査行列

	inline static thread_local std::array<std::array<fptype,C>,W*S/C> alphabeta;
	inline static thread_local std::array<std::array<fptype*,W>,C-S> alphabetap;
	static void init();
public:
	explicit Iterative_decoding(const CheckMatrix_regular<S,C,W> &H);
	void decode_init();//decodeで使用する変数の初期化
	template<std::floating_point U, boxplusclass P>
	bool iterate(std::array<U,C> &LPR, const std::array<U,C> &LLR, const P &bp);
};

////////////////////////////////////////////////////////////////
//                                                            //
//                  class Iterative_decoding                  //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
void Iterative_decoding<T>::init(){
	constexpr std::size_t Hsize = C-S;
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
Iterative_decoding<T>::Iterative_decoding(const T &Hp){
	static bool init;
	if(!init){
		H = Hp;
		init = true;
	}
}

template<CheckMatrix T>
void Iterative_decoding<T>::decode_init(){
	static thread_local bool initialized;
	if(!initialized){
		init();
		initialized = true;
	}
	for(auto &[ai, bi]: alphabeta) for(auto &bij: bi) bij = 0;
}

template<CheckMatrix T>
template<std::floating_point U, boxplusclass P>
bool Iterative_decoding<T>::iterate(std::array<U,C> &LPR, const std::array<U,C> &LLR, const P &bp){
	//apply LLR
	for(auto &[ai, bi]: alphabeta) for(std::size_t j=0; j<C; ++j) bi[j] += LLR[j];
	//row update
	for(auto &abpi: alphabetap){
		accumlator_t<P,fptype> acc;
		for(const auto [apij,bpij]: abpi) acc += bp.forward(*bpij);
		for(const auto [apij,bpij]: abpi) *apij = bp.backward(acc-bp.forward(*bpij));
	}
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

////////////////////////////////////////////////////////////////
//                                                            //
//       class Iterative_decoding<CheckMatrix_regular>        //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C, std::size_t W>
void Iterative_decoding<CheckMatrix_regular<S,C,W>>::init(){
	constexpr std::size_t Hsize = C-S;
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
			auto &abi = alphabeta[k];
			abpi[j] = &abi[hij];
		}
	}
}

template<std::size_t S, std::size_t C, std::size_t W>
Iterative_decoding<CheckMatrix_regular<S,C,W>>::Iterative_decoding(const CheckMatrix_regular<S,C,W> &Hp){
	static bool init;
	if(!init){
		H = Hp;
		init = true;
	}
}

template<std::size_t S, std::size_t C, std::size_t W>
void Iterative_decoding<CheckMatrix_regular<S,C,W>>::decode_init(){
	static thread_local bool initialized;
	if(!initialized){
		init();
		initialized = true;
	}
	for(auto &bi: alphabeta) for(auto &bij: bi) bij = 0;
}

template<std::size_t S, std::size_t C, std::size_t W>
template<std::floating_point U, boxplusclass P>
bool Iterative_decoding<CheckMatrix_regular<S,C,W>>::iterate(std::array<U,C> &LPR, const std::array<U,C> &LLR, const P &bp){
	//apply LLR
	for(auto &bi: alphabeta) for(std::size_t j=0; j<C; ++j) bi[j] += LLR[j];
	//row update
	for(auto &bi: alphabeta) for(auto &bij: bi) bij = bp.forward(bij);
	for(auto &abpi: alphabetap){
		accumlator_t<P,fptype> acc;
		for(const auto bpij: abpi) acc += *bpij;
		for(const auto abpij: abpi) *abpij = acc-*abpij;
	}
	for(auto &ai: alphabeta) for(auto &aij: ai) aij = bp.backward(aij);
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

// template<CheckMatrix T>
// void Iterative_decoding<T>::MinSum::rowupdate(){
// 	using signtype = std::uint32_t;
// 	static_assert(sizeof(fptype)==sizeof(signtype));
// 	static constexpr signtype signmask = 1u<<(sizeof(signtype)*8u-1u);
// 	for(auto &abpi: alphabetap){
// 		signtype signprod = 0;
// 		for(const auto [apij,bpij]: abpi){
// 			signprod ^= std::bit_cast<signtype>(*bpij);
// 		}
// 		for(std::size_t j=0, jend=abpi.size(); j<jend; ++j){
// 			auto [apij,bpij] = abpi[j];
// 			auto min = std::numeric_limits<fptype>::infinity();
// 			for(std::size_t k=0u, kend=abpi.size(); k<kend; ++k) if(j != k){
// 				auto temp = std::fabs(*abpi[k].second);
// 				if(temp<min) min = temp;
// 			}
// 			auto sign = (std::bit_cast<const signtype>(*bpij)^signprod)&signmask;
// 			*apij = std::bit_cast<fptype>(sign|std::bit_cast<signtype>(min));
// 		}
// 	}
// }

}

#endif