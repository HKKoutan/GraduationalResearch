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
class Sumproduct_decoding {
	using fptype = float;
	static constexpr std::size_t S = T::sourcesize();
	static constexpr std::size_t C = T::codesize();

	const T H;//検査行列
	std::vector<std::pair<std::array<fptype,C>,std::array<fptype,C>>> alphabeta;
	std::array<std::vector<std::pair<fptype*,const fptype*>>,C-S> alphabetap;

	static auto alphabeta_size(const T &H);
	static auto alphabetap_init(const T &H, std::vector<std::pair<std::array<fptype,C>,std::array<fptype,C>>> &alphabeta);
public:
	explicit Sumproduct_decoding(const T &H);
	void decode_init();//decodeで使用する変数の初期化
	template<boxplusclass P>
	bool iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp);
};

template<std::size_t S, std::size_t C, std::size_t W>
class Sumproduct_decoding<CheckMatrix_regular<S,C,W>> {
	using fptype = float;
	static constexpr std::size_t Hsize = C-S;
	static constexpr std::size_t Hones = W*Hsize;
	static constexpr std::size_t VW = Hones/C;//列重み

	const CheckMatrix_regular<S,C,W> H;//検査行列
	std::unique_ptr<fptype[][C]> alphabeta;
	std::unique_ptr<fptype*[][W]> alphabetap;
	// std::array<std::array<fptype,C>,VW> alphabeta;
	// std::array<std::array<fptype*,W>,Hsize> alphabetap;

	void alphabetap_init();
public:
	explicit Sumproduct_decoding(const CheckMatrix_regular<S,C,W> &H);
	void decode_init();//decodeで使用する変数の初期化
	template<boxplusclass P>
	bool iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp);
};

////////////////////////////////////////////////////////////////
//                                                            //
//                  class Sumproduct_decoding                 //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
auto Sumproduct_decoding<T>::alphabeta_size(const T &H){
	constexpr std::size_t Hsize = C-S;
	std::array<std::size_t,C> Hheight{};
	for(std::size_t i=0; i<Hsize; ++i) for(auto j: H[i]) ++Hheight[j];
	return std::ranges::max(Hheight);
}

template<CheckMatrix T>
auto Sumproduct_decoding<T>::alphabetap_init(const T &H, std::vector<std::pair<std::array<fptype,C>,std::array<fptype,C>>> &alphabeta){
	constexpr std::size_t Hsize = C-S;
	std::array<std::vector<std::pair<fptype*,const fptype*>>,C-S> alphabetap;

	// std::array<std::vector<std::uint64_t>,C> HT{};//Hの転置
	// for(std::size_t i=0; i<Hsize; ++i) for(auto j: H[i]) HT[j].push_back(i);

	for(std::size_t i=0; i<Hsize; ++i){
		auto &Hi = H[i];
		auto &abpi = alphabetap[i];
		//Hとalphabetapの要素数を揃える
		abpi.resize(Hi.size());
		//alpha<-alphap beta<-betap
		for(std::size_t j=0, jend=Hi.size(); j<jend; ++j){
			auto &hij = Hi[j];
			auto &Hj = H.T[hij];
			std::size_t k=0;
			while(Hj[k]!=i) ++k;
			auto &[ai, bi] = alphabeta[k];
			abpi[j] = std::make_pair(&ai[hij],&bi[hij]);
		}
	}
	return alphabetap;
}

template<CheckMatrix T>
Sumproduct_decoding<T>::Sumproduct_decoding(const T &H):
	H(H),
	alphabeta(alphabeta_size(H)),
	alphabetap(alphabetap_init(H,alphabeta))
{}

template<CheckMatrix T>
void Sumproduct_decoding<T>::decode_init(){
	for(auto &[ai, bi]: alphabeta) for(auto &bij: bi) bij = 0;
}

template<CheckMatrix T>
template<boxplusclass P>
bool Sumproduct_decoding<T>::iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp){
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
	constexpr std::size_t Hsize = C-S;
	for(std::size_t i=0; i<Hsize; ++i){
		auto parity = false;
		for(const auto &j : H[i]) parity ^= LPR[j]<0;
		if(parity) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////
//                                                            //
//       class Sumproduct_decoding<CheckMatrix_regular>       //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C, std::size_t W>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::alphabetap_init(){
	// std::array<std::array<std::size_t,VW>,C> HT;//Hの転置
	// {
	// 	std::array<std::size_t,C> HTc = {};
	// 	for(std::size_t i=0; i<Hsize; ++i) for(auto j: H[i]) HT[j][HTc[j]++] = i;
	// }
	for(std::size_t i=0; i<Hsize; ++i){
		auto &Hi = H[i];
		auto &abpi = alphabetap[i];
		//alpha<-alphap beta<-betap
		for(std::size_t j=0; j<W; ++j){
			auto &hij = Hi[j];
			auto &Hj = H.T[hij];
			std::size_t k=0;
			while(Hj[k]!=i) ++k;
			auto &abk = alphabeta[k];
			abpi[j] = &abk[hij];
		}
	}
}

template<std::size_t S, std::size_t C, std::size_t W>
Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::Sumproduct_decoding(const CheckMatrix_regular<S,C,W> &H):
	H(H),
	alphabeta(std::make_unique<fptype[][C]>(W)),
	alphabetap(std::make_unique<fptype*[][W]>(Hsize))
{
	alphabetap_init();
}

template<std::size_t S, std::size_t C, std::size_t W>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::decode_init(){
	// for(auto &bi: alphabeta) for(auto &bij: bi) bij = 0;
	for(auto i=0; i<VW; ++i) for(auto &bij: alphabeta[i]) bij = 0;
}

template<std::size_t S, std::size_t C, std::size_t W>
template<boxplusclass P>
bool Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp){
	//apply LLR
	for(auto i=0; i<VW; ++i){
		auto &bi = alphabeta[i];
		for(std::size_t j=0; j<C; ++j) bi[j] += LLR[j];
	}
	//row update
	for(auto i=0; i<VW; ++i) for(auto &bij: alphabeta[i]) bij = bp.forward(bij);
	for(auto i=0; i<Hsize; ++i){
		auto &abpi = alphabetap[i];
		accumlator_t<P,fptype> acc;
		for(const auto bpij: abpi) acc += *bpij;
		for(const auto abpij: abpi) *abpij = acc-*abpij;
	}
	for(auto i=0; i<VW; ++i) for(auto &aij: alphabeta[i]) aij = bp.backward(aij);
	//column update
	for(auto &lpj: LPR) lpj = 0;
	for(auto i=0; i<VW; ++i){
		auto &ai = alphabeta[i];
		for(std::size_t j=0; j<C; ++j) LPR[j] += ai[j];
	}
	for(auto i=0; i<VW; ++i){
		auto &abi = alphabeta[i];
		for(std::size_t j=0; j<C; ++j) abi[j] = LPR[j]-abi[j];
	}
	for(std::size_t j=0; j<C; ++j) LPR[j] += LLR[j];
	//parity check
	for(std::size_t i=0; i<Hsize; ++i){
		auto parity = false;
		for(const auto &j : H[i]) parity ^= LPR[j]<0;
		if(parity) return false;
	}
	return true;
}

// template<std::size_t S, std::size_t C, std::size_t W>
// template<std::floating_point U, boxplusclass P>
// bool Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::iterate(std::array<U,C> &LPR, const std::array<U,C> &LLR, const P &bp){
// 	//apply LLR
// 	for(auto &bi: alphabeta) for(std::size_t j=0; j<C; ++j) bi[j] += LLR[j];
// 	//row update
// 	for(auto &bi: alphabeta) for(auto &bij: bi) bij = bp.forward(bij);
// 	for(auto &abpi: alphabetap){
// 		accumlator_t<P,fptype> acc;
// 		for(const auto bpij: abpi) acc += *bpij;
// 		for(const auto abpij: abpi) *abpij = acc-*abpij;
// 	}
// 	for(auto &ai: alphabeta) for(auto &aij: ai) aij = bp.backward(aij);
// 	//column update
// 	for(auto &lpj: LPR) lpj = 0;
// 	for(auto &ai: alphabeta) for(std::size_t j=0; j<C; ++j) LPR[j] += ai[j];
// 	for(auto &abi: alphabeta) for(std::size_t j=0; j<C; ++j) abi[j] = LPR[j]-abi[j];
// 	for(std::size_t j=0; j<C; ++j) LPR[j] += LLR[j];
// 	//parity check
// 	constexpr std::size_t Hsize = C-S;
// 	for(std::size_t i=0; i<Hsize; ++i){
// 		auto parity = false;
// 		for(const auto &j : H[i]) parity ^= LPR[j]<0;
// 		if(parity) return false;
// 	}
// 	return true;
// }

// template<CheckMatrix T>
// void Sumproduct_decoding<T>::MinSum::rowupdate(){
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