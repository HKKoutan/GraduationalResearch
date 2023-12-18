#ifndef INCLUDE_GUARD_ldpc_LDPCdecoding
#define INCLUDE_GUARD_ldpc_LDPCdecoding

#include <algorithm>
#include "LDPCCheckMatrix.hpp"
#include "LDPCboxplus.hpp"

namespace code::LDPC {

namespace {
	using fptype = float;
	using uitype = unsigned int;
	static_assert(sizeof(fptype)==sizeof(uitype));
}

template<CheckMatrix T>
class Sumproduct_decoding {
	// using fptype = float;
	static constexpr std::uint32_t S = T::sourcesize();
	static constexpr std::uint32_t C = T::codesize();
	static constexpr std::uint32_t Hsize = C-S;

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
	template<boxplusclass P>
	void decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, int iterationlimit);
};

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
class Sumproduct_decoding<CheckMatrix_regular<S,C,W>> {
	using T = CheckMatrix_regular<S,C,W>;
	// using fptype = float;
	static constexpr std::uint32_t Hsize = C-S;
	static constexpr std::uint32_t Hones = W*Hsize;
	static constexpr std::uint32_t VW = Hones/C;//列重み
	static_assert(Hones<(1ui64<<32));

	const T H;//検査行列
	std::unique_ptr<fptype[]> alphabeta;
	std::unique_ptr<uitype[]> alphabetaidx;
	// std::unique_ptr<fptype*[]> alphabetap;
	// std::array<std::array<fptype,C>,VW> alphabeta;
	// std::array<std::array<fptype*,W>,Hsize> alphabetap;

	void alphabetap_init();
public:
	explicit Sumproduct_decoding(const T &H);
	void decode_init();//decodeで使用する変数の初期化
	template<boxplusclass P>
	bool iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp);
	template<boxplusclass P>
	void decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, int iterationlimit);
};

////////////////////////////////////////////////////////////////
//                                                            //
//                  class Sumproduct_decoding                 //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
auto Sumproduct_decoding<T>::alphabeta_size(const T &H){
	std::array<std::uint32_t,C> Hheight{};
	for(std::uint32_t i=0; i<C; ++i) Hheight[i] = H.colweight(i);
	return std::ranges::max(Hheight);
}

template<CheckMatrix T>
auto Sumproduct_decoding<T>::alphabetap_init(const T &H, std::vector<std::pair<std::array<fptype,C>,std::array<fptype,C>>> &alphabeta){
	std::array<std::vector<std::pair<fptype*,const fptype*>>,C-S> alphabetap;

	// std::array<std::vector<std::uint64_t>,C> HT{};//Hの転置
	// for(std::uint32_t i=0; i<Hsize; ++i) for(auto j: H[i]) HT[j].push_back(i);

	for(std::uint32_t i=0; i<Hsize; ++i){
		auto &Hi = H[i];
		auto &abpi = alphabetap[i];
		//Hとalphabetapの要素数を揃える
		abpi.resize(Hi.size());
		//alpha<-alphap beta<-betap
		for(std::uint32_t j=0, jend=Hi.size(); j<jend; ++j){
			std::uint32_t hij = Hi[j];
			auto &Hj = H.T[hij];
			std::uint32_t k=0;
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
	for(auto &[ai, bi]: alphabeta) for(std::uint32_t j=0; j<C; ++j) bi[j] += LLR[j];
	//row update
	for(auto &abpi: alphabetap){
		accumlator_t<P,fptype> acc;
		for(const auto [apij,bpij]: abpi) acc += bp.forward(*bpij);
		for(const auto [apij,bpij]: abpi) *apij = bp.backward(acc-bp.forward(*bpij));
	}
	//column update
	for(auto &lpj: LPR) lpj = 0;
	for(auto &[ai, bi]: alphabeta) for(std::uint32_t j=0; j<C; ++j) LPR[j] += ai[j];
	for(auto &[ai, bi]: alphabeta) for(std::uint32_t j=0; j<C; ++j) bi[j] = LPR[j]-ai[j];
	for(std::uint32_t j=0; j<C; ++j) LPR[j] += LLR[j];
	//parity check
	for(std::uint32_t i=0; i<Hsize; ++i){
		auto parity = false;
		for(const auto &j : H[i]) parity ^= LPR[j]<0;
		if(parity) return false;
	}
	return true;
}

template<CheckMatrix T>
template<boxplusclass P>
void Sumproduct_decoding<T>::decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, int iterationlimit){
	decode_init();
	for(int iter=0; !iterate(LPR, LLR, bp) && iter<iterationlimit; ++iter);
}

////////////////////////////////////////////////////////////////
//                                                            //
//       class Sumproduct_decoding<CheckMatrix_regular>       //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::alphabetap_init(){
	for(std::uint32_t i=0; i<Hsize; ++i){
		auto &Hi = H[i];
		auto abxi = alphabetaidx.get()+W*i;
		//alpha<-alphap beta<-betap
		for(std::uint32_t j=0; j<W; ++j){
			uitype hij = Hi[j];
			auto &Hj = H.T[hij];
			int k=0;
			while(Hj[k]!=i) ++k;
			abxi[j] = static_cast<uitype>(C)*k+hij;
		}
	}
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::Sumproduct_decoding(const T &H):
	H(H),
	alphabeta(std::make_unique<fptype[]>(Hones)),
	alphabetaidx(std::make_unique<uitype[]>(Hones))
{
	alphabetap_init();
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::decode_init(){
	// for(auto &bi: alphabeta) for(auto &bij: bi) bij = 0;
	for(auto i=0; i<Hones; ++i) alphabeta[i] = 0;
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
template<boxplusclass P>
bool Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp){
	//apply LLR
	for(auto i=0; i<VW; ++i){
		auto bi = alphabeta.get()+C*i;
		for(std::uint32_t j=0; j<C; ++j) bi[j] += LLR[j];
	}
	//row update
	for(auto i=0; i<Hones; ++i) alphabeta[i] = bp.forward(alphabeta[i]);
	for(auto i=0; i<Hsize; ++i){
		auto abxi = alphabetaidx.get()+W*i;
		accumlator_t<P,fptype> acc;
		for(std::uint32_t j=0; j<W; ++j) acc += alphabeta[abxi[j]];
		for(std::uint32_t j=0; j<W; ++j) alphabeta[abxi[j]] = acc-alphabeta[abxi[j]];
	}
	for(auto i=0; i<Hones; ++i) alphabeta[i] = bp.backward(alphabeta[i]);
	//column update
	for(auto &lpj: LPR) lpj = 0;
	for(auto i=0; i<VW; ++i){
		auto ai = alphabeta.get()+C*i;
		for(std::uint32_t j=0; j<C; ++j) LPR[j] += ai[j];
	}
	for(auto i=0; i<VW; ++i){
		auto abi = alphabeta.get()+C*i;
		for(std::uint32_t j=0; j<C; ++j) abi[j] = LPR[j]-abi[j];
	}
	for(std::uint32_t j=0; j<C; ++j) LPR[j] += LLR[j];
	//parity check
	for(std::uint32_t i=0; i<Hsize; ++i){
		auto parity = false;
		for(const auto &j : H[i]) parity ^= LPR[j]<0;
		if(parity) return false;
	}
	return true;
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
template<boxplusclass P>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, int iterationlimit){
	decode_init();
	for(int iter=0; !iterate(LPR, LLR, bp) && iter<iterationlimit; ++iter);
}

}

#endif