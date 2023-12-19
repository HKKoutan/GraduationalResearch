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
	static constexpr std::uint32_t S = T::sourcesize();
	static constexpr std::uint32_t C = T::codesize();
	static constexpr std::uint32_t Hsize = C-S;

	const T H;//検査行列
	const std::uint32_t Hones;
	const std::uint32_t VW;
	std::unique_ptr<fptype[]> alphabeta;
	inline static std::unique_ptr<uitype[]> alphabetaidx;

	void alphabetaidx_init();
public:
	explicit Sumproduct_decoding(const T &H);
	void decode_init();//decodeで使用する変数の初期化
	template<boxplusclass P>
	bool iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp);
	template<boxplusclass P>
	void decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, std::uint32_t iterationlimit);
};

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
class Sumproduct_decoding<CheckMatrix_regular<S,C,W>> {
	using T = CheckMatrix_regular<S,C,W>;
	static constexpr std::uint32_t Hsize = C-S;
	static constexpr std::uint32_t Hones = W*Hsize;
	static constexpr std::uint32_t VW = Hones/C;//列重み

	const T H;//検査行列
	std::unique_ptr<fptype[]> alphabeta;
	inline static std::unique_ptr<uitype[]> alphabetaidx;

	void alphabetaidx_init();
public:
	explicit Sumproduct_decoding(const T &H);
	void decode_init();//decodeで使用する変数の初期化
	template<boxplusclass P>
	bool iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp);
	template<boxplusclass P>
	void decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, std::uint32_t iterationlimit);
};

////////////////////////////////////////////////////////////////
//                                                            //
//                  class Sumproduct_decoding                 //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
void Sumproduct_decoding<T>::alphabetaidx_init(){
	alphabetaidx = std::make_unique<uitype[]>(Hones);

	for(std::uint32_t i=0; i<Hsize; ++i){
		auto &Hi = H[i];
		auto abpi = alphabetaidx.get()+H.headidxcol(i);
		//alpha<-alphap beta<-betap
		for(std::uint32_t j=0, jend=Hi.size(); j<jend; ++j){
			std::uint32_t hij = Hi[j];
			auto &Hj = H.T[hij];
			std::uint32_t k=0;
			while(Hj[k]!=i) ++k;
			abpi[j] = C*k+hij;
		}
	}
}

template<CheckMatrix T>
Sumproduct_decoding<T>::Sumproduct_decoding(const T &H):
	H(H),
	Hones(H.countones()),
	VW(H.weightcolmax()),
	alphabeta(std::make_unique<fptype[]>(C*VW))
{
	if(!alphabetaidx) alphabetaidx_init();
}

template<CheckMatrix T>
void Sumproduct_decoding<T>::decode_init(){
	for(std::uint32_t i=0, iend=C*VW; i<iend; ++i) alphabeta[i] = 0;
}

template<CheckMatrix T>
template<boxplusclass P>
bool Sumproduct_decoding<T>::iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp){
	//apply LLR
	for(std::uint32_t i=0; i<VW; ++i){
		auto bi = alphabeta.get()+C*i;
		for(std::uint32_t j=0; j<C; ++j) bi[j] += LLR[j];
	}
	//row update
	for(std::uint32_t i=0; i<VW; ++i){
		auto bi = alphabeta.get()+C*i;
		for(std::uint32_t j=0; j<C; ++j){
			if(i<H.weightcol(j)) bi[j] = bp.forward(bi[j]);
			else bi[j]=0;
		}
	}
	for(std::uint32_t i=0; i<Hsize; ++i){
		auto abxi = alphabetaidx.get()+H.headidxcol(i);
		auto abxiend = alphabetaidx.get()+H.headidxcol(i+1);
		accumlator_t<P,fptype> acc;
		for(auto abxij=abxi; abxij!=abxiend; ++abxij) acc += alphabeta[*abxij];
		for(auto abxij=abxi; abxij!=abxiend; ++abxij) alphabeta[*abxij] = acc-alphabeta[*abxij];
	}
	for(std::uint32_t i=0; i<VW; ++i){
		auto bi = alphabeta.get()+C*i;
		for(std::uint32_t j=0; j<C; ++j){
			if(i<H.weightcol(j)) bi[j] = bp.backward(bi[j]);
			else bi[j]=0;
		}
	}
	//column update
	for(auto &lpj: LPR) lpj = 0;
	for(std::uint32_t i=0; i<VW; ++i){
		auto ai = alphabeta.get()+C*i;
		for(std::uint32_t j=0; j<C; ++j) LPR[j] += ai[j];
	}
	for(std::uint32_t i=0; i<VW; ++i){
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

template<CheckMatrix T>
template<boxplusclass P>
void Sumproduct_decoding<T>::decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, std::uint32_t iterationlimit){
	decode_init();
	for(std::uint32_t iter=0; !iterate(LPR, LLR, bp) && iter<iterationlimit; ++iter);
}

////////////////////////////////////////////////////////////////
//                                                            //
//       class Sumproduct_decoding<CheckMatrix_regular>       //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::alphabetaidx_init(){
	alphabetaidx = std::make_unique<uitype[]>(Hones);

	for(std::uint32_t i=0; i<Hsize; ++i){
		auto &Hi = H[i];
		auto abxi = alphabetaidx.get()+H.headidxcol(i);
		//alpha<-alphap beta<-betap
		for(std::uint32_t j=0; j<W; ++j){
			std::uint32_t  hij = Hi[j];
			auto &Hj = H.T[hij];
			std::uint32_t  k=0;
			while(Hj[k]!=i) ++k;
			abxi[j] = C*k+hij;
		}
	}
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::Sumproduct_decoding(const T &H):
	H(H),
	alphabeta(std::make_unique<fptype[]>(Hones))
{
	if(!alphabetaidx) alphabetaidx_init();
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::decode_init(){
	// for(auto &bi: alphabeta) for(auto &bij: bi) bij = 0;
	for(std::uint32_t i=0; i<Hones; ++i) alphabeta[i] = 0;
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
template<boxplusclass P>
bool Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp){
	//apply LLR
	for(std::uint32_t i=0; i<VW; ++i){
		auto bi = alphabeta.get()+C*i;
		for(std::uint32_t j=0; j<C; ++j) bi[j] += LLR[j];
	}
	//row update
	for(std::uint32_t i=0; i<Hones; ++i) alphabeta[i] = bp.forward(alphabeta[i]);
	for(std::uint32_t i=0; i<Hsize; ++i){
		auto abxi = alphabetaidx.get()+H.headidxcol(i);
		accumlator_t<P,fptype> acc;
		for(std::uint32_t j=0; j<W; ++j) acc += alphabeta[abxi[j]];
		for(std::uint32_t j=0; j<W; ++j) alphabeta[abxi[j]] = acc-alphabeta[abxi[j]];
	}
	for(std::uint32_t i=0; i<Hones; ++i) alphabeta[i] = bp.backward(alphabeta[i]);
	//column update
	for(auto &lpj: LPR) lpj = 0;
	for(std::uint32_t i=0; i<VW; ++i){
		auto ai = alphabeta.get()+C*i;
		for(std::uint32_t j=0; j<C; ++j) LPR[j] += ai[j];
	}
	for(std::uint32_t i=0; i<VW; ++i){
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
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, std::uint32_t iterationlimit){
	decode_init();
	for(std::uint32_t iter=0; !iterate(LPR, LLR, bp) && iter<iterationlimit; ++iter);
}

}

#endif