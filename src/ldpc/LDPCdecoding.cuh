#ifndef INCLUDE_GUARD_ldpc_LDPCdecoding
#define INCLUDE_GUARD_ldpc_LDPCdecoding

#include <algorithm>
// #include <type_traits>
#include "LDPCCheckMatrix.cuh"
#include "LDPCboxplus.cuh"

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
	inline static std::unique_ptr<uitype[]> alphabetaidx;

	void alphabetaidx_init();
public:
	explicit Sumproduct_decoding(const T &H);
	std::size_t alphabetasize() const{return C*VW;}
	template<boxplusclass P>
	bool iterate(std::unique_ptr<fptype[]> &alphabeta, std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp);
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
	inline static std::unique_ptr<uitype[]> alphabetaidx;
	inline static std::unique_ptr<uitype[],util::cuda_delete> alphabetaidx_device;

	void alphabetaidx_init();
public:
	explicit Sumproduct_decoding(const T &H);
	static constexpr std::size_t alphabetasize(){return Hones;}
	template<boxplusclass P>
	bool iterate(std::unique_ptr<fptype[]> &alphabeta, std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp);
	template<boxplusclass P>
	void iterate(bool *parity, std::unique_ptr<fptype[],util::cuda_delete> &alphabeta, fptype *LPR, const fptype *LLR, const P &bp);
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
	VW(H.weightcolmax())
{
	if(!alphabetaidx) alphabetaidx_init();
}

template<CheckMatrix T>
template<boxplusclass P>
bool Sumproduct_decoding<T>::iterate(std::unique_ptr<fptype[]> &alphabeta, std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp){
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
	auto alphabeta = std::make_unique<fptype[]>(alphabetasize());
	for(std::uint32_t iter=0; !iterate(alphabeta, LPR, LLR, bp) && iter<iterationlimit; ++iter);
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

	alphabetaidx_device = util::make_cuda_unique<uitype[]>(Hones);
	util::check_cuda_error(cudaMemcpy(alphabetaidx_device.get(), alphabetaidx.get(), sizeof(uitype)*Hones, cudaMemcpyHostToDevice));
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::Sumproduct_decoding(const T &H):
	H(H)
{
	if(!alphabetaidx) alphabetaidx_init();
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
template<boxplusclass P>
bool Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::iterate(std::unique_ptr<fptype[]> &alphabeta, std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp){
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

namespace {
	template<boxplusclass P, CheckMatrix U>
	__global__ void parallel_decode(
		fptype *alphabeta,
		uitype *alphabetaidx,
		fptype *LPR,
		const fptype *LLR,
		const P &bp,
		typename P::internaldatatype bpd,
		U H,
		typename U::internaldatatype Hd,
		std::uint32_t iterationlimit
	){
		constexpr std::uint32_t Hsize = U::size();
		constexpr std::uint32_t W = U::weightrowmax(Hd);
		constexpr std::uint32_t C = U::codesize();
		constexpr std::uint32_t VW = U::weightcolmax(Hd);
		constexpr std::uint32_t Hones = U::countones();
		constexpr std::uint32_t tsize = 1024;
		const std::uint32_t tid = threadIdx.x;
		assert(blockDim.x == tsize);
		__shared__ bool parity;

		if(tid==0) parity = true;
		__syncthreads();

		for(std::uint32_t itr=0; itr<iterationlimit&&parity; ++itr){
			//apply LLR & boxplus forwarding
			for(auto i=tid; i<Hones; i+=tsize) alphabeta[i] = P::forward(alphabeta[i]+LLR[i%C],bpd);
			__syncthreads();
			//rowupdate
			for(auto i=tid; i<Hsize; i+=tsize){
				accumlator_t<P,fptype> acc;
				const auto jbegin = U::headidxcol(i,Hd), jend = U::headidxcol(i+1,Hd);
				for(auto j=jbegin; j<jend; ++j) acc += alphabeta[alphabetaidx[j]];
				for(auto j=jbegin; j<jend; ++j) alphabeta[alphabetaidx[j]] = acc-alphabeta[alphabetaidx[j]];
			}
			//余ったスレッドでLPR初期化
			for(auto i=tsize-tid-1; i<C; i+=tsize) LPR[i] = 0;
			__syncthreads();
			//boxplus backwarding
			for(auto i=tid; i<Hones; i+=tsize) alphabeta[i] = P::backward(alphabeta[i],bpd);
			__syncthreads();
			//column update
			for(auto i=tid; i<C; i+=tsize) for(auto j=i; j<Hones; j+=C) LPR[i] += alphabeta[j];
			__syncthreads();
			for(auto i=tid; i<Hones; i+=tsize) alphabeta[i] = LPR[i%C]-alphabeta[i];
			__syncthreads();
			for(auto i=tsize-tid-1; i<C; i+=tsize) LPR[i] += LLR[i];
			__syncthreads();
			//check parity
			if(tid==0) parity = 0;
			for(auto i=tid, iend=((Hsize-1)&0xfffffe00)+tsize; i<iend; i+=tsize){
				bool parity_row = false;
				if(i<Hsize) for(auto j: U::getrow(i,Hd)) parity_row^=LPR[j]<0;
				parity_row = __syncthreads_or(parity_row);
				if(tid==0) parity |= parity_row;
			}
			__syncthreads();
		}
	}
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
template<boxplusclass P>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, std::uint32_t iterationlimit){
	auto alphabeta = util::make_cuda_unique<fptype[]>(alphabetasize());
	util::check_cuda_error(::cudaMemset(alphabeta.get(),0,alphabetasize()));

	fptype *LPR_device;//対数事後確率比：列ごとのalphaの和+QLLR
	fptype *LLR_device;
	util::check_cuda_error(::cudaMalloc(&LPR_device, sizeof(fptype)*C));
	util::check_cuda_error(::cudaMalloc(&LLR_device, sizeof(fptype)*C));
	util::check_cuda_error(::cudaMemcpy(LLR_device,LLR.data(),sizeof(fptype)*C,cudaMemcpyHostToDevice));

	parallel_decode<<<1,1024>>>(alphabeta.get(), alphabetaidx_device.get(), LPR_device, LLR_device, bp, bp.data(), H, H.data(), iterationlimit);

	util::check_cuda_error(::cudaMemcpy(LPR.data(),LPR_device,sizeof(fptype)*C,cudaMemcpyDeviceToHost));
	util::check_cuda_error(::cudaFree(LPR_device));
	util::check_cuda_error(::cudaFree(LLR_device));
}

}

#endif