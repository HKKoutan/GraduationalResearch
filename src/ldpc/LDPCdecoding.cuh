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
	inline static std::unique_ptr<uitype[],util::cuda_delete> alphabetaidx;

	void alphabetaidx_init();
public:
	explicit Sumproduct_decoding(const T &H);
	static constexpr std::size_t alphabetasize(){return Hones;}
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

namespace {
	template<CheckMatrix U>
	__global__ void parallel_idx_init(uitype *idx, U H, typename U::internaldatatype Hd){
		std::uint32_t i = blockIdx.x*blockDim.x+threadIdx.x;
		std::uint32_t j = blockIdx.y*blockDim.y+threadIdx.y;
		std::uint32_t W = U::weightrowmax(Hd);
		if(i<H.size()&&j<W){
			auto &Hi = U::getrow(i,Hd);
			std::uint32_t hij = Hi[j];
			auto &Hj = U::getcol(hij,Hd);
			std::uint32_t k=0;
			while(Hj[k]!=i) ++k;
			idx[W*i+j] = H.codesize()*k+hij;
		}
	}
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::alphabetaidx_init(){
	alphabetaidx = util::make_cuda_unique<uitype[]>(Hones);

	const dim3 grid(((Hsize-1)/128+1),((W-1)/8+1),1);
	const dim3 thread(128,8,1);
	parallel_idx_init<<<grid,thread>>>(alphabetaidx.get(), H, H.data());
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::Sumproduct_decoding(const T &H):
	H(H)
{
	if(!alphabetaidx) alphabetaidx_init();
}

namespace {
	__global__ void parallel_broadcast_add(fptype *lhs, const fptype *rhs, std::size_t width, std::size_t height){//lhs[i][j]+=rhs[j] for all i<height, j<width
		std::uint32_t j = blockIdx.x*blockDim.x+threadIdx.x;
		std::uint32_t i = blockIdx.y*blockDim.y+threadIdx.y;
		if(i<height&&j<width) lhs[i*width+j] += rhs[j];
	}
	__global__ void parallel_broadcast_negsub(fptype *lhs, const fptype *rhs, std::size_t width, std::size_t height){//lhs[i][j]=rhs[j]-lhs[i][j] for all i<height, j<width
		std::uint32_t j = blockIdx.x*blockDim.x+threadIdx.x;
		std::uint32_t i = blockIdx.y*blockDim.y+threadIdx.y;
		if(i<height&&j<width) lhs[i*width+j] = rhs[j]-lhs[i*width+j];
	}
	__global__ void parallel_reduce_add(fptype *lhs, const fptype *rhs, std::size_t width, std::size_t height){//lhs[j]+=rhs[i][j] for all i<height, j<width
		std::uint32_t j = blockIdx.x*blockDim.x+threadIdx.x;
		std::uint32_t k = blockIdx.y*blockDim.y+threadIdx.y;
		if(k==0&&j<width){
			for(std::uint32_t i=0; i<height; ++i) lhs[j] += rhs[i*width+j];
		}
	}
	template<class P>
	__global__ void parallel_accumlate(fptype *arr, uitype *idx, std::size_t width, std::size_t height, const P &bp){
		std::uint32_t j = blockIdx.x*blockDim.x+threadIdx.x;
		std::uint32_t k = blockIdx.y*blockDim.y+threadIdx.y;
		if(k==0&&j<height){
			accumlator_t<P,fptype> acc;
			for(std::uint32_t i=0; i<width; ++i) acc += arr[idx[j*width+i]];
			for(std::uint32_t i=0; i<width; ++i) arr[idx[j*width+i]] = acc-arr[idx[j*width+i]];
		}
	}
	// template<typename T,CheckMatrix U>
	// __global__ void check_parity(T* ptr, U H, typename U::internaldatatype Hd){
	// 	std::uint32_t j = blockIdx.x*blockDim.x+threadIdx.x;
	// 	std::uint32_t k = blockIdx.y*blockDim.y+threadIdx.y;
	// }
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
template<boxplusclass P>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::iterate(bool *parity, std::unique_ptr<fptype[],util::cuda_delete> &alphabeta, fptype *LPR, const fptype *LLR, const P &bp){
	const dim3 grid(((int(C)-1)/256+1),((int(VW)-1)/4+1),1);
	const dim3 thread(256,4,1);
	//apply LLR
	parallel_broadcast_add<<<grid,thread>>>(alphabeta.get(), LLR, C, VW);
	//row update
	bp.forward_vec(alphabeta.get(), Hones);
	parallel_accumlate<<<grid,thread>>>(alphabeta.get(), alphabetaidx.get(), W, Hsize, bp);
	bp.backward_vec(alphabeta.get(), Hones);
	//column update
	cudaMemset(LPR, 0, sizeof(fptype)*C);
	parallel_reduce_add<<<grid,thread>>>(LPR, alphabeta.get(), C, VW);
	parallel_broadcast_negsub<<<grid,thread>>>(alphabeta.get(), LPR, C, VW);
	parallel_broadcast_add<<<grid,thread>>>(LPR, LLR, C, 1);
	//parity check
	// for(std::size_t i=0; i<Hsize; ++i){
	// 	auto parity = false;
	// 	for(const auto &j : H[i]) parity ^= LPR[j]<0;
	// 	if(parity) return false;
	// }
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
template<boxplusclass P>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, std::uint32_t iterationlimit){
	auto alphabeta = util::make_cuda_unique<fptype[]>(alphabetasize());

	fptype *LPR_device;//対数事後確率比：列ごとのalphaの和+QLLR
	fptype *LLR_device;
	auto errc = ::cudaMalloc(&LPR_device, sizeof(fptype)*C);
	if(errc!=0) throw std::runtime_error("CUDA Error");
	errc = ::cudaMalloc(&LLR_device, sizeof(fptype)*C);
	if(errc!=0) throw std::runtime_error("CUDA Error");
	errc = ::cudaMemcpy(LLR_device,LLR.data(),sizeof(fptype)*C,cudaMemcpyHostToDevice);
	if(errc!=0) throw std::runtime_error("CUDA Error");

	int itr = 0;
	bool *parity = nullptr;
	while(itr<iterationlimit){
		iterate(parity, alphabeta, LPR_device, LLR_device, bp);
		++itr;
	}

	errc = ::cudaMemcpy(LPR.data(),LPR_device,sizeof(fptype)*C,cudaMemcpyDeviceToHost);
	if(errc!=0) throw std::runtime_error("CUDA Error");

	errc = ::cudaFree(LPR_device);
	if(errc!=0) throw std::runtime_error("CUDA Error");
	errc = ::cudaFree(LLR_device);
	if(errc!=0) throw std::runtime_error("CUDA Error");
}

}

#endif