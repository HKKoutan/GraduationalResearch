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
	// using fptype = float;
	static constexpr std::size_t S = T::sourcesize();
	static constexpr std::size_t C = T::codesize();
	static constexpr std::size_t Hsize = C-S;

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
	using T = CheckMatrix_regular<S,C,W>;
	// using fptype = float;
	static constexpr std::size_t Hsize = C-S;
	static constexpr std::size_t Hones = W*Hsize;
	static constexpr std::size_t VW = Hones/C;//列重み
	static_assert(Hones<(1ui64<<32));

	const T H;//検査行列
	std::unique_ptr<fptype[][C],util::cuda_delete<fptype[][C]>> alphabeta;
	std::unique_ptr<uitype[][W],util::cuda_delete<uitype[][W]>> alphabetaidx;
	// std::unique_ptr<fptype*[][W],util::cuda_delete<fptype*[][W]>> alphabetap;
	// std::array<std::array<fptype,C>,VW> alphabeta;
	// std::array<std::array<fptype*,W>,Hsize> alphabetap;
public:
	explicit Sumproduct_decoding(const T &H);
	void decode_init();//decodeで使用する変数の初期化
	template<boxplusclass P>
	void iterate(bool *parity, fptype *LPR, fptype *LLR, const P &bp);
	template<boxplusclass P>
	void decode(fptype *LPR, fptype *LLR, const P &bp, int iterationlimit);
};

////////////////////////////////////////////////////////////////
//                                                            //
//                  class Sumproduct_decoding                 //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
auto Sumproduct_decoding<T>::alphabeta_size(const T &H){
	std::array<std::size_t,C> Hheight{};
	for(std::size_t i=0; i<C; ++i) Hheight[i] = H.colsize(i);
	return std::ranges::max(Hheight);
}

template<CheckMatrix T>
auto Sumproduct_decoding<T>::alphabetap_init(const T &H, std::vector<std::pair<std::array<fptype,C>,std::array<fptype,C>>> &alphabeta){
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

// template<std::size_t S, std::size_t C, std::size_t W>
// void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::alphabetap_init(){
// 	// std::array<std::array<std::size_t,VW>,C> HT;//Hの転置
// 	// {
// 	// 	std::array<std::size_t,C> HTc = {};
// 	// 	for(std::size_t i=0; i<Hsize; ++i) for(auto j: H[i]) HT[j][HTc[j]++] = i;
// 	// }
// 	for(std::size_t i=0; i<Hsize; ++i){
// 		auto &Hi = H[i];
// 		auto &abpi = alphabetap[i];
// 		//alpha<-alphap beta<-betap
// 		for(std::size_t j=0; j<W; ++j){
// 			auto &hij = Hi[j];
// 			auto &Hj = H.T[hij];
// 			std::size_t k=0;
// 			while(Hj[k]!=i) ++k;
// 			auto &abk = alphabeta[k];
// 			abpi[j] = &abk[hij];
// 		}
// 	}
// }

namespace {
	template<CheckMatrix U>
	__global__ void alphabetap_init(uitype *alphabetaidx, U H, typename U::internaldatatype Hd){
		int i = blockIdx.x*blockDim.x+threadIdx.x;
		int j = blockIdx.y*blockDim.y+threadIdx.y;
		if(i<H.size()){
			int W = H.colweight(i);
			if(j<W){
				auto Hi = U::getrow(i,Hd);
				auto hij = Hi[j];
				auto Hj = U::getcol(hij,Hd);
				int k=0;
				while(Hj[k]!=i) ++k;
				alphabetaidx[W*i+j] = H.codesize()*k+hij;
			}
		}
	}
}

template<std::size_t S, std::size_t C, std::size_t W>
Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::Sumproduct_decoding(const T &H):
	H(H),
	alphabeta(util::make_cuda_unique<fptype[][C]>(VW)),
	alphabetaidx(util::make_cuda_unique<uitype[][W]>(Hsize))
{
	const dim3 grid(((C-1)/128+1),((VW-1)/8+1),1);
	const dim3 thread(128,8,1);
	alphabetap_init<<<grid,thread>>>(&alphabetaidx[0][0], H, H.data());
}

template<std::size_t S, std::size_t C, std::size_t W>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::decode_init(){
	cudaMemset(&alphabeta[0][0], 0, sizeof(fptype)*Hones);
}

namespace {
	__global__ void parallel_broadcast_add(fptype *lhs, fptype *rhs, std::size_t width, std::size_t height){//lhs[i][j]+=rhs[j] for all i<height, j<width
		int j = blockIdx.x*blockDim.x+threadIdx.x;
		int i = blockIdx.y*blockDim.y+threadIdx.y;
		if(i<height&&j<width) lhs[i*width+j] += rhs[j];
	}
	__global__ void parallel_broadcast_negsub(fptype *lhs, fptype *rhs, std::size_t width, std::size_t height){//lhs[i][j]=rhs[j]-lhs[i][j] for all i<height, j<width
		int j = blockIdx.x*blockDim.x+threadIdx.x;
		int i = blockIdx.y*blockDim.y+threadIdx.y;
		if(i<height&&j<width) lhs[i*width+j] = rhs[j]-lhs[i*width+j];
	}
	__global__ void parallel_reduce_add(fptype *lhs, fptype *rhs, std::size_t width, std::size_t height){//lhs[j]+=rhs[i][j] for all i<height, j<width
		int j = blockIdx.x*blockDim.x+threadIdx.x;
		int k = blockIdx.y*blockDim.y+threadIdx.y;
		if(k==0&&j<width){
			for(int i=0; i<height; ++i) lhs[j] += rhs[i*width+j];
		}
	}
	template<class P>
	__global__ void parallel_accumlate(fptype *arr, uitype *idx, std::size_t width, std::size_t height, const P &bp){
		int j = blockIdx.x*blockDim.x+threadIdx.x;
		int k = blockIdx.y*blockDim.y+threadIdx.y;
		if(k==0&&j<height){
			accumlator_t<P,fptype> acc;
			for(int i=0; i<width; ++i) acc += arr[idx[j*width+i]];
			for(int i=0; i<width; ++i) arr[idx[j*width+i]] = acc-arr[idx[j*width+i]];
		}
	}
	// template<typename T,CheckMatrix U>
	// __global__ void check_parity(T* ptr, U H, typename U::internaldatatype Hd){
	// 	int j = blockIdx.x*blockDim.x+threadIdx.x;
	// 	int k = blockIdx.y*blockDim.y+threadIdx.y;
	// }
}

template<std::size_t S, std::size_t C, std::size_t W>
template<boxplusclass P>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::iterate(bool *parity, fptype *LPR, fptype *LLR, const P &bp){
	const dim3 grid(((int(C)-1)/256+1),((int(VW)-1)/4+1),1);
	const dim3 thread(256,4,1);
	//apply LLR
	// for(auto i=0; i<VW; ++i){
	// 	auto &bi = alphabeta[i];
	// 	for(std::size_t j=0; j<C; ++j) bi[j] += LLR[j];
	// }
	parallel_broadcast_add<<<grid,thread>>>(&alphabeta[0][0], LLR, C, VW);
	//row update
	// for(auto i=0; i<VW; ++i) for(auto &bij: alphabeta[i]) bij = bp.forward(bij);
	bp.forward_vec(&alphabeta[0][0], Hones);
	// for(auto i=0; i<Hsize; ++i){
	// 	auto &abpi = alphabetap[i];
	// 	accumlator_t<P,fptype> acc;
	// 	for(const auto bpij: abpi) acc += *bpij;
	// 	for(const auto abpij: abpi) *abpij = acc-*abpij;
	// }
	parallel_accumlate<<<grid,thread>>>(&alphabeta[0][0], &alphabetaidx[0][0], W, Hsize, bp);
	// for(auto i=0; i<VW; ++i) for(auto &aij: alphabeta[i]) aij = bp.backward(aij);
	bp.backward_vec(&alphabeta[0][0], Hones);
	//column update
	// for(std::size_t j=0; j<C; ++j) LPR[j] = 0;
	cudaMemset(LPR, 0, sizeof(fptype)*C);
	// for(auto i=0; i<VW; ++i){
	// 	auto &ai = alphabeta[i];
	// 	for(std::size_t j=0; j<C; ++j) LPR[j] += ai[j];
	// }
	parallel_reduce_add<<<grid,thread>>>(LPR, &alphabeta[0][0], C, VW);
	// for(auto i=0; i<VW; ++i){
	// 	auto &abi = alphabeta[i];
	// 	for(std::size_t j=0; j<C; ++j) abi[j] = LPR[j]-abi[j];
	// }
	parallel_broadcast_negsub<<<grid,thread>>>(&alphabeta[0][0], LPR, C, VW);
	// for(std::size_t j=0; j<C; ++j) LPR[j] += LLR[j];
	parallel_broadcast_add<<<grid,thread>>>(LPR, LLR, C, 1);
	//parity check
	// for(std::size_t i=0; i<Hsize; ++i){
	// 	auto parity = false;
	// 	for(const auto &j : H[i]) parity ^= LPR[j]<0;
	// 	if(parity) return false;
	// }
}

template<std::size_t S, std::size_t C, std::size_t W>
template<boxplusclass P>
void Sumproduct_decoding<CheckMatrix_regular<S,C,W>>::decode(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR, const P &bp, int iterationlimit){
	int itr = 0;
	bool *parity = nullptr;
	fptype *LPR_device = nullptr;//対数事後確率比：列ごとのalphaの和+QLLR
	fptype *LLR_device = nullptr;
	cudaMalloc(&LPR_device, sizeof(fptype)*C);
	cudaMalloc(&LLR_device, sizeof(fptype)*C);
	cudaMemcpy(LLR_device,LLR.data(),sizeof(fptype)*C,cudaMemcpyHostToDevice);
	// cudaMalloc(&parity, sizeof(bool));
	// cudaMemset(parity,0,sizeof(bool));
	while(itr<iterationlimit){
		iterate(parity, LPR_device, LLR_device, bp);
		++itr;
	}

	cudaMemcpy(LPR.data(),LPR_device,sizeof(fptype)*C,cudaMemcpyDeviceToHost);
	cudaFree(LPR_device);
	cudaFree(LLR_device);
}

}

#endif