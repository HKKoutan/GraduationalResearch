#ifndef INCLUDE_GUARD_common_memory
#define INCLUDE_GUARD_common_memory

#include <memory>
#include <exception>

namespace util{

void check_cuda_error(cudaError_t errc){
	if(errc != cudaError::cudaSuccess) throw std::runtime_error("CUDA Error");
}

struct cuda_delete {
	void operator()(void* p) const{check_cuda_error(::cudaFree(p));}
};

template<typename T>
typename std::enable_if<std::is_array<T>::value, std::unique_ptr<T,cuda_delete>>::type make_cuda_unique(const std::size_t n){
	using U = std::remove_extent<T>::type;
	U* p = nullptr;
	check_cuda_error(::cudaMalloc(reinterpret_cast<void**>(&p), sizeof(U)*n));
	return std::unique_ptr<T,cuda_delete>{p};
}

template<typename T>
typename std::unique_ptr<T,cuda_delete>::type make_cuda_unique(){
	T* p = nullptr;
	check_cuda_error(::cudaMalloc(reinterpret_cast<void**>(&p), sizeof(T)));
	return std::unique_ptr<T,cuda_delete>{p};
}

}
#endif