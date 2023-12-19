#include <memory>
#include <exception>

namespace util{

struct cuda_delete {
	void operator()(void* p) const{
		auto errc = ::cudaFree(p);
		if(errc!=0) throw std::runtime_error("CUDA Error");
	}
};

template<typename T>
typename std::enable_if<std::is_array<T>::value, std::unique_ptr<T,cuda_delete>>::type make_cuda_unique(const std::size_t n){
	using U = std::remove_extent<T>::type;
	U* p = nullptr;
	::cudaMalloc(reinterpret_cast<void**>(&p), sizeof(U)*n);
	return std::unique_ptr<T,cuda_delete>{p};
}

template<typename T>
typename std::unique_ptr<T,cuda_delete>::type make_cuda_unique(){
	T* p = nullptr;
	::cudaMalloc(reinterpret_cast<void**>(&p), sizeof(T));
	return std::unique_ptr<T,cuda_delete>{p};
}

}
