#include <memory>
#include <exception>

namespace util{

template<typename T>
struct cuda_delete {
	void operator()(void* p) const{
		auto errc = ::cudaFree(p);
		if(errc!=0) throw std::runtime_error("CUDA Error");
	}
};

template<typename T>
typename std::enable_if<std::is_array<T>::value, std::unique_ptr<T,cuda_delete<T>>>::type make_cuda_unique(const std::size_t n){
	using U = std::remove_extent<T>::type;
	U* p = nullptr;
	auto errc = ::cudaMalloc(reinterpret_cast<void**>(&p), sizeof(U)*n);
	if(errc!=0) throw std::runtime_error("CUDA Error");
	return std::unique_ptr<T,cuda_delete<T>>{p};
}

template<typename T>
typename std::unique_ptr<T,cuda_delete<T>>::type make_cuda_unique(){
	T* p = nullptr;
	auto errc = ::cudaMalloc(reinterpret_cast<void**>(&p), sizeof(T));
	if(errc!=0) throw std::runtime_error("CUDA Error");
	return std::unique_ptr<T,cuda_delete<T>>{p};
}

}
