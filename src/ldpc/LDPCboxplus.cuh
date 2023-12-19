#ifndef INCLUDE_GUARD_ldpc_LDPCboxplus
#define INCLUDE_GUARD_ldpc_LDPCboxplus

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cstdint>
#include <bit>
#include <concepts>
#include <cmath>
#include <limits>
#include <cassert>

namespace code::LDPC {

namespace{
	template<std::floating_point T>
	struct uint_of_length {};
	template<>
	struct uint_of_length<float> {
		using type = std::uint32_t;
	};
	template<>
	struct uint_of_length<double> {
		using type = std::uint64_t;
	};
	template<std::floating_point T>
	using uint_of_length_t = uint_of_length<T>::type;

	template<std::floating_point T>
	class sum_accumlator{
		T abssum;
		uint_of_length_t<T> signprod;
	public:
		__host__ __device__ sum_accumlator():abssum(0),signprod(0){}
		__host__ __device__ void operator+=(T rhs);
		__host__ __device__ void operator+=(sum_accumlator<T> rhs);
		__host__ __device__ T operator-(T rhs) const;
	};

	template<std::floating_point T>
	class min_accumlator{
		std::pair<T,T> absmin;
		uint_of_length_t<T> signprod;
	public:
		__host__ __device__ min_accumlator():absmin(std::numeric_limits<T>::infinity(),std::numeric_limits<T>::infinity()),signprod(0){}
		__host__ __device__ void operator+=(T rhs);
		__host__ __device__ T operator-(T rhs) const;
	};
}

struct minsum {
	template<std::floating_point T>
	static inline T forward(T x){return x;}
	template<std::floating_point T>
	static inline T backward(T x){return x;}
};

template<std::uint8_t precision=0>
class phi_calc {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto UPPER_BOUND = 0x1p6f;

	template<std::floating_point T>
	static T common(T x);
public:
	template<std::floating_point T>
	static inline T forward(T x){return common(x);}
	template<std::floating_point T>
	static inline T backward(T x){return common(x);}
};

template<std::uint8_t precision=0> class phi_table;

template<>
class phi_table<0> {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto LOWER_BOUND_U = std::bit_cast<uint32_t>(LOWER_BOUND);
	static constexpr auto UPPER_BOUND = 0x1p6f;
	static constexpr auto UPPER_BOUND_U = std::bit_cast<uint32_t>(UPPER_BOUND);
	static constexpr auto VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	static constexpr auto CACHE_FILENAME = "gallager_float.bin";

	inline static std::unique_ptr<float[]> values;
	inline static std::unique_ptr<float[],util::cuda_delete> values_device;

	static void values_init();//キャッシュファイルを読み込む。失敗したら、値を計算してキャッシュファイルに保存する。
	static bool read_values(std::unique_ptr<float[]> &vec);
	static bool write_values(const std::unique_ptr<float[]> &vec);
	inline float get(float x) const;
	__global__ friend void get_parallel(const phi_table<0> &p ,float *arr, std::uint32_t size, float *values);
public:
	phi_table();
	template<std::floating_point T>
	inline T forward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	template<std::floating_point T>
	inline T backward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	inline void forward_vec(float *arr, std::uint32_t size) const{get_parallel<<<((size-1)/1024+1),1024>>>(*this,arr,size,values_device.get());}
	inline void backward_vec(float *arr, std::uint32_t size) const{get_parallel<<<((size-1)/1024+1),1024>>>(*this,arr,size,values_device.get());}
};

template<>
class phi_table<1> {
	static constexpr auto LOWER_BOUND = 0x1p-10f;
	static constexpr auto LOWER_BOUND_U = 0x00000000u;
	static constexpr auto UPPER_BOUND = 0x1.ffdp4f;
	static constexpr auto UPPER_BOUND_U = 0x0000ffffu;
	static constexpr auto EXPONENT_BIAS = (11+0x00000080)<<23;//指数の内部表現の差
	static constexpr auto SHIFT_HALF_FLOAT = 11u;//仮数部の長さの差
	static constexpr auto VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	static constexpr auto CACHE_FILENAME = "gallager_half.bin";

	inline static std::unique_ptr<float[]> values;//インデックスは符号なし浮動小数点数[指数4bits(-10~+5)+仮数12bits(ケチ表現)]の内部表現

	static void values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	static bool read_values(std::unique_ptr<float[]> &vec);
	static bool write_values(const std::unique_ptr<float[]> &vec);
	inline float get(float x) const;
public:
	phi_table();
	template<std::floating_point T>
	inline T forward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	template<std::floating_point T>
	inline T backward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
};

template<>
class phi_table<2> {
	static constexpr auto LOWER_BOUND = 0x1p-10f;
	static constexpr auto LOWER_BOUND_U = 0x00000000u;
	static constexpr auto UPPER_BOUND = 0x1.ffdp4f;
	static constexpr auto UPPER_BOUND_U = 0x000000ffu;
	static constexpr auto EXPONENT_BIAS = (11+0x00000080)<<23;//指数の内部表現の差
	static constexpr auto SHIFT_MINI_FLOAT = 19u;//仮数部の長さの差
	static constexpr auto VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;

	inline static std::unique_ptr<float[]> values;//インデックスは符号なし浮動小数点数[指数4bits(-10~+5)+仮数12bits(ケチ表現)]の内部表現

	static void values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	inline float get(float x) const;
public:
	phi_table();
	template<std::floating_point T>
	inline T forward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	template<std::floating_point T>
	inline T backward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
};

template<class B, std::floating_point T> struct accumlator;

template<std::floating_point T>
struct accumlator<minsum,T> {
	using type = min_accumlator<T>;
};
template<std::floating_point T, std::uint8_t P>
struct accumlator<phi_calc<P>,T> {
	using type = sum_accumlator<T>;
};
template<std::floating_point T, std::uint8_t P>
struct accumlator<phi_table<P>,T> {
	using type = sum_accumlator<T>;
};

template<class B, std::floating_point T> using accumlator_t = accumlator<B,T>::type;

template<class T>
concept boxplusclass = requires(T x){
	typename accumlator_t<T,float>;
	{x.forward(0.f)} -> std::same_as<float>;
	{x.backward(0.f)} -> std::same_as<float>;
};

////////////////////////////////////////////////////////////////
//                                                            //
//                    class sum_accumlator                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::floating_point T>
__host__ __device__ inline void sum_accumlator<T>::operator+=(T rhs){
	abssum += std::fabs(rhs);
	#ifdef __CUDA_ARCH__
		signprod ^= reinterpret_cast<uint_of_length_t<T>&>(rhs);
	#else
		signprod ^= std::bit_cast<uint_of_length_t<T>>(rhs);
	#endif
}

template<std::floating_point T>
__host__ __device__ inline void sum_accumlator<T>::operator+=(sum_accumlator<T> rhs){
	abssum += rhs.abssum;
	signprod ^= rhs.signprod;
}

template<std::floating_point T>
__host__ __device__ inline T sum_accumlator<T>::operator-(T rhs) const{
	constexpr uint_of_length_t<T> signmask = 1u<<(sizeof(T)*8u-1u);
	#ifdef __CUDA_ARCH__
		auto diff = abssum-std::fabs(rhs);
		auto temp = (signprod^reinterpret_cast<uint_of_length_t<T>&>(rhs))&signmask|reinterpret_cast<uint_of_length_t<T>&>(diff);
		return reinterpret_cast<T&>(temp);
	#else
		return std::bit_cast<T>((signprod^std::bit_cast<uint_of_length_t<T>>(rhs))&signmask|std::bit_cast<uint_of_length_t<T>>(abssum-std::fabs(rhs)));
	#endif
}

////////////////////////////////////////////////////////////////
//                                                            //
//                    class min_accumlator                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::floating_point T>
__host__ __device__ inline void min_accumlator<T>::operator+=(T rhs){
	auto arhs = std::fabs(rhs);
	if(absmin.first>=arhs){
		absmin.second = absmin.first;
		absmin.first = arhs;
	}else if(absmin.second>arhs){
		absmin.second = arhs;
	}
	signprod ^= std::bit_cast<uint_of_length_t<T>>(rhs);
}

template<std::floating_point T>
__host__ __device__ inline T min_accumlator<T>::operator-(T rhs) const{
	constexpr uint_of_length_t<T> signmask = 1u<<(sizeof(T)*8u-1u);
	auto arhs = std::fabs(rhs);
	auto val = absmin.first==arhs?absmin.second:absmin.first;
	return std::bit_cast<T>((signprod^std::bit_cast<uint_of_length_t<T>>(rhs))&signmask|std::bit_cast<uint_of_length_t<T>>(val));
}

////////////////////////////////////////////////////////////////
//                                                            //
//                       class phi_calc                       //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint8_t precision>
template<std::floating_point T>
inline T phi_calc<precision>::common(T x){
	static_assert(precision<2,"invalid precision specification");
	auto y = std::fabs(x);
	//定義域を限定
	if(y<LOWER_BOUND) y = LOWER_BOUND;
	if(y>UPPER_BOUND) y = UPPER_BOUND;
	if constexpr(precision==0) y = static_cast<T>(std::log1p(static_cast<double>(2)/std::expm1(static_cast<double>(y))));
	if constexpr(precision==1) y = std::log1p(static_cast<T>(2)/std::expm1(y));
	y = std::bit_cast<T>(std::bit_cast<uint_of_length_t<T>>(y)|std::bit_cast<uint_of_length_t<T>>(x)&0x80000000);
	return y;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                     class phi_table<0>                     //
//                                                            //
////////////////////////////////////////////////////////////////

void phi_table<0>::values_init(){
	values = std::make_unique<float[]>(VALUE_RANGE);

	if(!read_values(values)){
		std::cerr<<"phi_table<0>: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<VALUE_RANGE; ++i){
			auto iu = i+LOWER_BOUND_U;
			values[i] = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(std::bit_cast<float>(iu)))));
		}
		if(!write_values(values)) std::cerr<<"phi_table<0>: Caching failed."<<std::endl;
	}

	values_device = util::make_cuda_unique<float[]>(VALUE_RANGE);
	auto errc = cudaMemcpy(values_device.get(), values.get(), sizeof(float)*VALUE_RANGE, cudaMemcpyHostToDevice);
	if(errc!=0) throw std::runtime_error("CUDA Error");
	cudaDeviceSynchronize();
}

bool phi_table<0>::read_values(std::unique_ptr<float[]> &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.get()), VALUE_RANGE*sizeof(float));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool phi_table<0>::write_values(const std::unique_ptr<float[]> &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.get()), VALUE_RANGE*sizeof(float));
	file.close();
	if(file.fail()) return false;
	return true;
}

phi_table<0>::phi_table(){
	if(!values) values_init();
}

inline float phi_table<0>::get(float x) const{
	auto xa = std::fabs(x);
	//定義域を限定
	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
	auto ya = values[std::bit_cast<std::uint32_t>(xa) - LOWER_BOUND_U];
	return std::bit_cast<float>(std::bit_cast<std::uint32_t>(ya)|std::bit_cast<std::uint32_t>(x)&0x80000000);
}

__global__ void get_parallel(const phi_table<0> &p ,float *arr, std::uint32_t size, float *values){
	std::uint32_t i = blockIdx.x*blockDim.x+threadIdx.x;
	if(i<size){
		auto x = arr[i];
		auto xa = std::fabs(x);
		//定義域を限定
		if(xa<phi_table<0>::LOWER_BOUND) xa = phi_table<0>::LOWER_BOUND;
		if(xa>phi_table<0>::UPPER_BOUND) xa = phi_table<0>::UPPER_BOUND;
		assert(reinterpret_cast<std::uint32_t&>(xa) - phi_table<0>::LOWER_BOUND_U < phi_table<0>::VALUE_RANGE);
		auto ya = values[reinterpret_cast<std::uint32_t&>(xa) - phi_table<0>::LOWER_BOUND_U];
		auto yu = reinterpret_cast<std::uint32_t&>(ya)|reinterpret_cast<std::uint32_t&>(x)&0x80000000;
		arr[i] = reinterpret_cast<float&>(yu);
	}
}

////////////////////////////////////////////////////////////////
//                                                            //
//                     class phi_table<1>                     //
//                                                            //
////////////////////////////////////////////////////////////////

void phi_table<1>::values_init(){
	values = std::make_unique<float[]>(VALUE_RANGE);

	if(!read_values(values)){
		std::cerr<<"phi_table<1>: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<VALUE_RANGE; ++i){
			auto x = std::bit_cast<float>(((i<<SHIFT_HALF_FLOAT)-EXPONENT_BIAS)&0x7fffffff);
			auto y = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(x))));
			values[i] = y;
		}
		if(!write_values(values)) std::cerr<<"phi_table<1>: Caching failed."<<std::endl;
	}
}

bool phi_table<1>::read_values(std::unique_ptr<float[]> &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.get()), VALUE_RANGE*sizeof(float));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool phi_table<1>::write_values(const std::unique_ptr<float[]> &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.get()), VALUE_RANGE*sizeof(float));
	file.close();
	if(file.fail()) return false;
	return true;
}

phi_table<1>::phi_table(){
	if(!values) values_init();
}

inline float phi_table<1>::get(float x) const{
	auto xa = std::fabs(x);
	//定義域を限定
	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
	auto xu = ((std::bit_cast<std::uint32_t>(xa)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT)&UPPER_BOUND_U;
	auto ya = values[xu];
	//線形補間
	constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_HALF_FLOAT);
	constexpr float ratiounit = 1.0f/(1<<SHIFT_HALF_FLOAT);
	auto y2 = values[xu+1];
	ya += (y2-ya)*static_cast<float>(std::bit_cast<std::uint32_t>(xa)&bottommask)*ratiounit;

	return std::bit_cast<float>(std::bit_cast<uint32_t>(ya)|std::bit_cast<uint32_t>(x)&0x80000000);
}

// inline float phi_table<1>::get(float x) const{
// 	auto xa = std::fabs(x);
// 	//定義域を限定
// 	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
// 	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
// 	#ifdef __CUDA_ARCH__
// 		auto xu = ((reinterpret_cast<uint32_t&>(xa)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT)&UPPER_BOUND_U;
// 		auto ya = values_device[xu];
// 		//線形補間
// 		constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_HALF_FLOAT);
// 		constexpr float ratiounit = 1.0f/(1<<SHIFT_HALF_FLOAT);
// 		auto y2 = values_device[xu+1];
// 		ya += (y2-ya)*static_cast<float>(reinterpret_cast<uint32_t&>(xa)&bottommask)*ratiounit;

// 		auto yu = reinterpret_cast<uint32_t&>(ya)|reinterpret_cast<uint32_t&>(x)&0x80000000;
// 		return reinterpret_cast<float&>(yu);
// 	#else
// 		auto xu = ((std::bit_cast<std::uint32_t>(xa)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT)&UPPER_BOUND_U;
// 		auto ya = values[xu];
// 		//線形補間
// 		constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_HALF_FLOAT);
// 		constexpr float ratiounit = 1.0f/(1<<SHIFT_HALF_FLOAT);
// 		auto y2 = values[xu+1];
// 		ya += (y2-ya)*static_cast<float>(std::bit_cast<std::uint32_t>(xa)&bottommask)*ratiounit;

// 		return std::bit_cast<float>(std::bit_cast<uint32_t>(ya)|std::bit_cast<uint32_t>(x)&0x80000000);
// 	#endif
// }

////////////////////////////////////////////////////////////////
//                                                            //
//                     class phi_table<2>                     //
//                                                            //
////////////////////////////////////////////////////////////////

void phi_table<2>::values_init(){
	values = std::make_unique<float[]>(VALUE_RANGE);

	for(auto i=0ui32; i<VALUE_RANGE; ++i){
		auto x = std::bit_cast<float>(((i<<SHIFT_MINI_FLOAT)-EXPONENT_BIAS)&0x7fffffff);
		auto y = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(x))));
		values[i] = y;
	}
}

phi_table<2>::phi_table(){
	if(!values) values_init();
}

inline float phi_table<2>::get(float x) const{
	auto xa = std::fabs(x);
	//定義域を限定
	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
	auto xu = ((std::bit_cast<std::uint32_t>(xa)+EXPONENT_BIAS)>>SHIFT_MINI_FLOAT)&UPPER_BOUND_U;
	auto ya = values[xu];
	//線形補間
	// constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_MINI_FLOAT);
	// constexpr float ratiounit = 1.0f/(1<<SHIFT_MINI_FLOAT);
	// auto y2 = values[xu+1];
	// ya += (y2-ya)*static_cast<float>(std::bit_cast<std::uint32_t>(xa)&bottommask)*ratiounit;
	//符号を復元
	return std::bit_cast<float>(std::bit_cast<std::uint32_t>(ya)|std::bit_cast<std::uint32_t>(x)&0x80000000);
}

// inline float phi_table<2>::get(float x) const{
// 	auto xa = std::fabs(x);
// 	//定義域を限定
// 	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
// 	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
// 	#ifdef __CUDA_ARCH__
// 		auto xu = ((reinterpret_cast<uint32_t&>(xa)+EXPONENT_BIAS)>>SHIFT_MINI_FLOAT)&UPPER_BOUND_U;
// 		auto ya = values_device[xu];
// 		//線形補間
// 		// constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_MINI_FLOAT);
// 		// constexpr float ratiounit = 1.0f/(1<<SHIFT_MINI_FLOAT);
// 		// auto y2 = values_device[xu+1];
// 		// ya += (y2-ya)*static_cast<float>(reinterpret_cast<uint32_t&>(xa)&bottommask)*ratiounit;

// 		auto yu = reinterpret_cast<uint32_t&>(ya)|reinterpret_cast<uint32_t&>(x)&0x80000000;
// 		return reinterpret_cast<float&>(yu);
// 	#else
// 		auto xu = ((std::bit_cast<std::uint32_t>(xa)+EXPONENT_BIAS)>>SHIFT_MINI_FLOAT)&UPPER_BOUND_U;
// 		auto ya = values[xu];
// 		//線形補間
// 		// constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_MINI_FLOAT);
// 		// constexpr float ratiounit = 1.0f/(1<<SHIFT_MINI_FLOAT);
// 		// auto y2 = values[xu+1];
// 		// ya += (y2-ya)*static_cast<float>(std::bit_cast<std::uint32_t>(xa)&bottommask)*ratiounit;

// 		return std::bit_cast<float>(std::bit_cast<uint32_t>(ya)|std::bit_cast<uint32_t>(x)&0x80000000);
// 	#endif
// }


}

#endif