#ifndef INCLUDE_GUARD_ldpc_LDPCboxplus
#define INCLUDE_GUARD_ldpc_LDPCboxplus

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <bit>
#include <concepts>
#include <cmath>
#include <limits>

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

	template<typename T>
	class device_table {
		std::vector<T> values;
		T *values_device = nullptr;
	public:
		device_table() = default;
		~device_table();
		T* devicedata(){return values_device;}
		operator bool(){return values.size();}
		void operator=(std::vector<T> &&rhs);
		inline T operator[](std::size_t i);
	};
}

struct minsum {
	template<std::floating_point T>
	__host__ __device__ static inline T forward(T x){return x;}
	template<std::floating_point T>
	__host__ __device__ static inline T backward(T x){return x;}
};

template<std::uint8_t precision=0>
class phi_calc {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto UPPER_BOUND = 0x1p6f;

	template<std::floating_point T>
	__host__ __device__ static T common(T x);
public:
	template<std::floating_point T>
	__host__ __device__ static inline T forward(T x){return common(x);}
	template<std::floating_point T>
	__host__ __device__ static inline T backward(T x){return common(x);}
};

template<std::uint8_t precision=0> class phi_table;

template<>
class phi_table<0> {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto LOWER_BOUND_U = std::bit_cast<uint32_t>(LOWER_BOUND);
	static constexpr auto UPPER_BOUND = 0x1p6f;
	static constexpr auto UPPER_BOUND_U = std::bit_cast<uint32_t>(UPPER_BOUND);
	static constexpr auto CACHE_FILENAME = "gallager_float.bin";

	inline static device_table<float> values;
	float *values_device;

	static std::vector<float> values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	static bool read_values(std::vector<float> &vec);
	static bool write_values(const std::vector<float> &vec);
	__host__ __device__ float get(float x) const;
public:
	phi_table();
	template<std::floating_point T>
	__host__ __device__ inline T forward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	template<std::floating_point T>
	__host__ __device__ inline T backward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
};

template<>
class phi_table<1> {
	static constexpr auto LOWER_BOUND = 0x1p-10f;
	static constexpr auto LOWER_BOUND_U = 0x00000000;
	static constexpr auto UPPER_BOUND = 0x1.ffdp4f;
	static constexpr auto UPPER_BOUND_U = 0x0000ffff;
	static constexpr auto EXPONENT_BIAS = (11+0x00000080)<<23;//指数の内部表現の差
	static constexpr auto SHIFT_HALF_FLOAT = 11;//仮数部の長さの差
	static constexpr auto CACHE_FILENAME = "gallager_half.bin";

	inline static device_table<float> values;//インデックスは符号なし浮動小数点数[指数4bits(-10~+5)+仮数12bits(ケチ表現)]の内部表現
	float *values_device;

	static std::vector<float> values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	static bool read_values(std::vector<float> &vec);
	static bool write_values(const std::vector<float> &vec);
	__host__ __device__ float get(float x) const;
public:
	phi_table();
	template<std::floating_point T>
	__host__ __device__ inline T forward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	template<std::floating_point T>
	__host__ __device__ inline T backward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
};

template<>
class phi_table<2> {
	static constexpr auto LOWER_BOUND = 0x1p-10f;
	static constexpr auto LOWER_BOUND_U = 0x00000000;
	static constexpr auto UPPER_BOUND = 0x1.ffdp4f;
	static constexpr auto UPPER_BOUND_U = 0x000000ff;
	static constexpr auto EXPONENT_BIAS = (11+0x00000080)<<23;//指数の内部表現の差
	static constexpr auto SHIFT_MINI_FLOAT = 19;//仮数部の長さの差
	static constexpr auto FG_VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;

	inline static device_table<float> values;//インデックスは符号なし浮動小数点数[指数4bits(-10~+5)+仮数4bits(ケチ表現)]の内部表現
	float *values_device;

	static std::vector<float> values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	__host__ __device__ float get(float x) const;
public:
	phi_table();
	template<std::floating_point T>
	__host__ __device__ inline T forward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	template<std::floating_point T>
	__host__ __device__ inline T backward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
};

template<class B, std::floating_point T> struct accumlator;

template<std::floating_point T>
struct accumlator<minsum,T> {
	using type = min_accumlator<T>;
	// using dtype = min_accumlator<T>;
};
template<std::floating_point T, std::uint8_t P>
struct accumlator<phi_calc<P>,T> {
	using type = sum_accumlator<T>;
	// using dtype = sum_accumlator<T>;
};
template<std::floating_point T, std::uint8_t P>
struct accumlator<phi_table<P>,T> {
	using type = sum_accumlator<T>;
	// using dtype = sum_accumlator<T>;
};

template<class B, std::floating_point T> using accumlator_t = accumlator<B,T>::type;
// template<class B, std::floating_point T> using deviceaccumlator_t = accumlator<B,T>::dtype;

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
	signprod ^= std::bit_cast<uint_of_length_t<T>>(rhs);
}

template<std::floating_point T>
__host__ __device__ inline T sum_accumlator<T>::operator-(T rhs) const{
	constexpr uint_of_length_t<T> signmask = 1u<<(sizeof(T)*8u-1u);
	return std::bit_cast<T>((signprod^std::bit_cast<uint_of_length_t<T>>(rhs))&signmask|std::bit_cast<uint_of_length_t<T>>(abssum-std::fabs(rhs)));
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
//                     class device_table                     //
//                                                            //
////////////////////////////////////////////////////////////////

template<typename T>
device_table<T>::~device_table(){
	cudaFree(values_device);
}

template<typename T>
void device_table<T>::operator=(std::vector<T> &&rhs){
	values = rhs;
	cudaMalloc(&values_device, sizeof(T)*values.size());
	cudaMemcpy(values_device, values.data(), sizeof(T)*values.size(), cudaMemcpyHostToDevice);
}

template<typename T>
inline T device_table<T>::operator[](std::size_t i){
	return values[i];
}

////////////////////////////////////////////////////////////////
//                                                            //
//                       class phi_calc                       //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint8_t precision>
template<std::floating_point T>
__host__ __device__ inline T phi_calc<precision>::common(T x){
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

std::vector<float> phi_table<0>::values_init(){
	constexpr auto FG_VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	std::vector<float> val(FG_VALUE_RANGE);

	if(!read_values(val)){
		std::cerr<<"phi_table<0>: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
			auto iu = i+LOWER_BOUND_U;
			val[i] = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(std::bit_cast<float>(iu)))));
		}
		if(!write_values(val)) std::cerr<<"phi_table<0>: Caching failed."<<std::endl;
	}
	return val;
}

bool phi_table<0>::read_values(std::vector<float> &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool phi_table<0>::write_values(const std::vector<float> &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

phi_table<0>::phi_table(){
	if(!values) values = values_init();
	values_device = values.devicedata();
}

__host__ __device__ inline float phi_table<0>::get(float x) const{
	auto xa = std::fabs(x);
	//定義域を限定
	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
	#ifdef __CUDA_ARCH__
		auto ya = values_device[reinterpret_cast<uint32_t&>(xa) - LOWER_BOUND_U];
		auto yu = reinterpret_cast<uint32_t&>(ya)|reinterpret_cast<uint32_t&>(x)&0x80000000;
		return reinterpret_cast<float&>(yu);
	#else
		auto ya = values[std::bit_cast<uint32_t>(xa) - LOWER_BOUND_U];
		return std::bit_cast<float>(std::bit_cast<uint32_t>(ya)|std::bit_cast<uint32_t>(x)&0x80000000);
	#endif
}

////////////////////////////////////////////////////////////////
//                                                            //
//                     class phi_table<1>                     //
//                                                            //
////////////////////////////////////////////////////////////////

std::vector<float> phi_table<1>::values_init(){
	constexpr auto FG_VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	std::vector<float> val(FG_VALUE_RANGE);

	if(!read_values(val)){
		std::cerr<<"phi_table<1>: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
			auto x = std::bit_cast<float>(((i<<SHIFT_HALF_FLOAT)-EXPONENT_BIAS)&0x7fffffff);
			auto y = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(x))));
			val[i] = y;
		}
		if(!write_values(val)) std::cerr<<"phi_table<1>: Caching failed."<<std::endl;
	}
	return val;
}

bool phi_table<1>::read_values(std::vector<float> &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool phi_table<1>::write_values(const std::vector<float> &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

phi_table<1>::phi_table(){
	if(!values) values = values_init();
	values_device = values.devicedata();
}

__host__ __device__ inline float phi_table<1>::get(float x) const{
	auto xa = std::fabs(x);
	//定義域を限定
	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
	#ifdef __CUDA_ARCH__
		auto xu = ((reinterpret_cast<uint32_t&>(xa)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT)&UPPER_BOUND_U;
		auto ya = values_device[xu];
		//線形補間
		constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_HALF_FLOAT);
		constexpr float ratiounit = 1.0f/(1<<SHIFT_HALF_FLOAT);
		auto y2 = values_device[xu+1];
		ya += (y2-ya)*static_cast<float>(reinterpret_cast<uint32_t&>(xa)&bottommask)*ratiounit;

		auto yu = reinterpret_cast<uint32_t&>(ya)|reinterpret_cast<uint32_t&>(x)&0x80000000;
		return reinterpret_cast<float&>(yu);
	#else
		auto xu = ((std::bit_cast<std::uint32_t>(xa)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT)&UPPER_BOUND_U;
		auto ya = values[xu];
		//線形補間
		constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_HALF_FLOAT);
		constexpr float ratiounit = 1.0f/(1<<SHIFT_HALF_FLOAT);
		auto y2 = values[xu+1];
		ya += (y2-ya)*static_cast<float>(std::bit_cast<std::uint32_t>(xa)&bottommask)*ratiounit;

		return std::bit_cast<float>(std::bit_cast<uint32_t>(ya)|std::bit_cast<uint32_t>(x)&0x80000000);
	#endif
}

////////////////////////////////////////////////////////////////
//                                                            //
//                     class phi_table<2>                     //
//                                                            //
////////////////////////////////////////////////////////////////

std::vector<float> phi_table<2>::values_init(){
	constexpr auto FG_VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	std::vector<float> val(FG_VALUE_RANGE);

	for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
		auto x = std::bit_cast<float>(((i<<SHIFT_MINI_FLOAT)-EXPONENT_BIAS)&0x7fffffff);
		auto y = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(x))));
		val[i] = y;
	}
	return val;
}

phi_table<2>::phi_table(){
	if(!values) values = values_init();
	values_device = values.devicedata();
}

__host__ __device__ inline float phi_table<2>::get(float x) const{
	auto xa = std::fabs(x);
	//定義域を限定
	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
	#ifdef __CUDA_ARCH__
		auto xu = ((reinterpret_cast<uint32_t&>(xa)+EXPONENT_BIAS)>>SHIFT_MINI_FLOAT)&UPPER_BOUND_U;
		auto ya = values_device[xu];
		//線形補間
		// constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_MINI_FLOAT);
		// constexpr float ratiounit = 1.0f/(1<<SHIFT_MINI_FLOAT);
		// auto y2 = values_device[xu+1];
		// ya += (y2-ya)*static_cast<float>(reinterpret_cast<uint32_t&>(xa)&bottommask)*ratiounit;

		auto yu = reinterpret_cast<uint32_t&>(ya)|reinterpret_cast<uint32_t&>(x)&0x80000000;
		return reinterpret_cast<float&>(yu);
	#else
		auto xu = ((std::bit_cast<std::uint32_t>(xa)+EXPONENT_BIAS)>>SHIFT_MINI_FLOAT)&UPPER_BOUND_U;
		auto ya = values[xu];
		//線形補間
		// constexpr std::uint32_t bottommask = ~static_cast<std::uint32_t>(0)>>(32-SHIFT_MINI_FLOAT);
		// constexpr float ratiounit = 1.0f/(1<<SHIFT_MINI_FLOAT);
		// auto y2 = values[xu+1];
		// ya += (y2-ya)*static_cast<float>(std::bit_cast<std::uint32_t>(xa)&bottommask)*ratiounit;

		return std::bit_cast<float>(std::bit_cast<uint32_t>(ya)|std::bit_cast<uint32_t>(x)&0x80000000);
	#endif
}


}

#endif