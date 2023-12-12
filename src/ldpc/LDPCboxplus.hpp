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
		sum_accumlator():abssum(0),signprod(0){}
		void operator+=(T rhs);
		T operator-(T rhs) const;
	};

	template<std::floating_point T>
	class min_accumlator{
		std::pair<T,T> absmin;
		uint_of_length_t<T> signprod;
	public:
		min_accumlator():absmin(std::numeric_limits<T>::infinity(),std::numeric_limits<T>::infinity()),signprod(0){}
		void operator+=(T rhs);
		T operator-(T rhs) const;
	};
}

struct minsum {
	template<std::floating_point T>
	using accumlator = min_accumlator<T>;
	template<std::floating_point T>
	static inline T forward(T x){return x;}
	template<std::floating_point T>
	static inline T backward(T x){return x;}
};

template<std::uint8_t precision=0>
class funcGallager_calc {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto UPPER_BOUND = 0x1p6f;

	template<std::floating_point T>
	static T common(T x);
public:
	template<std::floating_point T>
	using accumlator = sum_accumlator<T>;
	template<std::floating_point T>
	static inline T forward(T x){return common(x);}
	template<std::floating_point T>
	static inline T backward(T x){return common(x);}
};

template<std::uint8_t precision=0> class funcGallager_table {
	static_assert(precision<3,"invalid precision specification");
};

template<>
class funcGallager_table<0> {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto LOWER_BOUND_U = std::bit_cast<uint32_t>(LOWER_BOUND);
	static constexpr auto UPPER_BOUND = 0x1p6f;
	static constexpr auto UPPER_BOUND_U = std::bit_cast<uint32_t>(UPPER_BOUND);
	static constexpr auto CACHE_FILENAME = "gallager_float.bin";

	inline static std::vector<float> values;

	static decltype(values) values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	static bool read_values(decltype(values) &vec);
	static bool write_values(const decltype(values) &vec);
	float get(float x) const;
public:
	template<std::floating_point T>
	using accumlator = sum_accumlator<T>;
	funcGallager_table();
	template<std::floating_point T>
	inline T forward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	template<std::floating_point T>
	inline T backward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
};

template<>
class funcGallager_table<1> {
	static constexpr auto LOWER_BOUND = 0x1p-10f;
	static constexpr auto LOWER_BOUND_U = 0x00000000;
	static constexpr auto UPPER_BOUND = 0x1.ffdp4f;
	static constexpr auto UPPER_BOUND_U = 0x0000ffff;
	static constexpr auto EXPONENT_BIAS = (11+0x00000080)<<23;//指数の内部表現の差
	static constexpr auto SHIFT_HALF_FLOAT = 11;//仮数部の長さの差
	static constexpr auto CACHE_FILENAME = "gallager_half.bin";

	inline static std::vector<float> values;//インデックスは符号なし浮動小数点数[指数4bits(-10~+5)+仮数12bits(ケチ表現)]の内部表現

	static decltype(values) values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	static bool read_values(decltype(values) &vec);
	static bool write_values(const decltype(values) &vec);
	float get(float x) const;
public:
	template<std::floating_point T>
	using accumlator = sum_accumlator<T>;
	funcGallager_table();
	template<std::floating_point T>
	inline T forward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	template<std::floating_point T>
	inline T backward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
};

template<>
class funcGallager_table<2> {
	static constexpr auto LOWER_BOUND = 0x1p-10f;
	static constexpr auto LOWER_BOUND_U = 0x00000000;
	static constexpr auto UPPER_BOUND = 0x1.ffdp4f;
	static constexpr auto UPPER_BOUND_U = 0x000000ff;
	static constexpr auto EXPONENT_BIAS = (11+0x00000080)<<23;//指数の内部表現の差
	static constexpr auto SHIFT_MINI_FLOAT = 19;//仮数部の長さの差

	inline static std::vector<float> values;//インデックスは符号なし浮動小数点数[指数4bits(-10~+5)+仮数12bits(ケチ表現)]の内部表現

	static decltype(values) values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	float get(float x) const;
public:
	template<std::floating_point T>
	using accumlator = sum_accumlator<T>;
	funcGallager_table();
	template<std::floating_point T>
	inline T forward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
	template<std::floating_point T>
	inline T backward(T x) const{return static_cast<T>(get(static_cast<float>(x)));}
};

////////////////////////////////////////////////////////////////
//                                                            //
//                    class sum_accumlator                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::floating_point T>
inline void sum_accumlator<T>::operator+=(T rhs){
	abssum += std::fabs(rhs);
	signprod ^= std::bit_cast<uint_of_length_t<T>>(rhs);
}

template<std::floating_point T>
inline T sum_accumlator<T>::operator-(T rhs) const{
	constexpr uint_of_length_t<T> signmask = 1u<<(sizeof(T)*8u-1u);
	return std::bit_cast<T>((signprod^std::bit_cast<uint_of_length_t<T>>(rhs))&signmask|std::bit_cast<uint_of_length_t<T>>(abssum-std::fabs(rhs)));
}

////////////////////////////////////////////////////////////////
//                                                            //
//                    class min_accumlator                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::floating_point T>
inline void min_accumlator<T>::operator+=(T rhs){
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
inline T min_accumlator<T>::operator-(T rhs) const{
	constexpr uint_of_length_t<T> signmask = 1u<<(sizeof(T)*8u-1u);
	auto arhs = std::fabs(rhs);
	auto val = absmin.first==arhs?absmin.second:absmin.first;
	return std::bit_cast<T>((signprod^std::bit_cast<uint_of_length_t<T>>(rhs))&signmask|std::bit_cast<uint_of_length_t<T>>(val));
}

////////////////////////////////////////////////////////////////
//                                                            //
//                  class funcGallager_calc                   //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint8_t precision>
template<std::floating_point T>
inline T funcGallager_calc<precision>::common(T x){
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
//                class funcGallager_table<0>                 //
//                                                            //
////////////////////////////////////////////////////////////////

decltype(funcGallager_table<0>::values) funcGallager_table<0>::values_init(){
	constexpr auto FG_VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	decltype(values) val(FG_VALUE_RANGE);

	if(!read_values(val)){
		std::cerr<<"funcGallager_table: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
			auto iu = i+LOWER_BOUND_U;
			val[i] = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(std::bit_cast<float>(iu)))));
		}
		if(!write_values(val)){
			std::cerr<<"funcGallager_table: Caching failed."<<std::endl;
		}
	}
	return val;
}

bool funcGallager_table<0>::read_values(decltype(values) &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool funcGallager_table<0>::write_values(const decltype(values) &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

funcGallager_table<0>::funcGallager_table(){
	static bool init;
	if(!init){
		values = values_init();
		init = true;
	}
}

inline float funcGallager_table<0>::get(float x) const{
	auto xa = std::fabs(x);
	//定義域を限定
	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
	auto ya = values[std::bit_cast<uint32_t>(xa) - LOWER_BOUND_U];
	return std::bit_cast<float>(std::bit_cast<uint32_t>(ya)|std::bit_cast<uint32_t>(x)&0x80000000);
}

////////////////////////////////////////////////////////////////
//                                                            //
//                class funcGallager_table<1>                 //
//                                                            //
////////////////////////////////////////////////////////////////

decltype(funcGallager_table<1>::values) funcGallager_table<1>::values_init(){
	constexpr auto FG_VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	decltype(values) val(FG_VALUE_RANGE);

	if(!read_values(val)){
		std::cerr<<"funcGallager_halftable: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
			auto x = std::bit_cast<float>(((i<<SHIFT_HALF_FLOAT)-EXPONENT_BIAS)&0x7fffffff);
			auto y = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(x))));
			// auto value = static_cast<decltype(values)::value_type>((std::bit_cast<uint32_t>(y)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT);
			val[i] = y;
		}
		if(!write_values(val)){
			std::cerr<<"funcGallager_halftable: Caching failed."<<std::endl;
		}
	}
	return val;
}

bool funcGallager_table<1>::read_values(decltype(values) &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool funcGallager_table<1>::write_values(const decltype(values) &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

funcGallager_table<1>::funcGallager_table(){
	static bool init;
	if(!init){
		values = values_init();
		init = true;
	}
}

inline float funcGallager_table<1>::get(float x) const{
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

////////////////////////////////////////////////////////////////
//                                                            //
//                class funcGallager_table<2>                 //
//                                                            //
////////////////////////////////////////////////////////////////

decltype(funcGallager_table<2>::values) funcGallager_table<2>::values_init(){
	constexpr auto FG_VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	decltype(values) val(FG_VALUE_RANGE);

	for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
		auto x = std::bit_cast<float>(((i<<SHIFT_MINI_FLOAT)-EXPONENT_BIAS)&0x7fffffff);
		auto y = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(x))));
		val[i] = y;
	}
	return val;
}

funcGallager_table<2>::funcGallager_table(){
	static bool init;
	if(!init){
		values = values_init();
		init = true;
	}
}

inline float funcGallager_table<2>::get(float x) const{
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


}

#endif