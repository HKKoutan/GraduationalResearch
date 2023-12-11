#ifndef INCLUDE_GUARD_ldpc_LDPCboxplus
#define INCLUDE_GUARD_ldpc_LDPCboxplus

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <bit>
#include <concepts>
#include <cmath>

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
}

class funcGallager_calc {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto UPPER_BOUND = 0x1p6f;

	template<std::uint8_t precision, std::floating_point T>
	static T calc(T x);
public:
	template<std::uint8_t precision=0, std::floating_point T>
	static T forward(T x){return calc<precision>(x);}
	template<std::uint8_t precision=0, std::floating_point T>
	static T backward(T x){return calc<precision>(x);}

	template<std::floating_point T>
	class accumlator{
		T abssum;
		uint_of_length_t<T> signprod;
	public:
		accumlator():abssum(0),signprod(0){}
		void operator+=(T rhs);
		T operator-(T rhs) const;
	};
};

class funcGallager_table {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto LOWER_BOUND_U = std::bit_cast<uint32_t>(LOWER_BOUND);
	static constexpr auto UPPER_BOUND = 0x1p6f;
	static constexpr auto UPPER_BOUND_U = std::bit_cast<uint32_t>(UPPER_BOUND);
	static constexpr auto CACHE_FILENAME = "gallager_float.bin";

	inline static std::vector<float> values;

	static decltype(values) values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	static bool read_values(decltype(values) &vec);
	static bool write_values(const decltype(values) &vec);
public:
	funcGallager_table();
	template<std::floating_point T>
	T operator()(T x) const;
};

class funcGallager_halftable {
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
	template<std::floating_point T>
	static T interpolate(T x);
public:
	funcGallager_halftable();
	template<std::floating_point T>
	T operator()(T x) const;
};

////////////////////////////////////////////////////////////////
//                                                            //
//                  class funcGallager_calc                   //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint8_t precision, std::floating_point T>
T funcGallager_calc::calc(T x){
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

template<std::floating_point T>
void funcGallager_calc::accumlator<T>::operator+=(T rhs){
	abssum += std::fabs(rhs);
	signprod ^= std::bit_cast<uint_of_length_t<T>>(rhs);
}

template<std::floating_point T>
T funcGallager_calc::accumlator<T>::operator-(T rhs) const{
	constexpr uint_of_length_t<T> signmask = 1u<<(sizeof(T)*8u-1u);
	return std::bit_cast<T>((signprod^std::bit_cast<uint_of_length_t<T>>(rhs))&signmask|std::bit_cast<uint_of_length_t<T>>(abssum-std::fabs(rhs)));
}

////////////////////////////////////////////////////////////////
//                                                            //
//                  class funcGallager_table                  //
//                                                            //
////////////////////////////////////////////////////////////////

decltype(funcGallager_table::values) funcGallager_table::values_init(){
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

bool funcGallager_table::read_values(decltype(values) &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool funcGallager_table::write_values(const decltype(values) &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

funcGallager_table::funcGallager_table(){
	static bool init;
	if(!init){
		values = values_init();
		init = true;
	}
}

template<>
float funcGallager_table::operator()(float x) const{
	auto xa = std::fabs(x);
	//定義域を限定
	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
	auto ya = values[std::bit_cast<uint32_t>(xa) - LOWER_BOUND_U];
	return std::bit_cast<float>(std::bit_cast<uint32_t>(ya)|std::bit_cast<uint32_t>(x)&0x80000000);
}

////////////////////////////////////////////////////////////////
//                                                            //
//                class funcGallager_halftable                //
//                                                            //
////////////////////////////////////////////////////////////////

decltype(funcGallager_halftable::values) funcGallager_halftable::values_init(){
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

bool funcGallager_halftable::read_values(decltype(values) &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool funcGallager_halftable::write_values(const decltype(values) &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

template<>
float funcGallager_halftable::interpolate(float x){
	auto xu = ((std::bit_cast<std::uint32_t>(x)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT)&0x0000ffff;
	auto ratio = static_cast<float>(std::bit_cast<std::uint32_t>(x)&0x000007ff)*0x1p-11f;
	auto y1 = values[xu];
	auto y2 = values[xu+1];
	return y1 + (y1-y2)*ratio;
}

funcGallager_halftable::funcGallager_halftable(){
	static bool init;
	if(!init){
		values = values_init();
		init = true;
	}
}

template<>
float funcGallager_halftable::operator()(float x) const{
	auto xa = std::fabs(x);
	//定義域を限定
	if(xa<LOWER_BOUND) xa = LOWER_BOUND;
	if(xa>UPPER_BOUND) xa = UPPER_BOUND;
	// auto xu = ((std::bit_cast<std::uint32_t>(xa)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT)&0x0000ffff;
	// auto ya = values[xu];
	auto ya = interpolate(xa);
	return std::bit_cast<float>(std::bit_cast<uint32_t>(ya)|std::bit_cast<uint32_t>(x)&0x80000000);
}

}

#endif