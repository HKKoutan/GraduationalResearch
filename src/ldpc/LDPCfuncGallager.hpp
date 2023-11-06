#ifndef INCLUDE_GUARD_ldpc_LDPCfuncGallager
#define INCLUDE_GUARD_ldpc_LDPCfuncGallager

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <bit>
#include <concepts>
#include <cmath>

namespace code::LDPC {

class func_Gallager_double {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto UPPER_BOUND = 0x1p6f;
public:
	template<std::floating_point T>
	T operator()(T x) const;
};

class func_Gallager_single {
	static constexpr auto LOWER_BOUND = 0x1p-16f;
	static constexpr auto UPPER_BOUND = 0x1p6f;
public:
	template<std::floating_point T>
	T operator()(T x) const;
};

class func_Gallager_table {
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
	func_Gallager_table();
	template<std::floating_point T>
	T operator()(T x) const;
};

class func_Gallager_halftable {
	static constexpr auto LOWER_BOUND = 0x1p-10f;
	static constexpr auto LOWER_BOUND_U = 0x0000;
	static constexpr auto UPPER_BOUND = 0x1p5f;
	static constexpr auto UPPER_BOUND_U = 0xffff;
	static constexpr auto EXPONENT_BIAS = 11<<23;//指数の内部表現の差
	static constexpr auto EXPONENT_ONE = 0x40000000;//指数の内部表現の差
	static constexpr auto SHIFT_HALF_FLOAT = 11;//仮数部の長さの差
	static constexpr auto CACHE_FILENAME = "gallager_half.bin";

	inline static std::vector<float> values;//符号なし 指数4bits(-10~+5)+仮数12bits(ケチ表現)

	static decltype(values) values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	static bool read_values(decltype(values) &vec);
	static bool write_values(const decltype(values) &vec);

public:
	func_Gallager_halftable();
	template<std::floating_point T>
	T operator()(T x) const;
};

////////////////////////////////////////////////////////////////
//                                                            //
//                 class funcGallager_double                  //
//                                                            //
////////////////////////////////////////////////////////////////

template<>
float func_Gallager_double::operator()(float x) const{
	auto y = std::fabs(x);
	//定義域を限定
	if(y<LOWER_BOUND) y = LOWER_BOUND;
	if(y>UPPER_BOUND) y = UPPER_BOUND;
	y = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(y))));
	y = std::bit_cast<float>(std::bit_cast<uint32_t>(y)|std::bit_cast<uint32_t>(x)&0x80000000);
	return y;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                 class funcGallager_single                  //
//                                                            //
////////////////////////////////////////////////////////////////

template<>
float func_Gallager_single::operator()(float x) const{
	auto y = std::fabs(x);
	//定義域を限定
	if(y<LOWER_BOUND) y = LOWER_BOUND;
	if(y>UPPER_BOUND) y = UPPER_BOUND;
	y = std::log1p(2.0f/std::expm1(y));
	y = std::bit_cast<float>(std::bit_cast<uint32_t>(y)|std::bit_cast<uint32_t>(x)&0x80000000);
	return y;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                  class funcGallager_table                  //
//                                                            //
////////////////////////////////////////////////////////////////

decltype(func_Gallager_table::values) func_Gallager_table::values_init(){
	constexpr auto FG_VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	decltype(values) val(FG_VALUE_RANGE);

	if(!read_values(val)){
		std::cerr<<"func_Gallager_table: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
			auto iu = i+LOWER_BOUND_U;
			val[i] = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(std::bit_cast<float>(iu)))));
		}
		if(!write_values(val)){
			std::cerr<<"func_Gallager_table: Caching failed."<<std::endl;
		}
	}
	return val;
}

bool func_Gallager_table::read_values(decltype(values) &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool func_Gallager_table::write_values(const decltype(values) &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

func_Gallager_table::func_Gallager_table(){
	static bool init;
	if(!init){
		values = values_init();
		init = true;
	}
}

template<>
float func_Gallager_table::operator()(float x) const{
	auto y = std::fabs(x);
	//定義域を限定
	if(y<LOWER_BOUND) y = LOWER_BOUND;
	if(y>UPPER_BOUND) y = UPPER_BOUND;
	y = values[std::bit_cast<uint32_t>(y) - LOWER_BOUND_U];
	return std::bit_cast<float>(std::bit_cast<uint32_t>(y)|std::bit_cast<uint32_t>(x)&0x80000000);
}

////////////////////////////////////////////////////////////////
//                                                            //
//                class funcGallager_halftable                //
//                                                            //
////////////////////////////////////////////////////////////////

decltype(func_Gallager_halftable::values) func_Gallager_halftable::values_init(){
	constexpr auto FG_VALUE_RANGE = UPPER_BOUND_U - LOWER_BOUND_U + 1u;
	decltype(values) val(FG_VALUE_RANGE);

	if(!read_values(val)){
		std::cerr<<"func_Gallager_table: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
			auto x = std::bit_cast<float>((((i<<SHIFT_HALF_FLOAT)|EXPONENT_ONE)-EXPONENT_BIAS)&0x7fffffff);
			auto y = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(x))));
			// auto value = static_cast<decltype(values)::value_type>((std::bit_cast<uint32_t>(y)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT);
			val[i] = y;
		}
		if(!write_values(val)){
			std::cerr<<"func_Gallager_table: Caching failed."<<std::endl;
		}
	}
	return val;
}

bool func_Gallager_halftable::read_values(decltype(values) &vec){
	std::ifstream file(CACHE_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool func_Gallager_halftable::write_values(const decltype(values) &vec){
	std::ofstream file(CACHE_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

func_Gallager_halftable::func_Gallager_halftable(){
	static bool init;
	if(!init){
		values = values_init();
		init = true;
	}
}

template<>
float func_Gallager_halftable::operator()(float x) const{
	auto y = std::fabs(x);
	//定義域を限定
	if(y<LOWER_BOUND) y = LOWER_BOUND;
	if(y>UPPER_BOUND) y = UPPER_BOUND;
	auto xu = ((std::bit_cast<std::uint32_t>(y)+EXPONENT_BIAS)>>SHIFT_HALF_FLOAT)&0x0000ffff;
	auto yu = values[xu];
	return std::bit_cast<float>(std::bit_cast<uint32_t>(yu)|std::bit_cast<uint32_t>(x)&0x80000000);
}

}

#endif