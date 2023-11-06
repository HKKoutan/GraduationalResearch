#ifndef INCLUDE_GUARD_ldpc_LDPCfuncGallager
#define INCLUDE_GUARD_ldpc_LDPCfuncGallager

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <bit>
#include <cmath>

namespace code::LDPC {

class func_Gallager_std {
	static constexpr auto FG_LOWER_BOUND_F = 0x1p-10f;
	static constexpr auto FG_UPPER_BOUND_F = 0x1p5f;
public:
	template<std::floating_point T>
	T operator()(T x) const;
};

class func_Gallager_table {
	static constexpr auto FG_LOWER_BOUND_F = 0x1p-16f;
	static constexpr auto FG_LOWER_BOUND_U = std::bit_cast<uint32_t>(0x1p-16f);// FG_LOWER_BOUND_Fの内部表現
	static constexpr auto FG_UPPER_BOUND_F = 0x1p6f;
	static constexpr auto FG_UPPER_BOUND_U = std::bit_cast<uint32_t>(0x1p6f);// FG_UPPER_BOUND_Fの内部表現
	static constexpr auto FG_FILENAME = "gallager_float.bin";

	static std::vector<float> values_init();//キャッシュファイルを読み込み値を返す。失敗したら、値を計算してキャッシュファイルに保存する。
	static bool read_values(std::vector<float> &vec);
	static bool write_values(const std::vector<float> &vec);

	inline static std::vector<float> values;
public:
	func_Gallager_table();
	template<std::floating_point T>
	T operator()(T x) const;
};

////////////////////////////////////////////////////////////////
//                                                            //
//                   class funcGallager_std                   //
//                                                            //
////////////////////////////////////////////////////////////////

template<>
float func_Gallager_std::operator()(float x) const{
	auto y = std::fabs(x);
	//定義域を限定
	if(y<FG_LOWER_BOUND_F) y = FG_LOWER_BOUND_F;
	if(y>FG_UPPER_BOUND_F) y = FG_UPPER_BOUND_F;
	y = static_cast<float>(std::log1p(2.0f/std::expm1(y)));
	y = std::bit_cast<float>(std::bit_cast<uint32_t>(y)|std::bit_cast<uint32_t>(x)&0x80000000);
	return y;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                  class funcGallager_table                  //
//                                                            //
////////////////////////////////////////////////////////////////

std::vector<float> func_Gallager_table::values_init(){
	constexpr auto FG_VALUE_RANGE = FG_UPPER_BOUND_U - FG_LOWER_BOUND_U + 1u;
	std::vector<float> val(FG_VALUE_RANGE);

	if(!read_values(val)){
		std::cerr<<"func_Gallager_table: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
			auto iu = i+FG_LOWER_BOUND_U;
			val[i] = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(std::bit_cast<float>(iu)))));
		}
		if(!write_values(val)){
			std::cerr<<"func_Gallager_table: Caching failed."<<std::endl;
		}
	}
	return val;
}

bool func_Gallager_table::read_values(std::vector<float> &vec){
	std::ifstream file(FG_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool func_Gallager_table::write_values(const std::vector<float> &vec){
	std::ofstream file(FG_FILENAME, std::ios::out | std::ios::binary);

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
	if(y<FG_LOWER_BOUND_F) y = FG_LOWER_BOUND_F;
	if(y>FG_UPPER_BOUND_F) y = FG_UPPER_BOUND_F;
	y = values[std::bit_cast<uint32_t>(y) - FG_LOWER_BOUND_U];
	y = std::bit_cast<float>(std::bit_cast<uint32_t>(y)|std::bit_cast<uint32_t>(x)&0x80000000);
	return y;
}

}

#endif