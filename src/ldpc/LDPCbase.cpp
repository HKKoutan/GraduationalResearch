﻿#include "LDPCbase.hpp"

using code::LDPC::CheckMatrix, code::LDPC::func_Gallager;

////////////////////////////////////////////////////////////////
//                                                            //
//                     class CheckMatrix                      //
//                                                            //
////////////////////////////////////////////////////////////////

const char *CheckMatrix<252,504>::path = "H.txt";
const char *CheckMatrix<256,512>::path = "36_512.txt";
const char *CheckMatrix<512,1024>::path = "36_1024.txt";
const char *CheckMatrix<1024,2048>::path = "36_2048.txt";
const char *CheckMatrix<5000,10000>::path = "H2.txt";

////////////////////////////////////////////////////////////////
//                                                            //
//                    class func_Gallager                     //
//                                                            //
////////////////////////////////////////////////////////////////

const std::vector<float> func_Gallager::values = values_init();

std::vector<float> func_Gallager::values_init(){
	constexpr auto FG_VALUE_RANGE = FG_UPPER_BOUND_U - FG_LOWER_BOUND_U + 1u;
	std::vector<float> val(FG_VALUE_RANGE);

	if(!read_values(val)){
		std::cerr<<"func_Gallager: Cache file not found."<<std::endl;
		for(auto i=0ui32; i<FG_VALUE_RANGE; ++i){
			auto iu = i+FG_LOWER_BOUND_U;
			val[i] = static_cast<float>(std::log1p(2.0/std::expm1(static_cast<double>(reinterpret_cast<float&>(iu)))));
		}
		if(!write_values(val)){
			std::cerr<<"func_Gallager: Caching failed."<<std::endl;
		}
	}
	return val;
}

bool func_Gallager::read_values(std::vector<float> &vec){
	std::ifstream file(FG_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool func_Gallager::write_values(const std::vector<float> &vec){
	std::ofstream file(FG_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

// double func_Gallager::operator()(double x) const{
// 	auto xf = std::fabs(static_cast<float>(x));
// 	if(xf<FG_LOWER_BOUND_F) return FG_MAX_VALUE;//定義域を限定
// 	else if(xf>=FG_UPPER_BOUND_F) return FG_MIN_VALUE;
// 	else return values[reinterpret_cast<uint32_t&>(xf) - FG_LOWER_BOUND_U];
// }

float func_Gallager::operator()(float x) const{
	auto y = std::fabs(x);
	//定義域を限定
	if(y<FG_LOWER_BOUND_F) y = FG_LOWER_BOUND_F;
	if(y>FG_UPPER_BOUND_F) y = FG_UPPER_BOUND_F;
	y = values[reinterpret_cast<uint32_t&>(y) - FG_LOWER_BOUND_U];
	// y = static_cast<float>(std::log1p(2.0/std::expm1(y)));
	reinterpret_cast<uint32_t&>(y) |= reinterpret_cast<uint32_t&>(x)&0x80000000;
	return y;
}
