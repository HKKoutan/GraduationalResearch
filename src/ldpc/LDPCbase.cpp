#include "LDPCbase.hpp"

using code::LDPC::CheckMatrix, code::LDPC::func_Gallager_std, code::LDPC::func_Gallager_table;

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
//                  class func_Gallager_std                   //
//                                                            //
////////////////////////////////////////////////////////////////

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
//                 class func_Gallager_table                  //
//                                                            //
////////////////////////////////////////////////////////////////

std::vector<float> func_Gallager_table::values::values_init(){
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

bool func_Gallager_table::values::read_values(std::vector<float> &vec){
	std::ifstream file(FG_FILENAME, std::ios::in | std::ios::binary);
	if(!file.is_open()) return false;

	file.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

bool func_Gallager_table::values::write_values(const std::vector<float> &vec){
	std::ofstream file(FG_FILENAME, std::ios::out | std::ios::binary);

	file.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(vec.front()));
	file.close();
	if(file.fail()) return false;
	return true;
}

// double func_Gallager_table::operator()(double x) const{
// 	auto xf = std::fabs(static_cast<float>(x));
// 	if(xf<FG_LOWER_BOUND_F) return FG_MAX_VALUE;//定義域を限定
// 	else if(xf>=FG_UPPER_BOUND_F) return FG_MIN_VALUE;
// 	else return values[reinterpret_cast<uint32_t&>(xf) - FG_LOWER_BOUND_U];
// }

float func_Gallager_table::operator()(float x) const{
	auto y = std::fabs(x);
	//定義域を限定
	if(y<FG_LOWER_BOUND_F) y = FG_LOWER_BOUND_F;
	if(y>FG_UPPER_BOUND_F) y = FG_UPPER_BOUND_F;
	y = table[std::bit_cast<uint32_t>(y) - FG_LOWER_BOUND_U];
	y = std::bit_cast<float>(std::bit_cast<uint32_t>(y)|std::bit_cast<uint32_t>(x)&0x80000000);
	return y;
}
