#ifndef INCLUDE_GUARD_ldpc_checkmatrix
#define INCLUDE_GUARD_ldpc_checkmatrix

#include <fstream>
#include <charconv>
#include <array>
#include <vector>
#include <string>
#include <iterator>
#include <ranges>
#include <type_traits>

namespace code::LDPC {

//以下で定義するクラスが満たすべきconcept
template<class T>
concept CheckMatrix = requires(std::remove_reference_t<T>& x){
	{T::codesize()} -> std::same_as<std::size_t>;
	{T::sourcesize()} -> std::same_as<std::size_t>;
	{x.size()} -> std::same_as<std::size_t>;
	{x.begin()} -> std::random_access_iterator;
	{x.cbegin()} -> std::random_access_iterator;
	{x.end()} -> std::random_access_iterator;
	{x.cend()} -> std::random_access_iterator;
	{x[0]} -> std::ranges::random_access_range;
};

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class CheckMatrix_irregular{
	inline static std::array<std::vector<std::uint64_t>,C-S> pos1;//行ごとに検査行列の1がある列番号を列挙
	static const char *path;
	static void readCheckMatrix();//pos1はインスタンスを生成して初めて初期化される
public:
	CheckMatrix_irregular();
	static constexpr auto codesize() noexcept{return C;}
	static constexpr auto sourcesize() noexcept{return S;}
	static constexpr auto size() noexcept{return C-S;}
	constexpr auto begin() noexcept{return pos1.cbegin();}
	constexpr auto cbegin() noexcept{return pos1.cbegin();}
	constexpr auto end() noexcept{return pos1.cend();}
	constexpr auto cend() noexcept{return pos1.cend();}
	constexpr const auto &operator[](std::size_t x) const noexcept{return pos1[x];}
};

template<std::size_t S, std::size_t C, std::size_t W>//S:Source length, C:Code length, W:row weight
class CheckMatrix_regular{
	inline static std::array<std::array<std::uint64_t,W>,C-S> pos1;//行ごとに検査行列の1がある列番号を列挙
	static const char *path;
	static void readCheckMatrix();//pos1はインスタンスを生成して初めて初期化される
public:
	CheckMatrix_regular();
	static constexpr auto codesize() noexcept{return C;}
	static constexpr auto sourcesize() noexcept{return S;}
	static constexpr auto size() noexcept{return C-S;}
	constexpr auto begin() noexcept{return pos1.cbegin();}
	constexpr auto cbegin() noexcept{return pos1.cbegin();}
	constexpr auto end() noexcept{return pos1.cend();}
	constexpr auto cend() noexcept{return pos1.cend();}
	constexpr const auto &operator[](std::size_t x) const noexcept{return pos1[x];}
};

template<std::size_t S, std::size_t C> struct validCheckMatrixType {using type = void;};
template<std::size_t S, std::size_t C> using validCheckMatrixType_t = validCheckMatrixType<S,C>::type;//S,Cに応じて適切なCheckMatrixの型を返す

//特殊化の定義
#define SET_CheckMatrix_regular(P,S,C,W) const char *CheckMatrix_regular<S,C,W>::path = P; template<> struct validCheckMatrixType<S,C> {using type = CheckMatrix_regular<S,C,W>;};
#define SET_CheckMatrix_irregular(P,S,C) const char *CheckMatrix_irregular<S,C>::path = P; template<> struct validCheckMatrixType<S,C> {using type = CheckMatrix_irregular<S,C>;};

SET_CheckMatrix_regular("H.txt",252,504,6)
SET_CheckMatrix_regular("36_512.txt",256,512,6)
SET_CheckMatrix_regular("36_1024.txt",512,1024,6)
SET_CheckMatrix_regular("36_2048.txt",1024,2048,6)
SET_CheckMatrix_regular("H2.txt",5000,10000,6)
SET_CheckMatrix_irregular("H3.txt",4999,10000)

#undef SET_CheckMatrix_regular
#undef SET_CheckMatrix_irregular

////////////////////////////////////////////////////////////////
//                                                            //
//                class CheckMatrix_irregular                 //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C>
void CheckMatrix_irregular<S,C>::readCheckMatrix(){
	std::ifstream file(path, std::ios_base::in);
	if(!file.is_open()) throw std::runtime_error("LDPC: cannot open file.");

	std::string buf;
	std::size_t colsize = 0;
	for(auto &row: pos1){
		if(!std::getline(file, buf)) throw std::runtime_error("Conflict between index data and file content detected.");

		const char *bi = &buf.front();
		const char *const bend = &buf.back();
		while(bi!=bend){
			uint64_t val;
			auto r = std::from_chars(bi, bend, val);
			//読み込みに失敗した場合
			if(r.ec!=std::errc{}) throw std::runtime_error("LDPC: invalid text format.");
			//行列の列数を記録
			if(val>colsize) colsize=val;
			row.push_back(val-1);
			bi=r.ptr;
			while(bi!=bend&&*bi==' ') bi++;//空白を読み飛ばす
		}
	}
	if(std::getline(file, buf)||colsize!=C) throw std::runtime_error("Conflict between index data and file content detected.");
	file.close();
}

template<std::size_t S, std::size_t C>
CheckMatrix_irregular<S,C>::CheckMatrix_irregular(){
	static bool read;
	if(!read){
		readCheckMatrix();
		read = true;
	}
}

////////////////////////////////////////////////////////////////
//                                                            //
//                 class CheckMatrix_regular                  //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C, std::size_t W>
void CheckMatrix_regular<S,C,W>::readCheckMatrix(){
	std::ifstream file(path, std::ios_base::in);
	if(!file.is_open()) throw std::runtime_error("LDPC: cannot open file.");

	std::string buf;
	std::size_t colsize = 0;
	for(auto &row: pos1){
		if(!std::getline(file, buf)) throw std::runtime_error("Conflict between index data and file content detected.");

		const char *bi = &buf.front();
		const char *const bend = &buf.back();
		std::size_t i=0;
		while(bi!=bend){
			uint64_t val;
			if(i==W) throw std::runtime_error("Conflict between index data and file content detected [].");//行重みは固定
			auto r = std::from_chars(bi, bend, val);
			//読み込みに失敗した場合
			if(r.ec!=std::errc{}) throw std::runtime_error("LDPC: invalid text format.");
			//行列の列数を記録
			if(val>colsize) colsize=val;
			row[i++] = val-1;
			bi=r.ptr;
			while(bi!=bend&&*bi==' ') bi++;//空白を読み飛ばす
		}
	}
	if(std::getline(file, buf)||colsize!=C) throw std::runtime_error("Conflict between index data and file content detected.");
	file.close();
}

template<std::size_t S, std::size_t C, std::size_t W>
CheckMatrix_regular<S,C,W>::CheckMatrix_regular(){
	static bool read;
	if(!read){
		readCheckMatrix();
		read = true;
	}
}

}

#endif