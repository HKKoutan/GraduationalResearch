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
	{x.codesize()} -> std::same_as<std::size_t>;
	{x.sourcesize()} -> std::same_as<std::size_t>;
	{x.size()} -> std::same_as<std::size_t>;
	{x.begin()} -> std::random_access_iterator;
	{x.cbegin()} -> std::random_access_iterator;
	{x.end()} -> std::random_access_iterator;
	{x.cend()} -> std::random_access_iterator;
	{x[0]} -> std::ranges::random_access_range;
};

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class CheckMatrix_irregular{
	static std::array<std::vector<std::uint64_t>,C-S> pos1;//行ごとに検査行列の1がある列番号を列挙
	static const char *path;
	static void readCheckMatrix();
public:
	// using reference = typename decltype(pos1)::reference;
	// using const_reference = typename decltype(pos1)::const_reference;
	// using iterator = typename decltype(pos1)::iterator;
	// using const_iterator = typename decltype(pos1)::const_iterator;
	// using value_type = typename decltype(pos1)::value_type;
	CheckMatrix_irregular();
	static constexpr auto codesize() noexcept{return C;}
	static constexpr auto sourcesize() noexcept{return S;}
	static constexpr auto size() noexcept{return C-S;}
	static constexpr auto begin() noexcept{return pos1.cbegin();}
	static constexpr auto cbegin() noexcept{return pos1.cbegin();}
	static constexpr auto end() noexcept{return pos1.cend();}
	static constexpr auto cend() noexcept{return pos1.cend();}
	constexpr const auto &operator[](std::size_t x) const noexcept{return pos1[x];}
};

template<std::size_t S, std::size_t C, std::size_t W>//S:Source length, C:Code length, W:row weight
class CheckMatrix_regular{
	static std::array<std::array<std::uint64_t,W>,C-S> pos1;//行ごとに検査行列の1がある列番号を列挙
	static const char *path;
	static void readCheckMatrix();
public:
	// using reference = typename decltype(pos1)::reference;
	// using const_reference = typename decltype(pos1)::const_reference;
	// using iterator = typename decltype(pos1)::iterator;
	// using const_iterator = typename decltype(pos1)::const_iterator;
	// using value_type = typename decltype(pos1)::value_type;
	CheckMatrix_regular();
	static constexpr auto codesize() noexcept{return C;}
	static constexpr auto sourcesize() noexcept{return S;}
	static constexpr auto size() noexcept{return C-S;}
	static constexpr auto begin() noexcept{return pos1.cbegin();}
	static constexpr auto cbegin() noexcept{return pos1.cbegin();}
	static constexpr auto end() noexcept{return pos1.cend();}
	static constexpr auto cend() noexcept{return pos1.cend();}
	constexpr const auto &operator[](std::size_t x) const noexcept{return pos1[x];}
};

template<std::size_t S, std::size_t C> auto getCheckMatrix();//S,Cに応じて適切なCheckMatrixのインスタンスを返す

template<std::size_t S, std::size_t C> struct CheckMatrixType {using type = std::invoke_result_t<getCheckMatrix>;};
template<std::size_t S, std::size_t C> using CheckMatrixType_t = CheckMatrixType<S,C>::type;//S,Cに応じて適切なCheckMatrixの型を返す

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