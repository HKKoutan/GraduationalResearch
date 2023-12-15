#ifndef INCLUDE_GUARD_ldpc_LDPCCheckMatrix
#define INCLUDE_GUARD_ldpc_LDPCCheckMatrix

#include <fstream>
#include <charconv>
#include <array>
#include <vector>
#include <string>
#include <ranges>
#include <memory>
#include <cassert>

namespace code::LDPC {

//以下で定義するクラスが満たすべきconcept
template<class T>
concept CheckMatrix = requires(T x){
	{T::codesize()} -> std::same_as<std::size_t>;
	{T::sourcesize()} -> std::same_as<std::size_t>;
	{x.size()} -> std::same_as<std::size_t>;
};

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class CheckMatrix_irregular{
	static constexpr std::size_t Size = C-S;
	static const char *path;

	inline static std::unique_ptr<std::size_t[]> col1;//検査行列の1がある列番号
	inline static std::unique_ptr<std::size_t[]> rowidx;//検査行列各行の先頭にあたるcol1のインデックス 大きさSize+1

	static void readCheckMatrix();//pos1はインスタンスを生成して初めて初期化される
public:
	class sized_ptr {
		std::size_t* const ptr;
		std::size_t Size;
	public:
		sized_ptr(std::size_t* ptr, std::size_t size):ptr(ptr),Size(size){}
		const std::size_t* begin() const noexcept{return ptr;}
		const std::size_t* end() const noexcept{return ptr+Size;}
		const std::size_t* data() const noexcept{return ptr;}
		std::size_t size() const noexcept{return Size;}
		inline const std::size_t& operator[](std::size_t i) const{
			assert(i<Size);
			return *(ptr+i);
		}
	};

	CheckMatrix_irregular();
	static constexpr bool is_regular() noexcept{return false;}
	static constexpr std::size_t codesize() noexcept{return C;}
	static constexpr std::size_t sourcesize() noexcept{return S;}
	static constexpr std::size_t size() noexcept{return Size;}
	std::size_t countones() const noexcept{return rowidx[Size];}
	// const std::size_t* begin() const noexcept{return col1.get();}
	// const std::size_t* end() const noexcept{return col1.get()+rowidx[Size];}
	// const std::size_t* data() const noexcept{return col1.get();}
	const sized_ptr operator[](std::size_t i) const{
		assert(i<Size);
		return sized_ptr(col1.get()+rowidx[i],rowidx[i+1]-rowidx[i]);
	}
};

template<std::size_t S, std::size_t C, std::size_t W>//S:Source length, C:Code length, W:row weight
class CheckMatrix_regular{
	static constexpr std::size_t Size = C-S;
	static constexpr std::size_t Ones = W*Size;
	static constexpr std::size_t VW = Ones/C;
	static const char *path;

	inline static std::unique_ptr<std::size_t[]> col1;//検査行列の1がある列番号 大きさOnes 幅W

	static void readCheckMatrix();//pos1はインスタンスを生成して初めて初期化される
public:
	template<std::size_t Size>
	class sized_ptr {
		const std::size_t* ptr;
	public:
		sized_ptr(const std::size_t* ptr):ptr(ptr){}
		const std::size_t* begin() const noexcept{return ptr;}
		const std::size_t* end() const noexcept{return ptr+Size;}
		const std::size_t* data() const noexcept{return ptr;}
		static constexpr std::size_t size(){return Size;}
		// void operator++(){ptr+=Size;}
		// bool operator==(sized_ptr &rhs) const{return ptr==rhs.ptr;}
		inline const std::size_t& operator[](std::size_t i) const{
			assert(i<Size);
			return *(ptr+i);
		}
	};

	CheckMatrix_regular();
	static constexpr bool is_regular() noexcept{return true;}
	static constexpr std::size_t codesize() noexcept{return C;}
	static constexpr std::size_t sourcesize() noexcept{return S;}
	static constexpr std::size_t size() noexcept{return Size;}
	static constexpr std::size_t countones() noexcept{return W*Size;}
	// const std::size_t* begin() const noexcept{return sized_ptr<Size>(col1.get());}
	// const std::size_t* end() const noexcept{return sized_ptr<Size>(col1.get()+W*Size);}
	// const std::size_t* data() const noexcept{return sized_ptr<Size>(col1.get());}
	const sized_ptr<W> operator[](std::size_t i) const{
		assert(i<Size);
		return sized_ptr<W>(col1.get()+W*i);
	}
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
SET_CheckMatrix_irregular("H3.txt",5001,10000)

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

	rowidx = std::make_unique<std::size_t[]>(Size+1);
	std::vector<std::size_t> data;

	std::string buf;
	bool colsizecheck = false;
	for(std::size_t i=0; i<Size; ++i){
		rowidx[i] = data.size();
		if(!std::getline(file, buf)) throw std::runtime_error("Conflict between index data and file content detected.");

		const char *bi = &buf.front();
		const char *const bend = &buf.back();
		while(bi!=bend){
			std::size_t val;
			auto r = std::from_chars(bi, bend, val);
			if(r.ec!=std::errc{}) throw std::runtime_error("LDPC: invalid text format.");//読み込みに失敗した場合
			if(val==C) colsizecheck = true;//行列の列数が正しいことを確認
			if(val>C) throw std::runtime_error("Conflict between index data and file content detected.");//列数が多すぎる場合
			data.push_back(val-1);
			bi=r.ptr;
			while(bi!=bend&&*bi==' ') bi++;//空白を読み飛ばす
		}
		
	}
	if(std::getline(file, buf)||!colsizecheck) throw std::runtime_error("Conflict between index data and file content detected.");
	file.close();

	rowidx[Size] = data.size();
	col1 = std::make_unique<std::size_t[]>(rowidx[Size]);
	for(std::size_t i=0, iend=rowidx[Size]; i<iend; ++i) col1[i] = data[i];
}

template<std::size_t S, std::size_t C>
CheckMatrix_irregular<S,C>::CheckMatrix_irregular(){
	if(!col1) readCheckMatrix();
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

	col1 = std::make_unique<std::size_t[]>(Ones);

	std::string buf;
	bool colsizecheck = false;
	std::size_t i=0;
	for(std::size_t j=0; j<Size; ++j){
		if(!std::getline(file, buf)) throw std::runtime_error("Conflict between index data and file content detected.");

		const char *bi = &buf.front();
		const char *const bend = &buf.back();
		while(bi!=bend){
			uint64_t val;
			if(i==W*(j+1)) throw std::runtime_error("Conflict between index data and file content detected.");//行重みは固定
			auto r = std::from_chars(bi, bend, val);
			if(r.ec!=std::errc{}) throw std::runtime_error("LDPC: invalid text format.");//読み込みに失敗した場合
			if(val==C) colsizecheck = true;//行列の列数が正しいことを確認
			if(val>C) throw std::runtime_error("Conflict between index data and file content detected.");//列数が多すぎる場合
			col1[i++] = val-1;
			bi=r.ptr;
			while(bi!=bend&&*bi==' ') bi++;//空白を読み飛ばす
		}
		if(i!=W*(j+1)) throw std::runtime_error("Conflict between index data and file content detected.");//行重みは固定
	}
	if(std::getline(file, buf)||!colsizecheck) throw std::runtime_error("Conflict between index data and file content detected.");
	file.close();
}

template<std::size_t S, std::size_t C, std::size_t W>
CheckMatrix_regular<S,C,W>::CheckMatrix_regular(){
	if(!col1) readCheckMatrix();
}

}

#endif