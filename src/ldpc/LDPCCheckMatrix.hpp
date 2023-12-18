#ifndef INCLUDE_GUARD_ldpc_LDPCCheckMatrix
#define INCLUDE_GUARD_ldpc_LDPCCheckMatrix

#include <fstream>
#include <charconv>
#include <array>
#include <vector>
#include <string>
#include <cstdint>
#include <memory>
#include <cassert>

namespace code::LDPC {

//以下で定義するクラスが満たすべきconcept
template<class T>
concept CheckMatrix = requires(T x){
	{T::codesize()} -> std::same_as<std::uint32_t>;
	{T::sourcesize()} -> std::same_as<std::uint32_t>;
	{T::size()} -> std::same_as<std::uint32_t>;
	{x.countones()} -> std::same_as<std::uint32_t>;
	{x.colweight(0)} -> std::same_as<std::uint32_t>;
	{x[0][0]} -> std::same_as<const std::uint32_t&>;
	{x.T[0][0]} -> std::same_as<const std::uint32_t&>;
};

template<std::uint32_t S, std::uint32_t C>//S:Source length, C:Code length
class CheckMatrix_irregular{
	static constexpr std::uint32_t Size = C-S;
	static const char *path;

	inline static std::unique_ptr<std::uint32_t[]> col1;//検査行列の1がある列番号
	inline static std::unique_ptr<std::uint32_t[]> idxrow;//検査行列各行の先頭にあたるcol1のインデックス 大きさSize+1
	inline static std::unique_ptr<std::uint32_t[]> row1;//検査行列の1がある行番号
	inline static std::unique_ptr<std::uint32_t[]> idxcol;//検査行列各列の先頭にあたるrow1のインデックス 大きさC+1

	static void readCheckMatrix();//pos1はインスタンスを生成して初めて初期化される
public:
	class sized_ptr {
		std::uint32_t* const ptr;
		std::uint32_t Size;
	public:
		sized_ptr(std::uint32_t* ptr, std::uint32_t size):ptr(ptr),Size(size){}
		inline const std::uint32_t* begin() const noexcept{return ptr;}
		inline const std::uint32_t* end() const noexcept{return ptr+Size;}
		inline const std::uint32_t* data() const noexcept{return ptr;}
		inline std::uint32_t size() const noexcept{return Size;}
		inline const std::uint32_t& operator[](std::uint32_t i) const{
			assert(i<Size);
			return *(ptr+i);
		}
	};

	struct col_ref {
		inline const sized_ptr operator[](std::uint32_t i) const{
			assert(i<C);
			return sized_ptr(row1.get()+idxcol[i],idxcol[i+1]-idxcol[i]);
		}
	} T;

	CheckMatrix_irregular();
	static constexpr bool is_regular() noexcept{return false;}
	static constexpr std::uint32_t codesize() noexcept{return C;}
	static constexpr std::uint32_t sourcesize() noexcept{return S;}
	static constexpr std::uint32_t size() noexcept{return Size;}
	inline std::uint32_t countones() const noexcept{return idxrow[Size];}
	// const std::uint32_t* begin() const noexcept{return col1.get();}
	// const std::uint32_t* end() const noexcept{return col1.get()+idxrow[Size];}
	// const std::uint32_t* data() const noexcept{return col1.get();}
	inline std::uint32_t colweight(std::uint32_t i) const{
		assert(i<C);
		return idxcol[i+1]-idxcol[i];
	}
	const sized_ptr operator[](std::uint32_t i) const{
		assert(i<Size);
		return sized_ptr(col1.get()+idxrow[i],idxrow[i+1]-idxrow[i]);
	}
};

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>//S:Source length, C:Code length, W:row weight
class CheckMatrix_regular{
	static constexpr std::uint32_t Size = C-S;
	static constexpr std::uint32_t Ones = W*Size;
	static constexpr std::uint32_t VW = Ones/C;
	static const char *path;

	inline static std::unique_ptr<std::uint32_t[]> col1;//検査行列の1がある列番号 大きさOnes 幅W
	inline static std::unique_ptr<std::uint32_t[]> row1;//検査行列の1がある列番号 大きさOnes 幅VW

	static void readCheckMatrix();//pos1はインスタンスを生成して初めて初期化される
public:
	template<std::uint32_t Size>
	class sized_ptr {
		const std::uint32_t* ptr;
	public:
		sized_ptr(const std::uint32_t* ptr):ptr(ptr){}
		inline const std::uint32_t* begin() const noexcept{return ptr;}
		inline const std::uint32_t* end() const noexcept{return ptr+Size;}
		inline const std::uint32_t* data() const noexcept{return ptr;}
		static constexpr std::uint32_t size(){return Size;}
		// void operator++(){ptr+=Size;}
		// bool operator==(sized_ptr &rhs) const{return ptr==rhs.ptr;}
		inline const std::uint32_t& operator[](std::uint32_t i) const{
			assert(i<Size);
			return *(ptr+i);
		}
	};

	struct col_ref {
		inline const sized_ptr<VW> operator[](std::uint32_t i) const{
			assert(i<C);
			return sized_ptr<VW>(row1.get()+VW*i);
		}
	} T;

	CheckMatrix_regular();
	static constexpr bool is_regular() noexcept{return true;}
	static constexpr std::uint32_t codesize() noexcept{return C;}
	static constexpr std::uint32_t sourcesize() noexcept{return S;}
	static constexpr std::uint32_t size() noexcept{return Size;}
	static constexpr std::uint32_t countones() noexcept{return W*Size;}
	// const std::uint32_t* begin() const noexcept{return sized_ptr<Size>(col1.get());}
	// const std::uint32_t* end() const noexcept{return sized_ptr<Size>(col1.get()+W*Size);}
	// const std::uint32_t* data() const noexcept{return sized_ptr<Size>(col1.get());}
	static constexpr std::uint32_t colweight(std::uint32_t i){
		assert(i<C);
		return VW;
	}
	const sized_ptr<W> operator[](std::uint32_t i) const{
		assert(i<Size);
		return sized_ptr<W>(col1.get()+W*i);
	}
};

template<std::uint32_t S, std::uint32_t C> struct validCheckMatrixType {using type = void;};
template<std::uint32_t S, std::uint32_t C> using validCheckMatrixType_t = validCheckMatrixType<S,C>::type;//S,Cに応じて適切なCheckMatrixの型を返す

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

template<std::uint32_t S, std::uint32_t C>
void CheckMatrix_irregular<S,C>::readCheckMatrix(){
	idxrow = std::make_unique<std::uint32_t[]>(Size+1);
	idxcol = std::make_unique<std::uint32_t[]>(C+1);

	std::ifstream file(path, std::ios_base::in);
	if(!file.is_open()) throw std::runtime_error("LDPC: cannot open file.");

	std::string buf;
	bool colsizecheck = false;
	std::array<std::uint32_t,C> colidxcount = {};
	std::vector<std::uint32_t> data;
	for(std::uint32_t i=0; i<Size; ++i){
		idxrow[i] = static_cast<std::uint32_t>(data.size());
		if(!std::getline(file, buf)) throw std::runtime_error("Conflict between index data and file content detected.");

		const char *bi = &buf.front();
		const char *const bend = &buf.back();
		while(bi!=bend){
			std::uint32_t val;
			auto r = std::from_chars(bi, bend, val);
			if(r.ec!=std::errc{}) throw std::runtime_error("LDPC: invalid text format.");//読み込みに失敗した場合
			if(val==C) colsizecheck = true;//行列の列数が正しいことを確認
			if(val>C) throw std::runtime_error("Conflict between index data and file content detected.");//列数が多すぎる場合
			data.push_back(val-1);
			++colidxcount[val-1];
			bi=r.ptr;
			while(bi!=bend&&*bi==' ') bi++;//空白を読み飛ばす
		}
	}
	if(std::getline(file, buf)||!colsizecheck) throw std::runtime_error("Conflict between index data and file content detected.");
	file.close();

	idxrow[Size] = static_cast<std::uint32_t>(data.size());
	col1 = std::make_unique<std::uint32_t[]>(idxrow[Size]);
	row1 = std::make_unique<std::uint32_t[]>(idxrow[Size]);
	for(std::uint32_t i=0, iend=idxrow[Size]; i<iend; ++i) col1[i] = data[i];

	idxcol[0] = 0;
	for(std::uint32_t i=0; i<C; ++i) idxcol[i+1] = idxcol[i] + colidxcount[i];

	std::array<std::uint32_t,C> colcount = {};
	for(std::uint32_t i=0; i<Size; ++i) for(std::uint32_t j=idxrow[i], jend=idxrow[i+1]; j<jend; ++j){
		auto k = col1[j];
		row1[idxcol[k]+colcount[k]++] = i;
	}
}

template<std::uint32_t S, std::uint32_t C>
CheckMatrix_irregular<S,C>::CheckMatrix_irregular(){
	if(!idxrow) readCheckMatrix();
}

////////////////////////////////////////////////////////////////
//                                                            //
//                 class CheckMatrix_regular                  //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
void CheckMatrix_regular<S,C,W>::readCheckMatrix(){
	col1 = std::make_unique<std::uint32_t[]>(Ones);
	row1 = std::make_unique<std::uint32_t[]>(Ones);

	std::ifstream file(path, std::ios_base::in);
	if(!file.is_open()) throw std::runtime_error("LDPC: cannot open file.");

	std::string buf;
	bool colsizecheck = false;
	for(std::uint32_t i=0, j=0; i<Size; ++i){
		if(!std::getline(file, buf)) throw std::runtime_error("Conflict between index data and file content detected.");

		const char *bi = &buf.front();
		const char *const bend = &buf.back();
		while(bi!=bend){
			std::uint32_t val;
			if(j==W*(i+1)) throw std::runtime_error("Conflict between index data and file content detected.");//行重みは固定
			auto r = std::from_chars(bi, bend, val);
			if(r.ec!=std::errc{}) throw std::runtime_error("LDPC: invalid text format.");//読み込みに失敗した場合
			if(val==C) colsizecheck = true;//行列の列数が正しいことを確認
			if(val>C) throw std::runtime_error("Conflict between index data and file content detected.");//列数が多すぎる場合
			col1[j++] = val-1;
			bi=r.ptr;
			while(bi!=bend&&*bi==' ') bi++;//空白を読み飛ばす
		}
		if(j!=W*(i+1)) throw std::runtime_error("Conflict between index data and file content detected.");//行重みは固定
	}
	if(std::getline(file, buf)||!colsizecheck) throw std::runtime_error("Conflict between index data and file content detected.");
	file.close();

	std::array<std::uint32_t,C> colcount = {};
	for(std::uint32_t i=0; i<Size; ++i) for(std::uint32_t j=W*i, jend=W*(i+1); j<jend; ++j){
		auto k = col1[j];
		row1[VW*k+colcount[k]++] = i;
	}
}

template<std::uint32_t S, std::uint32_t C, std::uint32_t W>
CheckMatrix_regular<S,C,W>::CheckMatrix_regular(){
	if(!col1) readCheckMatrix();
}

}

#endif