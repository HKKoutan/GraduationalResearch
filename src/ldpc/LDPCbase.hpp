#ifndef INCLUDE_GUARD_ldpc_LDPCbase
#define INCLUDE_GUARD_ldpc_LDPCbase

#include <iostream>
#include <bitset>
#include <cstdint>
#include <bit>
#include <algorithm>
#include <memory>
#include <cmath>
#include <limits>
#include <exception>
#include "LDPCCheckMatrix.hpp"

namespace code::LDPC {

using fptype = float;

//decoder class

template<CheckMatrix T>
class I_LDPC_Decoding {
protected:
	static constexpr std::size_t S = T::sourcesize();
	static constexpr std::size_t C = T::codesize();
	const T H;//検査行列
	std::vector<std::pair<std::array<fptype,C>,std::array<fptype,C>>> alphabeta;
	std::array<std::pair<std::vector<fptype*>,std::vector<const fptype*>>,C-S> alphapbetap;

	virtual void rowupdate() = 0;
private:
	//メンバ初期化関数
	auto alphabeta_init() const;
	auto alphapbetap_init();
public:
	explicit I_LDPC_Decoding(const T &H) noexcept;
	virtual ~I_LDPC_Decoding(){}
	void decode_init();//decodeで使用する変数の初期化
	bool iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR);
	auto estimate(const std::array<fptype,C> &LEVR) const;//推定符号語を求める
};

class func_Gallager_std {
	static constexpr auto FG_LOWER_BOUND_F = 0x1p-16f;
	static constexpr auto FG_UPPER_BOUND_F = 0x1p6f;
public:
	// double operator()(double x) const;
	float operator()(float x) const;
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

	static const std::vector<float> values;
public:
	// double operator()(double x) const;
	float operator()(float x) const;
};

inline const std::vector<float> func_Gallager_table::values = func_Gallager_table::values_init();

template<CheckMatrix T>
class SumProduct_Decoding : public I_LDPC_Decoding<T> {
	const func_Gallager_table fg;

	using signtype = std::uint32_t;
	static_assert(sizeof(fptype)==sizeof(signtype));
	static constexpr signtype signmask = 1u<<(sizeof(signtype)*8u-1u);
public:
	explicit SumProduct_Decoding(const T &H): I_LDPC_Decoding<T>(H), fg(){}
	void rowupdate() override;//alpha,betaを更新する(行処理)
};

template<CheckMatrix T>
class MinSum_Decoding : public I_LDPC_Decoding<T> {
	using signtype = std::uint32_t;
	static_assert(sizeof(fptype)==sizeof(signtype));
	static constexpr signtype signmask = 1u<<(sizeof(signtype)*8u-1u);
public:
	explicit MinSum_Decoding(const T &H): I_LDPC_Decoding<T>(H){}
	void rowupdate() override;//alpha,betaを更新する(行処理)
};

////////////////////////////////////////////////////////////////
//                                                            //
//                   class I_LDPC_Decoding                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
auto I_LDPC_Decoding<T>::alphabeta_init() const{
	//HT(temporary variable)
	std::array<std::size_t,C> Hheight{};
	for(const auto &Hi: H) for(auto j: Hi) ++Hheight[j];
	return std::vector(std::ranges::max(Hheight), decltype(alphabeta)::value_type());
}

template<CheckMatrix T>
auto I_LDPC_Decoding<T>::alphapbetap_init(){
	constexpr auto isize = C-S;
	decltype(alphapbetap) apbp{};
	//HT(temporary variable)
	std::array<std::vector<std::uint64_t>,C> HT{};
	for(std::size_t i=0u, iend=isize; i<iend; ++i) for(auto j: H[i]) HT[j].push_back(i);

	for(std::size_t i=0u, iend=isize; i<iend; ++i){
		auto &Hi = H[i];
		auto &[api, bpi] = apbp[i];
		//Hとalphabetapの要素数を揃える
		api.resize(Hi.size());
		bpi.resize(Hi.size());
		//alpha<-alphap beta<-betap
		for(std::size_t j=0u, jend=Hi.size(); j<jend; ++j){
			auto &hij = Hi[j];
			auto &Hj = HT[hij];
			auto &[ai, bi] = alphabeta[std::ranges::find(Hj, i)-Hj.begin()];
			api[j] = &ai[hij];
			bpi[j] = &bi[hij];
		}
	}
	return apbp;
}

template<CheckMatrix T>
I_LDPC_Decoding<T>::I_LDPC_Decoding(const T &H) noexcept:
	H(H),
	alphabeta(alphabeta_init()),
	alphapbetap(alphapbetap_init())
{}

template<CheckMatrix T>
void I_LDPC_Decoding<T>::decode_init(){
	for(auto &[ai, bi]: alphabeta) for(auto &bij: bi) bij = 0;
}

template<CheckMatrix T>
bool I_LDPC_Decoding<T>::iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR){
	//apply LLR
	for(auto &[ai, bi]: alphabeta) for(std::size_t j=0u, jend=C; j<jend; ++j) bi[j] += LLR[j];
	//row update
	rowupdate();
	//column update
	for(auto &lpj: LPR) lpj = 0;
	for(auto &[ai, bi]: alphabeta) for(std::size_t j=0u, jend=C; j<jend; ++j) LPR[j] += ai[j];
	for(auto &[ai, bi]: alphabeta) for(std::size_t j=0u, jend=C; j<jend; ++j) bi[j] = LPR[j]-ai[j];
	for(std::size_t j=0u, jend=C; j<jend; ++j) LPR[j] += LLR[j];
	//parity check
	for(const auto &Hi : H){
		auto parity = false;
		for(const auto &j : Hi) parity ^= LPR[j]<0;
		if(parity) return false;
	}
	return true;
}

template<CheckMatrix T>
auto I_LDPC_Decoding<T>::estimate(const std::array<fptype,C> &LEVR) const{
	std::bitset<S> est;
	for(std::size_t j=0u, jend=S; j<jend; ++j) est[j] = LEVR[j]<0;
	return est;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                 class SumProduct_Decoding                  //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
void SumProduct_Decoding<T>::rowupdate(){
	for(auto &[ai, bi]: this->alphabeta) for(auto &bij: bi) bij = fg(bij);
	for(auto &[api, bpi]: this->alphapbetap){
		fptype abssum = 0;
		signtype signprod = 0;
		for(auto bpij: bpi){
			auto bij = *bpij;
			abssum += std::fabs(bij);
			signprod ^= std::bit_cast<signtype>(bij);
		}
		for(std::size_t j=0u, jend=api.size(); j<jend; ++j){
			auto bij = *bpi[j];
			auto absval = fg(abssum-std::fabs(bij));
			auto sign = (std::bit_cast<signtype>(bij)^signprod)&signmask;
			*api[j] = std::bit_cast<fptype>(sign|std::bit_cast<signtype>(absval));
		}
	}
}

////////////////////////////////////////////////////////////////
//                                                            //
//                   class MinSum_Decoding                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
void MinSum_Decoding<T>::rowupdate(){
	for(auto &[api, bpi]: this->alphapbetap){
		signtype signprod = 0;
		for(auto bpij: bpi){
			signprod ^= std::bit_cast<signtype>(*bpij);
		}
		for(std::size_t j=0u, jend=api.size(); j<jend; ++j){
			auto min = std::numeric_limits<fptype>::infinity();
			for(std::size_t k=0u, kend=api.size(); k<kend; ++k) if(j != k){
				auto temp = std::fabs(*bpi[k]);
				if(temp<min) min = temp;
			}
			auto sign = (std::bit_cast<const signtype>(*bpi[j])^signprod)&signmask;
			*api[j] = std::bit_cast<fptype>(sign|std::bit_cast<signtype>(min));
		}
	}
}

}

#endif