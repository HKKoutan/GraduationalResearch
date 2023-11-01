#ifndef INCLUDE_GUARD_ldpc_LDPCbase
#define INCLUDE_GUARD_ldpc_LDPCbase

#include <iostream>
#include <fstream>
#include <charconv>
#include <array>
#include <bitset>
#include <vector>
#include <string>
#include <cstdint>
#include <bit>
#include <algorithm>
#include <ranges>
#include <memory>
#include <cmath>
#include <limits>
#include <exception>
#include "checkmatrix.hpp"

namespace code::LDPC {

using fptype = float;

//encorder class

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class I_LDPC_Encoding {
protected:
	std::vector<std::pair<std::size_t,std::size_t>> Q;//行置換(組織符号にするため)
public:
	I_LDPC_Encoding(){}
	virtual ~I_LDPC_Encoding(){}
	virtual std::bitset<C-S> systematic_redundancy(const std::bitset<S> &information) const = 0;//return parity bits
	auto substitution(std::array<fptype,C> vec) const;//Qに従って置換
	auto substitution(std::bitset<C> vec) const;//Qに従って置換
	auto inverse_substitution(std::array<fptype,C> vec) const;//置換を戻す
	auto inverse_substitution(std::bitset<C> vec) const;//置換を戻す
};

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class Generation_Matrix_Encoding : public I_LDPC_Encoding<S,C> {
	std::vector<std::bitset<S>> GT;

	auto GT_product(const std::bitset<S> &vec) const;//matpos1とvecの積を取る
public:
	explicit Generation_Matrix_Encoding(const CheckMatrixType_t<S,C> &H);
	std::bitset<C-S> systematic_redundancy(const std::bitset<S> &information) const override;
};

//decoder class

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class I_LDPC_Decoding {
protected:
	const CheckMatrixType_t<S,C> H;//検査行列
	std::vector<std::pair<std::array<fptype,C>,std::array<fptype,C>>> alphabeta;
	std::array<std::pair<std::vector<fptype*>,std::vector<const fptype*>>,C-S> alphapbetap;

	virtual void rowupdate() = 0;
private:
	//メンバ初期化関数
	auto alphabeta_init() const;
	auto alphapbetap_init();
public:
	explicit I_LDPC_Decoding(const CheckMatrixType_t<S,C> &H) noexcept;
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

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class SumProduct_Decoding : public I_LDPC_Decoding<S,C> {
	const func_Gallager_table fg;

	using signtype = std::uint32_t;
	static_assert(sizeof(fptype)==sizeof(signtype));
	static constexpr signtype signmask = 1u<<(sizeof(signtype)*8u-1u);
public:
	explicit SumProduct_Decoding(decltype(I_LDPC_Decoding<S,C>::H) H): I_LDPC_Decoding<S,C>(H), fg(){}
	void rowupdate() override;//alpha,betaを更新する(行処理)
};

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class MinSum_Decoding : public I_LDPC_Decoding<S,C> {
	using signtype = std::uint32_t;
	static_assert(sizeof(fptype)==sizeof(signtype));
	static constexpr signtype signmask = 1u<<(sizeof(signtype)*8u-1u);
public:
	explicit MinSum_Decoding(decltype(I_LDPC_Decoding<S,C>::H) H): I_LDPC_Decoding<S,C>(H){}
	void rowupdate() override;//alpha,betaを更新する(行処理)
};

////////////////////////////////////////////////////////////////
//                                                            //
//                   class I_LDPC_Encoding                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C>
auto I_LDPC_Encoding<S,C>::substitution(std::array<fptype,C> vec) const{
	for(auto [qa, qb]: this->Q){
		auto temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<std::size_t S, std::size_t C>
auto I_LDPC_Encoding<S,C>::substitution(std::bitset<C> vec) const{
	for(auto [qa, qb]: this->Q){
		bool temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<std::size_t S, std::size_t C>
auto I_LDPC_Encoding<S,C>::inverse_substitution(std::array<fptype,C> vec) const{
	for(auto [qa, qb]: this->Q|std::views::reverse){
		auto temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<std::size_t S, std::size_t C>
auto I_LDPC_Encoding<S,C>::inverse_substitution(std::bitset<C> vec) const{
	for(auto [qa, qb]: this->Q|std::views::reverse){
		bool temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

////////////////////////////////////////////////////////////////
//                                                            //
//              class Generation_Matrix_Encoding              //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C>
Generation_Matrix_Encoding<S,C>::Generation_Matrix_Encoding(const CheckMatrixType_t<S,C> &H):
	I_LDPC_Encoding<S,C>(),
	GT(C-S)
{
	constexpr auto isize = H.size();
	//Hをbitset形式に変換
	std::vector<std::bitset<C>> Hb(isize);
	for(std::size_t i=0u, iend=isize; i<iend; ++i) for(const auto &j: H[i]) Hb[i].set(j);

	//Hbを掃き出し、右側を単位行列にする
	auto Hp = Hb.rbegin();
	const auto Hend = Hb.rend();
	for(std::size_t pivot=C-1; Hp!=Hend; --pivot,++Hp){
		std::size_t j=pivot;
		auto Hi=Hp;

		do{//pivotとその左を捜索
			for(Hi=Hp; Hi!=Hend && !Hi->test(j); ++Hi);//同じ行に無かったら上の行から探す
		}while(--j<C && Hi==Hend);

		if(Hi!=Hend){//最後まで走査する前に1を見つけられた場合
			if(Hp!=Hi){//1が見つかった行とpivot行を交換
				auto temp = *Hp;
				*Hp = *Hi;
				*Hi = temp;
			}
			if(++j!=pivot){//1が見つかった列とpivot列が異なる場合　++はdo-while内の--jとの辻褄合わせ
				for(auto &Hk: Hb){//1が見つかった列とpivot列を交換
					bool temp = Hk[pivot];
					Hk[pivot] = Hk[j];
					Hk[j] = temp;
				}
				this->Q.push_back({pivot,j});
			}
			//掃き出し
			for(auto Hk=Hb.rbegin(); Hk!=Hend; ++Hk) if(Hk->test(pivot) && Hk!=Hp) *Hk^=*Hp;
		}else{//1が見つからなかった場合ランク落ちで不適格
			throw std::runtime_error("LDPC: invalid check matrix.");
		}
	}

	//単位行列部分は省略して生成行列を作成　Hbの左側の転置=生成行列の右側
	for(std::size_t i=0u, iend=isize; i<iend; ++i){
		auto &Hbi=Hb[i];
		auto &GTi=GT[i];
		for(std::size_t j=0u, jsize=S; j<jsize; ++j) GTi[j] = Hbi[j];
	}
}

template<std::size_t S, std::size_t C>
auto Generation_Matrix_Encoding<S,C>::GT_product(const std::bitset<S> &vec) const{
	std::bitset<C-S> sol;
	for(std::size_t i=0u, iend=C-S; i<iend; ++i){
		std::bitset<S> prod = vec&GT[i];
		auto num = prod.count();
		sol.set(i,static_cast<bool>(num&1));
	}
	return sol;
}

template<std::size_t S, std::size_t C>
std::bitset<C-S> Generation_Matrix_Encoding<S,C>::systematic_redundancy(const std::bitset<S> &information) const{
	return GT_product(information);
}

////////////////////////////////////////////////////////////////
//                                                            //
//                   class I_LDPC_Decoding                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C>
auto I_LDPC_Decoding<S,C>::alphabeta_init() const{
	//HT(temporary variable)
	std::array<std::size_t,C> Hheight{};
	for(const auto &Hi: H) for(auto j: Hi) ++Hheight[j];
	return std::vector(std::ranges::max(Hheight), decltype(alphabeta)::value_type());
}

template<std::size_t S, std::size_t C>
auto I_LDPC_Decoding<S,C>::alphapbetap_init(){
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

template<std::size_t S, std::size_t C>
I_LDPC_Decoding<S,C>::I_LDPC_Decoding(const CheckMatrixType_t<S,C> &H) noexcept:
	H(H),
	alphabeta(alphabeta_init()),
	alphapbetap(alphapbetap_init())
{}

template<std::size_t S, std::size_t C>
void I_LDPC_Decoding<S,C>::decode_init(){
	for(auto &[ai, bi]: alphabeta) for(auto &bij: bi) bij = 0;
}

template<std::size_t S, std::size_t C>
bool I_LDPC_Decoding<S,C>::iterate(std::array<fptype,C> &LPR, const std::array<fptype,C> &LLR){
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

template<std::size_t S, std::size_t C>
auto I_LDPC_Decoding<S,C>::estimate(const std::array<fptype,C> &LEVR) const{
	std::bitset<S> est;
	for(std::size_t j=0u, jend=S; j<jend; ++j) est[j] = LEVR[j]<0;
	return est;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                 class SumProduct_Decoding                  //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C>
void SumProduct_Decoding<S,C>::rowupdate(){
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

template<std::size_t S, std::size_t C>
void MinSum_Decoding<S,C>::rowupdate(){
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