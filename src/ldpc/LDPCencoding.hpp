#ifndef INCLUDE_GUARD_ldpc_LDPCencoding
#define INCLUDE_GUARD_ldpc_LDPCencoding

#include <iostream>
#include <bitset>
#include <cstdint>
#include <algorithm>
#include <memory>
#include <cmath>
#include <limits>
#include <exception>
#include "checkmatrix.hpp"

namespace code::LDPC {

template<CheckMatrix T>
class GenerationMatrix_encoding {
	static constexpr std::size_t S = T::sourcesize();
	static constexpr std::size_t C = T::codesize();
	std::vector<std::pair<std::size_t,std::size_t>> Q;//行置換(組織符号にするため)
	std::vector<std::bitset<S>> GT;
	auto GT_product(const std::bitset<S> &vec) const;//matpos1とvecの積を取る
public:
	explicit GenerationMatrix_encoding(const T &H);
	auto systematic_redundancy(const std::bitset<S> &information) const;

	auto substitution(std::bitset<C> vec) const;//Qに従って置換
	template<typename U>
	auto substitution(std::array<U,C> vec) const;//Qに従って置換
	template<std::ranges::random_access_range R>
	auto substitution(R vec) const;//Qに従って置換
	auto inverse_substitution(std::bitset<C> vec) const;//置換を戻す
	template<typename U>
	auto inverse_substitution(std::array<U,C> vec) const;//置換を戻す
	template<std::ranges::random_access_range R>
	auto inverse_substitution(R vec) const;//置換を戻す
};

template<CheckMatrix T>
GenerationMatrix_encoding<T>::GenerationMatrix_encoding(const T &H):
	GT(C-S)
{
	//Hをbitset形式に変換
	std::vector<std::bitset<C>> Hb(C-S);
	for(std::size_t i=0, iend=C-S; i<iend; ++i) for(const auto &j: H[i]) Hb[i].set(j);

	//Hbを掃き出し、右側を単位行列にする
	auto Hp = Hb.rbegin();
	const auto Hend = Hb.rend();
	for(std::size_t pivot=C-1; Hp!=Hend; --pivot,++Hp){
		std::size_t j=pivot;
		auto Hi=Hp;

		do{//pivotとその左を列ごとに捜索
			for(Hi=Hp; Hi!=Hend && !Hi->test(j); ++Hi);//列内でpivotから上の行を下から順に探す
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
			throw std::invalid_argument("LDPC_encoding: invalid check matrix.");
		}
	}

	//単位行列部分は省略して生成行列を作成　Hbの左側の転置=生成行列の右側
	for(std::size_t i=0, iend=C-S; i<iend; ++i){
		auto &Hbi=Hb[i];
		auto &GTi=GT[i];
		for(std::size_t j=0; j<S; ++j) GTi[j] = Hbi[j];
	}
}

template<CheckMatrix T>
auto GenerationMatrix_encoding<T>::GT_product(const std::bitset<S> &vec) const{
	std::bitset<C-S> sol;
	for(std::size_t i=0, iend=C-S; i<iend; ++i){
		std::bitset<S> prod = vec&GT[i];
		auto num = prod.count();
		sol.set(i,static_cast<bool>(num&1));
	}
	return sol;
}

template<CheckMatrix T>
auto GenerationMatrix_encoding<T>::systematic_redundancy(const std::bitset<S> &information) const{
	return GT_product(information);
}

template<CheckMatrix T>
auto GenerationMatrix_encoding<T>::substitution(std::bitset<C> vec) const{
	for(auto [qa, qb]: this->Q){
		bool temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<CheckMatrix T>
template<typename U>
auto GenerationMatrix_encoding<T>::substitution(std::array<U,C> vec) const{
	for(auto [qa, qb]: this->Q){
		U temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<CheckMatrix T>
template<std::ranges::random_access_range R>
auto GenerationMatrix_encoding<T>::substitution(R vec) const{
	if(vec.size()!=T::codesize()) throw std::invalid_argument("invalid input length");
	for(auto [qa, qb]: this->Q){
		typename R::value_type temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<CheckMatrix T>
auto GenerationMatrix_encoding<T>::inverse_substitution(std::bitset<C> vec) const{
	for(auto [qa, qb]: this->Q|std::views::reverse){
		bool temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<CheckMatrix T>
template<typename U>
auto GenerationMatrix_encoding<T>::inverse_substitution(std::array<U,C> vec) const{
	for(auto [qa, qb]: this->Q|std::views::reverse){
		U temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<CheckMatrix T>
template<std::ranges::random_access_range R>
auto GenerationMatrix_encoding<T>::inverse_substitution(R vec) const{
	if(vec.size()!=T::codesize()) throw std::invalid_argument("invalid input length");
	for(auto [qa, qb]: this->Q|std::views::reverse){
		typename R::value_type temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

}

#endif