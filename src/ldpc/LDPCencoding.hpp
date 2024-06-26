﻿#ifndef INCLUDE_GUARD_ldpc_LDPCencoding
#define INCLUDE_GUARD_ldpc_LDPCencoding

#include <iostream>
#include <bitset>
#include <algorithm>
#include <ranges>
#include <cmath>
#include <exception>
#include "LDPCCheckMatrix.hpp"

namespace code::LDPC {

template<CheckMatrix T>
class GenerationMatrix_encoding {
	static constexpr std::uint32_t S = T::sourcesize();
	static constexpr std::uint32_t C = T::codesize();
	static constexpr std::uint32_t Hsize = C-S;

	inline static std::unique_ptr<std::bitset<S>[]> G;//組織符号の生成行列の右側の転置
	inline static std::vector<std::pair<std::uint32_t,std::uint32_t>> Q;//Gに対応する行置換

	static void init(const T &H);
	auto GT_product(const std::bitset<S> &vec) const;//matpos1とvecの積を取る
public:
	explicit GenerationMatrix_encoding(const T &H);
	auto systematic_redundancy(const std::bitset<S> &information) const;
	auto systematic_encode(const std::bitset<S> &information) const;

	auto substitution(std::bitset<C> vec) const;//非組織符号->組織符号
	template<typename U>
	auto substitution(std::array<U,C> vec) const;//非組織符号->組織符号
	auto inverse_substitution(std::bitset<C> vec) const;//組織符号->非組織符号
	template<typename U>
	auto inverse_substitution(std::array<U,C> vec) const;//組織符号->非組織符号
};

////////////////////////////////////////////////////////////////
//                                                            //
//              class GenerationMatrix_encoding               //
//                                                            //
////////////////////////////////////////////////////////////////

template<CheckMatrix T>
void GenerationMatrix_encoding<T>::init(const T &H){
	//Hをbitset形式に変換
	std::vector<std::bitset<C>> Hb(Hsize);
	for(std::uint32_t i=0; i<Hsize; ++i) for(const auto &j: H[i]) Hb[i].set(j);

	{//Hbを掃き出し、右側を単位行列にする
		auto Hp = Hb.rbegin();
		const auto Hend = Hb.rend();
		for(std::uint32_t pivot=C-1; Hp!=Hend; --pivot, ++Hp){

			std::uint32_t j=pivot;
			decltype(Hp) Hi;
			do{//pivotとその左から列ごとに1を探す
				Hi=Hp;
				while(Hi!=Hend && !Hi->test(j)) ++Hi;//列内でpivotから上の行を下から順に探す
			}while(Hi==Hend && --j<C);

			//1が見つからなかった場合ランク落ちで不適格
			if(Hi==Hend) throw std::invalid_argument("LDPC_encoding: invalid check matrix.");

			if(Hp!=Hi){//1が見つかった行とpivot行を交換
				auto temp = *Hp;
				*Hp = *Hi;
				*Hi = temp;
			}
			if(j!=pivot){//1が見つかった列とpivot列を交換
				for(auto &Hk: Hb){
					bool temp = Hk[pivot];
					Hk[pivot] = Hk[j];
					Hk[j] = temp;
				}
				Q.push_back({pivot,j});
			}
			//掃き出し
			for(auto Hk=Hb.rbegin(); Hk!=Hend; ++Hk) if(Hk->test(pivot) && Hk!=Hp) *Hk^=*Hp;
		}
	}

	//単位行列部分は省略して生成行列を作成　Hbの左側の転置=生成行列の右側
	G = std::make_unique<std::bitset<S>[]>(Hsize);
	for(std::uint32_t i=0; i<Hsize; ++i){
		auto &Hbi=Hb[i];
		auto &GTi=G[i];
		for(std::uint32_t j=0; j<S; ++j) GTi[j] = Hbi[j];
	}
}


template<CheckMatrix T>
auto GenerationMatrix_encoding<T>::GT_product(const std::bitset<S> &vec) const{
	std::bitset<Hsize> sol;
	for(std::uint32_t i=0; i<Hsize; ++i){
		std::bitset<S> prod = vec&G[i];
		auto num = prod.count();
		sol.set(i,static_cast<bool>(num&1));
	}
	return sol;
}

template<CheckMatrix T>
GenerationMatrix_encoding<T>::GenerationMatrix_encoding(const T &H){
	if(!G) init(H);
}

template<CheckMatrix T>
auto GenerationMatrix_encoding<T>::systematic_redundancy(const std::bitset<S> &information) const{
	return GT_product(information);
}

template<CheckMatrix T>
auto GenerationMatrix_encoding<T>::systematic_encode(const std::bitset<S> &information) const{
	std::bitset<C> code;
	auto parity = GT_product(information);
	for(std::uint32_t i=0, iend=S; i<iend; ++i) code[i]=information[i];
	for(std::uint32_t i=S, iend=C; i<iend; ++i) code[i]=parity[i-S];
	return code;
}

template<CheckMatrix T>
auto GenerationMatrix_encoding<T>::substitution(std::bitset<C> vec) const{
	for(auto [qa, qb]: Q){
		bool temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<CheckMatrix T>
template<typename U>
auto GenerationMatrix_encoding<T>::substitution(std::array<U,C> vec) const{
	for(auto [qa, qb]: Q){
		U temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<CheckMatrix T>
auto GenerationMatrix_encoding<T>::inverse_substitution(std::bitset<C> vec) const{
	for(auto [qa, qb]: Q|std::views::reverse){
		bool temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

template<CheckMatrix T>
template<typename U>
auto GenerationMatrix_encoding<T>::inverse_substitution(std::array<U,C> vec) const{
	for(auto [qa, qb]: Q|std::views::reverse){
		U temp = vec[qa];
		vec[qa] = vec[qb];
		vec[qb] = temp;
	}
	return vec;
}

}

#endif