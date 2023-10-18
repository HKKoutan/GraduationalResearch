#ifndef __channel_SEQUENCER__
#define __channel_SEQUENCER__

#include <array>
#include <cmath>
#include <random>
#include <exception>
#include "DNASnttype.hpp"

namespace channel{

template<std::uint8_t ATGC=0x1B>
class Nanopore_Sequencing{
	std::mt19937_64 mt;
	std::uniform_real_distribution<> uniform;
	static constexpr std::uint8_t A = (ATGC>>6)&3, T = (ATGC>>4)&3, G = (ATGC>>2)&3, C = ATGC&3;
	static constexpr double A0 = (A>>1?-1.0:1.0), A1 = (A&1?-1.0:1.0), T0 = (T>>1?-1.0:1.0), T1 = (T&1?-1.0:1.0), G0 = (G>>1?-1.0:1.0), G1 = (G&1?-1.0:1.0), C0 = (C>>1?-1.0:1.0), C1 = (C&1?-1.0:1.0); 
	const std::array<double,4> error_rate;//error_rate[0]~error_rate[3] = p1~p4
	const std::array<double,4> non_error_rate;//ATGCそれぞれ誤らない確率
	const std::array<std::array<double,4>,4> condprob;//condprob[a][b] P(a|b) 条件付き確率
	auto condprob_init() const;

public:
	Nanopore_Sequencing(double alpha, std::int64_t seed);
	Nanopore_Sequencing(double alpha);

	template<std::size_t L>
	auto noise(const std::array<code::DNAS::nucleotide_t,L> &in);

	template<typename FPTYPE = typename float, std::size_t S, std::size_t R>
	auto VLRLL_LLR(const std::array<code::DNAS::nucleotide_t,S> &cm, const std::array<code::DNAS::nucleotide_t,R> &cr, code::DNAS::nucleotide_t initial_state) const;
};

auto Nanopore_Sequencing<0x1B>::condprob_init() const{
	return std::array<std::array<double,4>,4>{
		non_error_rate[0],error_rate[1],error_rate[3],error_rate[2],
		error_rate[1],non_error_rate[1],error_rate[2],error_rate[0],
		error_rate[3],error_rate[2],non_error_rate[2],error_rate[1],
		error_rate[2],error_rate[0],error_rate[1],non_error_rate[3]
	};
}

template<std::uint8_t ATGC>
Nanopore_Sequencing<ATGC>::Nanopore_Sequencing(double alpha) : Nanopore_Sequencing(alpha, 0){}

template<std::uint8_t ATGC>
Nanopore_Sequencing<ATGC>::Nanopore_Sequencing(double alpha, std::int64_t seed):
	mt(seed),
	uniform(0,1),
	error_rate{4*alpha, alpha, 0.01, 0},
	non_error_rate{1-error_rate[1]-error_rate[2]-error_rate[3], 1-error_rate[0]-error_rate[1]-error_rate[2], 1-error_rate[1]-error_rate[2]-error_rate[3], 1-error_rate[0]-error_rate[1]-error_rate[2]},//ATGCそれぞれの誤らない確率
	condprob(condprob_init())
{
	static_assert((A!=T)&&(A!=G)&&(A!=C)&&(T!=G)&&(T!=C)&&(G!=C));
	if(alpha>0.198||alpha<0) throw std::out_of_range("alpha out of range");
}

template<std::uint8_t ATGC>
template<std::size_t L>
auto Nanopore_Sequencing<ATGC>::noise(const std::array<code::DNAS::nucleotide_t,L> &in){
	auto out = in;
	for(auto &i: out){
		code::DNAS::nucleotide_t j=0;
		double rand = uniform(mt)-condprob[i][j];
		while(rand>=condprob[i][j]){
			j+=1;
			rand-=condprob[i][j];
		}
		i=j;
	}
	return out;
}

template<std::uint8_t ATGC>
template<typename FPTYPE, std::size_t S, std::size_t R>
auto Nanopore_Sequencing<ATGC>::VLRLL_LLR(const std::array<code::DNAS::nucleotide_t,S> &cm, const std::array<code::DNAS::nucleotide_t,R> &cr, code::DNAS::nucleotide_t initial_state) const{
	std::array<FPTYPE,S*2+R*2> LLR;
	std::size_t j=0u;
	//情報部
	for(const auto &current: cm){
		switch(current){
		case A:
			LLR[j++] = static_cast<FPTYPE>(A0*std::log((non_error_rate[0]+error_rate[1])/(error_rate[3]+error_rate[2])));
			LLR[j++] = static_cast<FPTYPE>(A1*std::log((non_error_rate[0]+error_rate[3])/(error_rate[1]+error_rate[2])));
			break;

		case T:
			LLR[j++] = static_cast<FPTYPE>(T0*std::log((non_error_rate[1]+error_rate[1])/(error_rate[0]+error_rate[2])));
			LLR[j++] = static_cast<FPTYPE>(T1*std::log((non_error_rate[1]+error_rate[0])/(error_rate[1]+error_rate[2])));
			break;

		case G:
			LLR[j++] = static_cast<FPTYPE>(G0*std::log((non_error_rate[2]+error_rate[1])/(error_rate[3]+error_rate[2])));
			LLR[j++] = static_cast<FPTYPE>(G1*std::log((non_error_rate[2]+error_rate[3])/(error_rate[1]+error_rate[2])));
			break;

		case C:
			LLR[j++] = static_cast<FPTYPE>(C0*std::log((non_error_rate[3]+error_rate[1])/(error_rate[0]+error_rate[2])));
			LLR[j++] = static_cast<FPTYPE>(C1*std::log((non_error_rate[3]+error_rate[0])/(error_rate[1]+error_rate[2])));
			break;
		}
	}
	//冗長部
	constexpr double p3_to_11 = 1.0/22.0;
	constexpr double p3_to_10 = 21.0/22.0;

	for(auto previous = initial_state; const auto &current: cr){
		double P0X =//遷移語が1 or 2になる組み合わせ
			condprob[previous][previous] * (condprob[previous+1][current] + condprob[previous+2][current]) +
			condprob[previous-1][previous] * (condprob[previous][current] + condprob[previous+1][current]) +
			condprob[previous-2][previous] * (condprob[previous-1][current] + condprob[previous][current]) +
			condprob[previous-3][previous] * (condprob[previous-2][current] + condprob[previous-1][current]);
		double P1X =//遷移語が0 or 3になる組み合わせ
			condprob[previous][previous] * (condprob[previous][current] + condprob[previous+3][current]) +
			condprob[previous-1][previous] * (condprob[previous-1][current] + condprob[previous+2][current]) +
			condprob[previous-2][previous] * (condprob[previous-2][current] + condprob[previous+1][current]) +
			condprob[previous-3][previous] * (condprob[previous-3][current] + condprob[previous][current]);
		double PX0 =//遷移語が1 or 3(*p3_to_10)になる組み合わせ
			(condprob[current-1][previous] + p3_to_10*condprob[current-3][previous]) * condprob[current][current] +
			(condprob[current][previous] + p3_to_10*condprob[current-2][previous]) * condprob[current+1][current] +
			(condprob[current+1][previous] + p3_to_10*condprob[current-1][previous]) * condprob[current+2][current] +
			(condprob[current+2][previous] + p3_to_10*condprob[current][previous]) * condprob[current+3][current];
		double PX1 =//遷移語が0 or 2 or 3(*p3_to_11)になる組み合わせ
			(condprob[current][previous] + condprob[current-2][previous] + p3_to_11*condprob[current-3][previous]) * condprob[current][current] +
			(condprob[current+1][previous] + condprob[current-1][previous] + p3_to_11*condprob[current-2][previous]) * condprob[current+1][current] +
			(condprob[current+2][previous] + condprob[current][previous] + p3_to_11*condprob[current-1][previous]) * condprob[current+2][current] +
			(condprob[current+3][previous] + condprob[current+1][previous] + p3_to_11*condprob[current][previous]) * condprob[current+3][current];

		previous = current;
		LLR[j++] = static_cast<FPTYPE>(std::log(P0X/P1X));
		LLR[j++] = static_cast<FPTYPE>(std::log(PX0/PX1));
	}
	return LLR;
}

}

#endif