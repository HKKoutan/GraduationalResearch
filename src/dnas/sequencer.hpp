#ifndef __channel_SEQUENCER__
#define __channel_SEQUENCER__

// #include <iostream>
// #include <vector>
#include <array>
#include <cstdint>
#include <cmath>
#include <random>
#include <exception>
#include "DNASnttype.hpp"

namespace channel{

class Nanopore_Sequencing{
	std::mt19937_64 mt;
	std::uniform_real_distribution<> uniform;
	const std::array<double,4> error_rate;//error_rate[0]~error_rate[3] = p1~p4
	const std::array<double,4> non_error_rate;//ATGCそれぞれ誤らない確率
	const std::array<std::array<double,4>,4> condprob;//condprob[a][b] P(a|b) 条件付き確率

public:
	Nanopore_Sequencing(double alpha, std::int64_t seed);
	Nanopore_Sequencing(double alpha);

	template<std::size_t L>
	void noise(std::array<code::DNAS::nucleotide_t,L> &c);

	template<std::size_t S>
	void message_LLR(const std::array<code::DNAS::nucleotide_t,S/2> &cm, std::array<double,S> &LLRm) const;
	template<std::size_t R>
	void redundancy_LLR(const std::array<code::DNAS::nucleotide_t,R/2> &cr, std::array<double,R> &LLRr, code::DNAS::nucleotide_t initial_state) const;
};

Nanopore_Sequencing::Nanopore_Sequencing(double alpha) : Nanopore_Sequencing(alpha, 0){}

Nanopore_Sequencing::Nanopore_Sequencing(double alpha, std::int64_t seed)
	: mt(seed), uniform(0,1), error_rate{4*alpha, alpha, 0.01, 0}, non_error_rate{1-error_rate[1]-error_rate[2]-error_rate[3], 1-error_rate[0]-error_rate[1]-error_rate[2], 1-error_rate[1]-error_rate[2]-error_rate[3], 1-error_rate[0]-error_rate[1]-error_rate[2]},
	condprob{
	non_error_rate[0],error_rate[1],error_rate[3],error_rate[2],
	error_rate[1],non_error_rate[1],error_rate[2],error_rate[0],
	error_rate[3],error_rate[2],non_error_rate[2],error_rate[1],
	error_rate[2],error_rate[0],error_rate[1],non_error_rate[3]
	}{
	if(alpha>0.198||alpha<0){
		throw std::out_of_range("alpha out of range");
	}
}

template<std::size_t L>
void Nanopore_Sequencing::noise(std::array<code::DNAS::nucleotide_t,L> &c){
	for(auto &i: c){
		//std::cout<<int(i)<<endl;
		nucleotide_t j=0;
		double rand = uniform(mt)-condprob[i][j];
		while(rand>=condprob[i][j]){
			j+=1;
			rand-=condprob[i][j];
		}
		i=j;
	}
}

template<std::size_t S>
void Nanopore_Sequencing::message_LLR(const std::array<code::DNAS::nucleotide_t,S/2> &cm, std::array<double,S> &LLRm) const{
	for(std::size_t i=0u, j=0u; i<S/2; ++i){
		switch(cm[i]){
		case 0:
			LLRm[j++] = log((non_error_rate[0]+error_rate[1])/(error_rate[3]+error_rate[2]));
			LLRm[j++] = log((non_error_rate[0]+error_rate[3])/(error_rate[1]+error_rate[2]));
			break;

		case 1:
			LLRm[j++] = log((non_error_rate[1]+error_rate[1])/(error_rate[0]+error_rate[2]));
			LLRm[j++] = -log((non_error_rate[1]+error_rate[0])/(error_rate[1]+error_rate[2]));
			break;

		case 2:
			LLRm[j++] = -log((non_error_rate[2]+error_rate[1])/(error_rate[3]+error_rate[2]));
			LLRm[j++] = log((non_error_rate[2]+error_rate[3])/(error_rate[1]+error_rate[2]));
			break;

		case 3:
			LLRm[j++] = -log((non_error_rate[3]+error_rate[1])/(error_rate[0]+error_rate[2]));
			LLRm[j++] = -log((non_error_rate[3]+error_rate[0])/(error_rate[1]+error_rate[2]));
			break;
		}
	}
}

template<std::size_t R>
void Nanopore_Sequencing::redundancy_LLR(const std::array<code::DNAS::nucleotide_t,R/2> &cr, std::array<double,R> &LLRr, code::DNAS::nucleotide_t initial_state) const{
	constexpr double p3_to_11 = 1.0/22.0;
	constexpr double p3_to_10 = 21.0/22.0;

	nucleotide_t previous = initial_state;

	for(std::size_t i=0u, iend=S/2 j=0u; i<iend; ++i){
/*	condprob{
	non_error_rate[0],error_rate[1],error_rate[3],error_rate[2],
	error_rate[1],non_error_rate[1],error_rate[2],error_rate[0],
	error_rate[3],error_rate[2],non_error_rate[2],error_rate[1],
	error_rate[2],error_rate[0],error_rate[1],non_error_rate[3]}
*/
		double P0X =//遷移語が1 or 2になる組み合わせ
			condprob[previous][previous] * (condprob[(previous+1)&3][current] + condprob[(previous+2)&3][current]) +
			condprob[(previous-1)&3][previous] * (condprob[previous][current] + condprob[(previous+1)&3][current]) +
			condprob[(previous-2)&3][previous] * (condprob[(previous-1)&3][current] + condprob[previous][current]) +
			condprob[(previous-3)&3][previous] * (condprob[(previous-2)&3][current] + condprob[(previous-1)&3][current]);
		double P1X =//遷移語が0 or 3になる組み合わせ
			condprob[previous][previous] * (condprob[previous][current] + condprob[(previous+3)&3][current]) +
			condprob[(previous-1)&3][previous] * (condprob[(previous-1)&3][current] + condprob[(previous+2)&3][current]) +
			condprob[(previous-2)&3][previous] * (condprob[(previous-2)&3][current] + condprob[(previous+1)&3][current]) +
			condprob[(previous-3)&3][previous] * (condprob[(previous-3)&3][current] + condprob[previous][current]);
		double PX0 =//遷移語が1 or 3(*p3_to_10)になる組み合わせ
			(condprob[(current-1)&3][previous] + p3_to_10*condprob[(current-3)&3][previous]) * condprob[current][current] +
			(condprob[current][previous] + p3_to_10*condprob[(current-2)&3][previous]) * condprob[(current+1)&3][current] +
			(condprob[(current+1)&3][previous] + p3_to_10*condprob[(current-1)&3][previous]) * condprob[(current+2)&3][current] +
			(condprob[(current+2)&3][previous] + p3_to_10*condprob[current][previous]) * condprob[(current+3)&3][current];
		double PX1 =//遷移語が0 or 2 or 3(*p3_to_11)になる組み合わせ
			(condprob[current][previous] + condprob[(current-2)&3][previous] + p3_to_11*condprob[(current-3)&3][previous]) * condprob[current][current] +
			(condprob[(current+1)&3][previous] + condprob[(current-1)&3][previous] + p3_to_11*condprob[(current-2)&3][previous]) * condprob[(current+1)&3][current] +
			(condprob[(current+2)&3][previous] + condprob[current][previous] + p3_to_11*condprob[(current-1)&3][previous]) * condprob[(current+2)&3][current] +
			(condprob[(current+3)&3][previous] + condprob[(current+1)&3][previous] + p3_to_11*condprob[current][previous]) * condprob[(current+3)&3][current];

		previous = current;

		LLRr[j++] = log(P0X/P1X);
		LLRr[j++] = log(PX0/PX1);
	}
}

}

#endif