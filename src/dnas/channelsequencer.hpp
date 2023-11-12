#ifndef INCLUDE_GUARD_dnas_channelsequencer
#define INCLUDE_GUARD_dnas_channelsequencer

#include <array>
#include <cmath>
#include <random>
#include <concepts>
#include <exception>
#include "DNASnttype.hpp"

namespace channel{

template<std::uint8_t ATGC=0x1B>
class Nanopore_Sequencing{
	std::mt19937_64 mt;
	std::uniform_real_distribution<> uniform;
	static constexpr std::uint8_t nA = (ATGC>>6)&3, nT = (ATGC>>4)&3, nG = (ATGC>>2)&3, nC = ATGC&3;
	static_assert(nA!=nT && nA!=nG && nA!=nC && nT!=nG && nT!=nC && nG!=nC);//ATGCに重複がない
	static constexpr double A0 = (nA>>1?-1.0:1.0), A1 = (nA&1?-1.0:1.0), T0 = (nT>>1?-1.0:1.0), T1 = (nT&1?-1.0:1.0), G0 = (nG>>1?-1.0:1.0), G1 = (nG&1?-1.0:1.0), C0 = (nC>>1?-1.0:1.0), C1 = (nC&1?-1.0:1.0); 
	const std::array<double,4> error_rate;//error_rate[0]~error_rate[3] = p1~p4
	const std::array<double,4> non_error_rate;//ATGCそれぞれ誤らない確率
	const std::array<std::array<double,4>,4> condprob;//condprob[a][b] P(a|b) 条件付き確率
	auto condprob_init() const;

public:
	Nanopore_Sequencing(double alpha, std::int64_t seed);
	Nanopore_Sequencing(double alpha);

	template<std::size_t L>
	auto noise(const std::array<code::DNAS::nucleotide_t<ATGC>,L> &in);
	template<std::size_t L>
	auto likelihood(const std::array<code::DNAS::nucleotide_t<ATGC>,L> &in);

	// template<std::size_t S, std::size_t R>
	// auto VLRLL_LLR(const std::array<code::DNAS::nucleotide_t<ATGC>,S> &cm, const std::array<code::DNAS::nucleotide_t<ATGC>,R> &cr, code::DNAS::nucleotide_t<ATGC> initial_state) const;
	// template<std::size_t S>
	// auto differential_LLR(const std::array<code::DNAS::nucleotide_t<ATGC>,S> &code, code::DNAS::nucleotide_t<ATGC> initial_state = 0) const;
	// template<std::size_t BS = 0, std::size_t S>
	// auto division_balancing_LLR(const std::array<code::DNAS::nucleotide_t<ATGC>,S> &code, code::DNAS::nucleotide_t<ATGC> initial_state = 0) const;
};

template<>
auto Nanopore_Sequencing<0x1B>::condprob_init() const{
	return std::array<std::array<double,4>,4>{
		non_error_rate[0],error_rate[1],error_rate[3],error_rate[2],
		error_rate[1],non_error_rate[1],error_rate[2],error_rate[0],
		error_rate[3],error_rate[2],non_error_rate[2],error_rate[1],
		error_rate[2],error_rate[0],error_rate[1],non_error_rate[3]
	};
}

template<>
auto Nanopore_Sequencing<0x27>::condprob_init() const{
	return std::array<std::array<double,4>,4>{
		non_error_rate[0],error_rate[3],error_rate[1],error_rate[2],
		error_rate[3],non_error_rate[2],error_rate[2],error_rate[1],
		error_rate[1],error_rate[2],non_error_rate[1],error_rate[0],
		error_rate[2],error_rate[1],error_rate[0],non_error_rate[3]
	};
}

template<std::uint8_t ATGC>
Nanopore_Sequencing<ATGC>::Nanopore_Sequencing(double alpha) : Nanopore_Sequencing(alpha, 0){}

template<std::uint8_t ATGC>
Nanopore_Sequencing<ATGC>::Nanopore_Sequencing(double alpha, std::int64_t seed):
	mt(seed),
	uniform(0,1),
	error_rate{4*alpha, alpha, 0.01, 0},
	non_error_rate{1.0-error_rate[1]-error_rate[2]-error_rate[3], 1.0-error_rate[0]-error_rate[1]-error_rate[2], 1.0-error_rate[1]-error_rate[2]-error_rate[3], 1.0-error_rate[0]-error_rate[1]-error_rate[2]},//ATGCそれぞれの誤らない確率
	condprob(condprob_init())
{
	static_assert((nA!=nT)&&(nA!=nG)&&(nA!=nC)&&(nT!=nG)&&(nT!=nC)&&(nG!=nC));
	if(alpha>0.198||alpha<0) throw std::out_of_range("alpha out of range");
}

template<std::uint8_t ATGC>
template<std::size_t L>
auto Nanopore_Sequencing<ATGC>::noise(const std::array<code::DNAS::nucleotide_t<ATGC>,L> &in){
	auto out = in;
	for(auto &i: out){
		code::DNAS::nucleotide_t<ATGC> j=0;
		double rand = uniform(mt)-condprob[i][j];
		while(rand>=0){
			j+=1;
			rand-=condprob[i][j];
		}
		i=j;
	}
	return out;
}

template<std::uint8_t ATGC>
template<std::size_t L>
auto Nanopore_Sequencing<ATGC>::likelihood(const std::array<code::DNAS::nucleotide_t<ATGC>,L> &in){
	std::array<code::DNAS::nucleotide_p<ATGC>,L> likelihoods;
	for(std::size_t i=0; i<L; ++i){
		auto &li = likelihoods[i];
		auto &ii = in[i];
		for(std::uint8_t n=0; n<4; ++n) li.lh(code::DNAS::nucleotide_t<ATGC>(n)) = condprob[ii][n];
	}
	return likelihoods;
}

// template<std::uint8_t ATGC>
// template<std::size_t S, std::size_t R>
// auto Nanopore_Sequencing<ATGC,float>::VLRLL_LLR(const std::array<code::DNAS::nucleotide_t<ATGC>,S> &cm, const std::array<code::DNAS::nucleotide_t<ATGC>,R> &cr, code::DNAS::nucleotide_t<ATGC> initial_state) const{
// 	std::array<float,S*2+R*2> LLR;
// 	std::size_t j=0u;
// 	//情報部
// 	for(const auto &current: cm){
// 		switch(current){
// 		case nA:
// 			LLR[j++] = static_cast<float>(A0*std::log(non_error_rate[0]+error_rate[1])-std::log(error_rate[3]+error_rate[2]));
// 			LLR[j++] = static_cast<float>(A1*std::log(non_error_rate[0]+error_rate[3])-std::log(error_rate[1]+error_rate[2]));
// 			break;

// 		case nT:
// 			LLR[j++] = static_cast<float>(T0*std::log(non_error_rate[1]+error_rate[1])-std::log(error_rate[0]+error_rate[2]));
// 			LLR[j++] = static_cast<float>(T1*std::log(non_error_rate[1]+error_rate[0])-std::log(error_rate[1]+error_rate[2]));
// 			break;

// 		case nG:
// 			LLR[j++] = static_cast<float>(G0*std::log(non_error_rate[2]+error_rate[1])-std::log(error_rate[3]+error_rate[2]));
// 			LLR[j++] = static_cast<float>(G1*std::log(non_error_rate[2]+error_rate[3])-std::log(error_rate[1]+error_rate[2]));
// 			break;

// 		case nC:
// 			LLR[j++] = static_cast<float>(C0*std::log(non_error_rate[3]+error_rate[1])-std::log(error_rate[0]+error_rate[2]));
// 			LLR[j++] = static_cast<float>(C1*std::log(non_error_rate[3]+error_rate[0])-std::log(error_rate[1]+error_rate[2]));
// 			break;
// 		}
// 	}
// 	//冗長部
// 	constexpr double p3_to_11 = 1.0/22.0;
// 	constexpr double p3_to_10 = 21.0/22.0;

// 	for(auto previous = initial_state; const auto &current: cr){
// 		double P0X =//遷移語が1 or 2になる組み合わせ
// 			condprob[previous][previous] * (condprob[previous+1][current] + condprob[previous+2][current]) +
// 			condprob[previous-1][previous] * (condprob[previous][current] + condprob[previous+1][current]) +
// 			condprob[previous-2][previous] * (condprob[previous-1][current] + condprob[previous][current]) +
// 			condprob[previous-3][previous] * (condprob[previous-2][current] + condprob[previous-1][current]);
// 		double P1X =//遷移語が0 or 3になる組み合わせ
// 			condprob[previous][previous] * (condprob[previous][current] + condprob[previous+3][current]) +
// 			condprob[previous-1][previous] * (condprob[previous-1][current] + condprob[previous+2][current]) +
// 			condprob[previous-2][previous] * (condprob[previous-2][current] + condprob[previous+1][current]) +
// 			condprob[previous-3][previous] * (condprob[previous-3][current] + condprob[previous][current]);
// 		double PX0 =//遷移語が1 or 3(*p3_to_10)になる組み合わせ
// 			(condprob[current-1][previous] + p3_to_10*condprob[current-3][previous]) * condprob[current][current] +
// 			(condprob[current][previous] + p3_to_10*condprob[current-2][previous]) * condprob[current+1][current] +
// 			(condprob[current+1][previous] + p3_to_10*condprob[current-1][previous]) * condprob[current+2][current] +
// 			(condprob[current+2][previous] + p3_to_10*condprob[current][previous]) * condprob[current+3][current];
// 		double PX1 =//遷移語が0 or 2 or 3(*p3_to_11)になる組み合わせ
// 			(condprob[current][previous] + condprob[current-2][previous] + p3_to_11*condprob[current-3][previous]) * condprob[current][current] +
// 			(condprob[current+1][previous] + condprob[current-1][previous] + p3_to_11*condprob[current-2][previous]) * condprob[current+1][current] +
// 			(condprob[current+2][previous] + condprob[current][previous] + p3_to_11*condprob[current-1][previous]) * condprob[current+2][current] +
// 			(condprob[current+3][previous] + condprob[current+1][previous] + p3_to_11*condprob[current][previous]) * condprob[current+3][current];

// 		previous = current;
// 		LLR[j++] = static_cast<float>(std::log(P0X)-std::log(P1X));
// 		LLR[j++] = static_cast<float>(std::log(PX0)-std::log(PX1));
// 	}
// 	return LLR;
// }

// template<std::uint8_t ATGC>
// template<std::size_t S>
// auto Nanopore_Sequencing<ATGC,float>::differential_LLR(const std::array<code::DNAS::nucleotide_t<ATGC>,S> &code, code::DNAS::nucleotide_t<ATGC> initial_state) const{
// 	std::array<float,S*2> LLR;
// 	auto previous = initial_state;

// 	for(std::size_t j=0u; const auto &current: code){
// 		double P0X =//遷移語が0 or 1になる組み合わせ
// 			condprob[previous][previous] * (condprob[previous][current] + condprob[previous+1][current]) +
// 			condprob[previous-1][previous] * (condprob[previous-1][current] + condprob[previous][current]) +
// 			condprob[previous-2][previous] * (condprob[previous-2][current] + condprob[previous-1][current]) +
// 			condprob[previous-3][previous] * (condprob[previous-3][current] + condprob[previous-2][current]);
// 		double P1X =//遷移語が2 or 3になる組み合わせ
// 			condprob[previous][previous] * (condprob[previous+2][current] + condprob[previous+3][current]) +
// 			condprob[previous-1][previous] * (condprob[previous+1][current] + condprob[previous+2][current]) +
// 			condprob[previous-2][previous] * (condprob[previous][current] + condprob[previous+1][current]) +
// 			condprob[previous-3][previous] * (condprob[previous-1][current] + condprob[previous][current]);
// 		double PX0 =//遷移語が0 or 2になる組み合わせ
// 			condprob[previous][previous] * (condprob[previous][current] + condprob[previous+2][current]) +
// 			condprob[previous-1][previous] * (condprob[previous-1][current] + condprob[previous+1][current]) +
// 			condprob[previous-2][previous] * (condprob[previous-2][current] + condprob[previous][current]) +
// 			condprob[previous-3][previous] * (condprob[previous-3][current] + condprob[previous-1][current]);
// 		double PX1 =//遷移語が1 or 3になる組み合わせ
// 			condprob[previous][previous] * (condprob[previous+1][current] + condprob[previous+3][current]) +
// 			condprob[previous-1][previous] * (condprob[previous][current] + condprob[previous+2][current]) +
// 			condprob[previous-2][previous] * (condprob[previous-1][current] + condprob[previous+1][current]) +
// 			condprob[previous-3][previous] * (condprob[previous-2][current] + condprob[previous][current]);

// 		previous = current;
// 		LLR[j++] = static_cast<float>(std::log(P0X)-std::log(P1X));
// 		LLR[j++] = static_cast<float>(std::log(PX0)-std::log(PX1));
// 	}
// 	return LLR;
// }

// template<>
// template<std::size_t BS, std::size_t S>
// auto Nanopore_Sequencing<0x1B,float>::division_balancing_LLR(const std::array<code::DNAS::nucleotide_t<0x1B>,S> &code, code::DNAS::nucleotide_t<0x1B> initial_state) const{
// 	constexpr std::size_t block_size = BS==0?S:BS;
// 	constexpr std::size_t div_size = block_size>>1;
// 	constexpr double div_prob = 1.0/div_size;
// 	constexpr double non_div_prob = 1.0-div_prob;
// 	static_assert(block_size%2==0&&S%block_size==0);
// 	std::array<float,S*2> LLR;
// 	auto previous = initial_state;

// 	for(std::size_t i=0u, iend=S/block_size, k=0u; i<iend; ++i){
// 		std::size_t block_head = i*block_size, block_tail = block_head+block_size;
// 		for(std::size_t j=block_head; j<block_tail; ++j){
// 			auto current = code[j];
// 			double P0X =//遷移語が0 or 1になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-2][current-2]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-1][current-2]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-3][current-2]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-2][current-2]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous][current-2]*div_prob) +
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-3][current-2]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-3][current]*non_div_prob + condprob[previous-1][current-2]*div_prob) +
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous][current-2]*div_prob));
// 			double P1X =//遷移語が2 or 3になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous][current-2]*div_prob) +
// 					(condprob[previous+3][current]*non_div_prob + condprob[previous+1][current-2]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-1][current-2]*div_prob) +
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous][current-2]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-2][current-2]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-1][current-2]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-3][current-2]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-2][current-2]*div_prob));
// 			double PX0 =//遷移語が0 or 2になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-2][current-2]*div_prob) +
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous][current-2]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-3][current-2]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-1][current-2]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous][current-2]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-2][current-2]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-3][current]*non_div_prob + condprob[previous-1][current-2]*div_prob) +
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-3][current-2]*div_prob));
// 			double PX1 =//遷移語が1 or 3になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-1][current-2]*div_prob) +
// 					(condprob[previous+3][current]*non_div_prob + condprob[previous+1][current-2]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-2][current-2]*div_prob) +
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous][current-2]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-3][current-2]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-1][current-2]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous][current-2]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-2][current-2]*div_prob));

// 			previous = current;
// 			LLR[k++] = static_cast<float>(std::log(P0X)-std::log(P1X));
// 			LLR[k++] = static_cast<float>(std::log(PX0)-std::log(PX1));			
// 		}
// 	}
// 	return LLR;
// }

// template<>
// template<std::size_t BS, std::size_t S>
// auto Nanopore_Sequencing<0x27,float>::division_balancing_LLR(const std::array<code::DNAS::nucleotide_t<0x27>,S> &code, code::DNAS::nucleotide_t<0x27> initial_state) const{
// 	constexpr std::size_t block_size = BS==0?S:BS;
// 	constexpr std::size_t div_size = block_size>>1;
// 	constexpr double div_prob = 1.0/div_size;
// 	constexpr double non_div_prob = 1.0-div_prob;
// 	static_assert(block_size%2==0&&S%block_size==0);
// 	std::array<float,S*2> LLR;
// 	auto previous = initial_state;

// 	for(std::size_t i=0u, iend=S/block_size, k=0u; i<iend; ++i){
// 		std::size_t block_head = i*block_size, block_tail = block_head+block_size;
// 		for(std::size_t j=block_head, jend=block_head+div_size; j<jend; ++j){
// 			auto current = code[j];
// 			double P0X =//遷移語が0 or 1になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-1][current-1]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous][current-1]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-2][current-1]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-1][current-1]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous-3][current-1]*div_prob) +
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-2][current-1]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-3][current]*non_div_prob + condprob[previous][current-1]*div_prob) +
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous-3][current-1]*div_prob));
// 			double P1X =//遷移語が2 or 3になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous+1][current-1]*div_prob) +
// 					(condprob[previous+3][current]*non_div_prob + condprob[previous+2][current-1]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous][current-1]*div_prob) +
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous+1][current-1]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-1][current-1]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous][current-1]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-2][current-1]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-1][current-1]*div_prob));
// 			double PX0 =//遷移語が0 or 2になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-1][current-1]*div_prob) +
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous+1][current-1]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-2][current-1]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous][current-1]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous-3][current-1]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-1][current-1]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-3][current]*non_div_prob + condprob[previous][current-1]*div_prob) +
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-2][current-1]*div_prob));
// 			double PX1 =//遷移語が1 or 3になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous][current-1]*div_prob) +
// 					(condprob[previous+3][current]*non_div_prob + condprob[previous+2][current-1]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-1][current-1]*div_prob) +
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous+1][current-1]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous-2][current-1]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous][current-1]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous-3][current-1]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-1][current-1]*div_prob));

// 			previous = current;
// 			LLR[k++] = static_cast<float>(std::log(P0X)-std::log(P1X));
// 			LLR[k++] = static_cast<float>(std::log(PX0)-std::log(PX1));			
// 		}
// 		for(std::size_t j=block_head+div_size; j<block_tail; ++j){
// 			auto current = code[j];
// 			double P0X =//遷移語が0 or 1になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-3][current-3]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-2][current-3]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous][current-3]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-3][current-3]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous-1][current-3]*div_prob) +
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous][current-3]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-3][current]*non_div_prob + condprob[previous-2][current-3]*div_prob) +
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous-1][current-3]*div_prob));
// 			double P1X =//遷移語が2 or 3になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous-1][current-3]*div_prob) +
// 					(condprob[previous+3][current]*non_div_prob + condprob[previous][current-3]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-2][current-3]*div_prob) +
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous-1][current-3]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-3][current-3]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-2][current-3]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous][current-3]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-3][current-3]*div_prob));
// 			double PX0 =//遷移語が0 or 2になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-3][current-3]*div_prob) +
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous-1][current-3]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous][current-3]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-2][current-3]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous-1][current-3]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-3][current-3]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-3][current]*non_div_prob + condprob[previous-2][current-3]*div_prob) +
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous][current-3]*div_prob));
// 			double PX1 =//遷移語が1 or 3になる組み合わせ
// 				condprob[previous][previous] * (
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-2][current-3]*div_prob) +
// 					(condprob[previous+3][current]*non_div_prob + condprob[previous][current-3]*div_prob)) +
// 				condprob[previous-1][previous] * (
// 					(condprob[previous][current]*non_div_prob + condprob[previous-3][current-3]*div_prob) +
// 					(condprob[previous+2][current]*non_div_prob + condprob[previous-1][current-3]*div_prob)) +
// 				condprob[previous-2][previous] * (
// 					(condprob[previous-1][current]*non_div_prob + condprob[previous][current-3]*div_prob) +
// 					(condprob[previous+1][current]*non_div_prob + condprob[previous-2][current-3]*div_prob)) +
// 				condprob[previous-3][previous] * (
// 					(condprob[previous-2][current]*non_div_prob + condprob[previous-1][current-3]*div_prob) +
// 					(condprob[previous][current]*non_div_prob + condprob[previous-3][current-3]*div_prob));

// 			previous = current;
// 			LLR[k++] = static_cast<float>(std::log(P0X)-std::log(P1X));
// 			LLR[k++] = static_cast<float>(std::log(PX0)-std::log(PX1));
// 		}
// 	}
// 	return LLR;
// }

}

#endif