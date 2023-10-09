#ifndef __channel_SEQUENCER__
#define __channel_SEQUENCER__

#include <iostream>
#include <vector>
#include <array>
#include <cstdint>
#include <cmath>
#include <random>
#include "DNASmytype.hpp"

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

		void noise(std::vector<code::DNAS::nucleotide_t> &c);

		void message_LLR(const std::vector<code::DNAS::nucleotide_t> &cm, std::vector<double> &LLRm) const;
		void redundancy_LLR(const std::vector<code::DNAS::nucleotide_t> &cr, std::vector<double> &LLRr, code::DNAS::nucleotide_t initial_state) const;
	};
}

#endif