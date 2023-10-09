#ifndef __code_CODEDNAS__
#define __code_CODEDNAS__

#include <iostream>
#include <cstdint>
#include <vector>
#include <exception>
#include "DNASmytype.hpp"

namespace code{
	namespace DNAS {
		nucleotide_t VLRLL_encode(const std::vector<bool> &source, std::vector<nucleotide_t> &code, const nucleotide_t initial_state = 0);
		std::vector<bool>::iterator VLRLL_encode(const std::vector<bool>::iterator &source_begin, const std::vector<bool>::iterator &source_end, const std::vector<nucleotide_t>::iterator &code_begin, const std::vector<nucleotide_t>::iterator &code_end, const nucleotide_t initial_state = 0);
		nucleotide_t modified_VLRLL_encode(const std::vector<bool> &source, std::vector<nucleotide_t> &code, const nucleotide_t initial_state=0);

		void interim_map(const std::vector<nucleotide_t> &source, std::vector<bool> &code);
		void interim_demap(const std::vector<bool> &source, std::vector<nucleotide_t> &code);

		void VLRLL_decode(const std::vector<nucleotide_t> &source, std::vector<bool> &decode, const nucleotide_t initial_state = 0);
		void modified_VLRLL_decode(const std::vector<nucleotide_t> &source, std::vector<bool> &decode, const nucleotide_t initial_state = 0);

		void nt_addequalizing_encode(const std::vector<nucleotide_t> &cr, std::vector<nucleotide_t> &crbar, std::vector<bool> &info, std::uint32_t qty_AT=0, std::uint32_t qty_GC=0);
		void nt_addequalizing_decode(const std::vector<nucleotide_t> &crbar, const std::vector<bool> &i, std::vector<nucleotide_t> &cr);
		void nt_qty_count(const std::vector<nucleotide_t> &c, std::uint32_t &qty_AT, std::uint32_t &qty_GC);
	}
}

#endif