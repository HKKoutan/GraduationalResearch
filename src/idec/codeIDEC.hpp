#ifndef __code_CODEIDEC__
#define __code_CODEIDEC__

#include "IDECbase.hpp"

namespace code {
	class Marker{//Marker code
		std::shared_ptr<const IDEC::MarkerCode> marker;
		std::unique_ptr<IDEC::Marker_Encoding> encoder;
		std::int32_t driftrange;
	public:
		Marker(const std::size_t sourcelength, const std::vector<bool> &marker, const std::size_t interval, const std::int32_t drift_range);

		std::vector<bool> encode(const std::vector<bool> &source) const;
		std::vector<bool> decode(const std::vector<bool> &y, const double rate_ins, const double rate_del, const double rate_sub, const std::unique_ptr<IDEC::BCJR_IDEC> &decoder) const;

		const auto create_decoder(){return std::unique_ptr<IDEC::BCJR_IDEC>(new IDEC::BCJR_IDEC(marker->codesize(), driftrange));}
	};
}

#endif