#include "codeIDEC.hpp"

using std::vector, std::string;
using std::int32_t, std::size_t, std::uint32_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;
using code::Marker;

////////////////////////////////////////////////////////////////
//                                                            //
//                        Marker_code                         //
//                                                            //
////////////////////////////////////////////////////////////////

Marker::Marker(const size_t sourcelength, const vector<bool> &marker, const size_t interval, const int32_t drift_range):
	marker(std::make_shared<const IDEC::MarkerCode>(marker, interval, sourcelength)),
	encoder(new IDEC::Marker_Encoding(this->marker)),
	driftrange(drift_range)
{}

vector<bool> Marker::encode(const vector<bool> &source) const{
	return encoder->insert(source);
}

vector<bool> Marker::decode(const vector<bool> &y, const double rate_ins, const double rate_del, const double rate_sub, const std::unique_ptr<IDEC::BCJR_IDEC> &decoder) const{
	auto Pr0 = encoder->Pr0_init();
	vector LLR(marker->codesize(), 0.0);
	decoder->decode_init(y);

	decoder->iterate(LLR, y, Pr0, rate_ins, rate_del, rate_sub);

	return decoder->estimate(encoder->extract(LLR));
}
