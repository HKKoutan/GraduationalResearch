#include "codeLDPCIDEC.hpp"

using std::vector, std::string;
using std::size_t, std::int32_t, std::uint64_t;
using std::cerr, std::endl;
using code::Marker_LDPC_SPD, code::Marker_LDPC_MPD;
using code::LDPC::fptype;

vector<fptype> code::dfc(const vector<double> &in){
	vector<fptype> out(in.size());
	for(size_t i=0u, iend=in.size(); i<iend; ++i) out[i]=static_cast<fptype>(in[i]);
	return out;
}

vector<double> code::fdc(const vector<fptype> &in){
	vector<double> out(in.size());
	for(size_t i=0u, iend=in.size(); i<iend; ++i) out[i]=static_cast<double>(in[i]);
	return out;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                      Marker_LDPC_SPD                       //
//                                                            //
////////////////////////////////////////////////////////////////

Marker_LDPC_SPD::Marker_LDPC_SPD(const string &path, const vector<bool> &marker, const size_t interval, const int32_t drift_range, uint64_t iterationlimit):
	H(std::make_shared<const LDPC::CheckMatrix>(LDPC::CheckMatrix(path))),
	marker(std::make_shared<const IDEC::MarkerCode>(marker, interval, H->codesize())),
	ldpcenc(new LDPC::Generation_Matrix_Encoding(*H)),
	markerenc(new IDEC::Marker_Encoding(this->marker)),
	driftrange(drift_range),
	iterationlimit(iterationlimit)
{}

vector<bool> Marker_LDPC_SPD::encode(const vector<bool> &information) const{
	auto ldpccode = information;
	auto parity = ldpcenc->systematic_parity(information);
	ldpccode.insert(ldpccode.end(), parity.begin(), parity.end());
	return markerenc->insert(ldpcenc->inverse_substitution(ldpccode));
}

vector<bool> Marker_LDPC_SPD::decode(const vector<bool> &y, const double rate_ins, const double rate_del, const double rate_sub, const std::pair<std::unique_ptr<LDPC::I_LDPC_Decoding>, std::unique_ptr<IDEC::BCJR_IDEC>> &decoder){
	const auto& [ldpcdec, markerdec] = decoder;

	auto Pr0 = markerenc->Pr0_init();
	vector<double> LLR(marker->codesize());
	markerdec->decode_init(y);

	markerdec->iterate(LLR, y, Pr0, rate_ins, rate_del, rate_sub);

	const auto QLLR = dfc(markerenc->extract(LLR));
	vector<fptype> QLPR(H->codesize());//列ごとのalphaの和+LLR

	ldpcdec->decode_init();
	for(auto iter=0ui64; !ldpcdec->iterate(QLPR, QLLR) && iter<iterationlimit; ++iter);

	return ldpcdec->estimate(ldpcenc->substitution(QLPR));
}

////////////////////////////////////////////////////////////////
//                                                            //
//                      Marker_LDPC_MPD                       //
//                                                            //
////////////////////////////////////////////////////////////////

Marker_LDPC_MPD::Marker_LDPC_MPD(const string &path, const vector<bool> &marker, const size_t interval, const int32_t drift_range, uint64_t iterationlimit):
	H(std::make_shared<const LDPC::CheckMatrix>(LDPC::CheckMatrix(path))),
	marker(std::make_shared<const IDEC::MarkerCode>(marker, interval, H->codesize())),
	ldpcenc(new LDPC::Generation_Matrix_Encoding(*H)),
	markerenc(new IDEC::Marker_Encoding(this->marker)),
	driftrange(drift_range),
	iterationlimit(iterationlimit)
{}

vector<bool> Marker_LDPC_MPD::encode(const vector<bool> &information) const{
	auto ldpccode = information;
	auto parity = ldpcenc->systematic_parity(information);
	ldpccode.insert(ldpccode.end(), parity.begin(), parity.end());
	return markerenc->insert(ldpcenc->inverse_substitution(ldpccode));
}

vector<bool> Marker_LDPC_MPD::decode(const vector<bool> &y, const double rate_ins, const double rate_del, const double rate_sub, const std::pair<std::unique_ptr<LDPC::I_LDPC_Decoding>, std::unique_ptr<IDEC::BCJR_IDEC>> &decoder){
	const auto& [ldpcdec, markerdec] = decoder;
	vector<double> LLR(marker->codesize());
	vector<fptype> QLLR(H->codesize());
	vector<fptype> QLPR(H->codesize());//列ごとのalphaの和+LLR

	markerdec->decode_init(y);
	ldpcdec->decode_init();
	auto iter = 0ui64;
	do{
		for(size_t i=0u, iend=QLPR.size(); i<iend; ++i) QLPR[i] -= QLLR[i];//対数事後確率比-対数事前確率比=対数尤度比
		auto Pr0 = markerenc->Pr0_init(fdc(QLPR));
		markerdec->iterate(LLR, y, Pr0, rate_ins, rate_del, rate_sub);
		QLLR = dfc(markerenc->extract(LLR));
	}while(!ldpcdec->iterate(QLPR, QLLR) && ++iter<iterationlimit);

	return ldpcdec->estimate(ldpcenc->substitution(QLPR));
}
