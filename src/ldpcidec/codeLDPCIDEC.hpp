#ifndef __code_CODELDPCIDEC__
#define __code_CODELDPCIDEC__

#include "../ldpc/LDPCbase.hpp"
#include "../idec/IDECbase.hpp"

namespace code{
	std::vector<LDPC::fptype> dfc(const std::vector<double> &in);
	std::vector<double> fdc(const std::vector<LDPC::fptype> &in);

	class Marker_LDPC_SPD {
		static constexpr std::uint64_t DEFAULT_ITERATION_LIMIT = 100;

		std::shared_ptr<const LDPC::CheckMatrix> H;
		std::shared_ptr<const IDEC::MarkerCode> marker;
		std::unique_ptr<const LDPC::I_LDPC_Encoding> ldpcenc;
		std::unique_ptr<const IDEC::Marker_Encoding> markerenc;
		std::int32_t driftrange;
		std::uint64_t iterationlimit;
	public:
		Marker_LDPC_SPD(const std::string &path, const std::vector<bool> &marker, const std::size_t interval, const std::int32_t drift_range, std::uint64_t iterationlimit);
		Marker_LDPC_SPD(const std::string &path, const std::vector<bool> &marker, const std::size_t interval, const std::int32_t drift_range): Marker_LDPC_SPD(path, marker, interval, drift_range, DEFAULT_ITERATION_LIMIT){}

		std::vector<bool> encode(const std::vector<bool> &information) const;
		std::vector<bool> decode(const std::vector<bool> &y, const double rate_ins, const double rate_del, const double rate_sub, const std::pair<std::unique_ptr<LDPC::I_LDPC_Decoding>, std::unique_ptr<IDEC::BCJR_IDEC>> &decoder);

		template<typename T> const auto create_decoder(){return std::make_pair(std::unique_ptr<LDPC::I_LDPC_Decoding>(new T(H)), std::unique_ptr<IDEC::BCJR_IDEC>(new IDEC::BCJR_IDEC(marker->codesize(), driftrange)));}

		auto sourcesize() const{return H->sourcesize();}
		auto codesize() const{return marker->codesize();}
	};

	class Marker_LDPC_MPD {
		static constexpr std::uint64_t DEFAULT_ITERATION_LIMIT = 100;

		std::shared_ptr<const LDPC::CheckMatrix> H;
		std::shared_ptr<const IDEC::MarkerCode> marker;
		std::unique_ptr<const LDPC::I_LDPC_Encoding> ldpcenc;
		std::unique_ptr<const IDEC::Marker_Encoding> markerenc;
		std::int32_t driftrange;
		std::uint64_t iterationlimit;
	public:
		Marker_LDPC_MPD(const std::string &path, const std::vector<bool> &marker, const std::size_t interval, const std::int32_t drift_range, std::uint64_t iterationlimit);
		Marker_LDPC_MPD(const std::string &path, const std::vector<bool> &marker, const std::size_t interval, const std::int32_t drift_range): Marker_LDPC_MPD(path, marker, interval, drift_range, DEFAULT_ITERATION_LIMIT){}

		std::vector<bool> encode(const std::vector<bool> &information) const;
		std::vector<bool> decode(const std::vector<bool> &y, const double rate_ins, const double rate_del, const double rate_sub, const std::pair<std::unique_ptr<LDPC::I_LDPC_Decoding>, std::unique_ptr<IDEC::BCJR_IDEC>> &decoder);

		template<typename T> const auto create_decoder(){return std::make_pair(std::unique_ptr<LDPC::I_LDPC_Decoding>(new T(H)), std::unique_ptr<IDEC::BCJR_IDEC>(new IDEC::BCJR_IDEC(marker->codesize(), driftrange)));}

		auto sourcesize() const{return H->sourcesize();}
		auto codesize() const{return marker->codesize();}
	};
}

#endif