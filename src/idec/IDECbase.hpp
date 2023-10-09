//Insertion Deletion Error Correction
#ifndef __code_IDECBASE__
#define __code_IDECBASE__

#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <memory>
#include <exception>

namespace code {
	namespace IDEC {
		class MarkerCode {
			const std::vector<bool> marker;
			const std::size_t interval;
			const std::size_t sourcelength;
			const std::size_t codelength;
		public:
			MarkerCode(const std::vector<bool> &marker, std::size_t interval, std::size_t sourcelength);

			auto intervalsize() const noexcept{return interval;}
			auto codesize() const noexcept{return codelength;}
			auto sourcesize() const noexcept{return sourcelength;}
			auto size() const noexcept{return marker.size();}

			const auto &ref_marker() const noexcept{return marker;}
			auto begin() const noexcept{return marker.cbegin();}
			auto cbegin() const noexcept{return marker.cbegin();}
			auto end() const noexcept{return marker.cend();}
			auto cend() const noexcept{return marker.cend();}
		};

		class Marker_Encoding {
			std::shared_ptr<const MarkerCode> marker;

			static double LVR_to_Pr0(double LVR);
		public:
			explicit Marker_Encoding(std::shared_ptr<const MarkerCode> marker): marker(marker){}

			std::vector<double> Pr0_init() const;
			std::vector<double> Pr0_init(const std::vector<double> &LEVR) const;
			std::vector<bool> insert(const std::vector<bool> &source) const;
			std::vector<bool> extract(const std::vector<bool> &code) const;
			std::vector<double> extract(const std::vector<double> &code) const;
		};

		class BCJR_IDEC {//Insertion Deletion Error Correction
			std::size_t length;
			std::int32_t drift_ceil;
			std::int32_t drift_floor;
			std::size_t smax;
			std::vector<std::vector<double>> alpha;
			std::vector<std::vector<double>> beta;
		public:
			BCJR_IDEC(std::size_t length, const std::int32_t drift_ceil, const std::int32_t drift_floor);
			BCJR_IDEC(std::size_t length, const std::int32_t drift_range): BCJR_IDEC(length, drift_range, drift_range){}

			void decode_init(const std::vector<bool> &y);

			void iterate(std::vector<double> &LLR, const std::vector<bool> &y, const std::vector<double> &Pr0, const double rate_ins, const double rate_del, const double rate_sub);
			static std::vector<bool> estimate(const std::vector<double> &LLR);
		};
	}
}

#endif