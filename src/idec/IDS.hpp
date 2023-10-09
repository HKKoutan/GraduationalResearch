//Insertion Deletion Substitution
#ifndef __channel_IDS__
#define __channel_IDS__

#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <random>

namespace channel {
	class IDS {
		std::mt19937_64 mt;
		std::uniform_real_distribution<> uniform;
		const double rate_ins;
		//const double deletion_rate;
		const double rate_indel;
		const double rate_sub;
		const std::int32_t drift_max;

	public:
		explicit IDS(double rate_ins, double rate_del, std::int32_t drift_max, double rate_sub, std::int64_t seed);
		explicit IDS(double rate_ins, double rate_del, double rate_sub, std::int64_t seed) : IDS(rate_ins, rate_del, 0, rate_sub, seed){}
		explicit IDS(double rate_ins, double rate_del, double rate_sub) : IDS(rate_ins, rate_del, 0, rate_sub, 0){}
		explicit IDS(double rate_ins, double rate_del, std::int32_t drift_max): IDS(rate_ins, rate_del, drift_max, 0.0, 0){}
		explicit IDS(double rate_ins, double rate_del): IDS(rate_ins, rate_del, 0, 0.0, 0){}
		explicit IDS(double rate_indel, std::int32_t drift_max, double rate_sub): IDS(rate_indel, rate_indel, drift_max, rate_sub, 0){}
		explicit IDS(double rate_indel, std::int32_t drift_max): IDS(rate_indel, rate_indel, drift_max, 0.0, 0){}

		std::vector<bool> noise(const std::vector<bool> &in);

		double insertion_rate() const{return rate_ins;}
		double deletion_rate() const{return rate_indel-rate_ins;}
		double substitution_rate() const{return rate_sub;}
		std::uint32_t drift_range() const{return drift_max;}
	};
}

#endif