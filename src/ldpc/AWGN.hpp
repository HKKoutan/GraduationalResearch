//Additive White Gaussian Noise
#ifndef __channel_AWGN__
#define __channel_AWGN__

#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <random>

namespace channel {

template<typename F, std::size_t L>
class AWGN {
};

template<std::size_t L>
class AWGN<double,L> {
	std::mt19937_64 mt;
	std::normal_distribution<double> norm;
	const double llrcoefficient;
public:
	explicit AWGN(double sigmasq, std::int64_t seed);
	explicit AWGN(double sigmasq): AWGN(sigmasq, 0){}

	auto noise(const std::bitset<L> &in);
	auto LLR(const std::array<double,L> &y) const;
};

template<std::size_t L>
class AWGN<float,L> {
	std::mt19937_64 mt;
	std::normal_distribution<double> norm;
	const float llrcoefficient;
public:
	explicit AWGN(double sigmasq, std::int64_t seed);
	explicit AWGN(double sigmasq): AWGN(sigmasq, 0){}

	auto noise(const std::bitset<L> &in);
	auto LLR(const std::array<float,L> &y) const;
};

template<std::size_t L>
AWGN<double,L>::AWGN(double sigmasq, std::int64_t seed):
	mt(seed),
	norm(0.0, std::abs(std::sqrt(sigmasq))),
	llrcoefficient(2.0/sigmasq)
{}

template< std::size_t L>
AWGN<float,L>::AWGN(double sigmasq, std::int64_t seed):
	mt(seed),
	norm(0.0, std::abs(std::sqrt(sigmasq))),
	llrcoefficient(static_cast<float>(2.0/sigmasq))
{}

template<std::size_t L>
auto AWGN<double,L>::noise(const std::bitset<L> &in){
	std::array<double,L> out(in.size());
	for(std::size_t i=0u, iend=in.size(); i<iend; ++i) out[i] = (in.test(i)?-1.0:1.0)+norm(mt);
	return out;
}

template<std::size_t L>
auto AWGN<float,L>::noise(const std::bitset<L> &in){
	std::array<float,L> out;
	for(std::size_t i=0u, iend=in.size(); i<iend; ++i) out[i] = static_cast<float>((in.test(i)?-1.0:1.0)+norm(mt));
	return out;
}

template<std::size_t L>
auto AWGN<double,L>::LLR(const std::array<double,L> &y) const{
	std::array<double,L> LLR(y.size());
	for(std::size_t i=0u, iend=y.size(); i<iend; ++i) LLR[i] = y[i]*llrcoefficient;
	return LLR;
}

template<std::size_t L>
auto AWGN<float,L>::LLR(const std::array<float,L> &y) const{
	std::array<float,L> LLR;
	for(std::size_t i=0u, iend=y.size(); i<iend; ++i) LLR[i] = y[i]*llrcoefficient;
	return LLR;
}

}

#endif