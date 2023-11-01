//Additive White Gaussian Noise
#ifndef INCLUDE_GUARD_ldpc_AWGN
#define INCLUDE_GUARD_ldpc_AWGN

#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <random>

namespace channel {

template<typename F>
class AWGN {
};

template<>
class AWGN<double> {
	std::mt19937_64 mt;
	std::normal_distribution<double> norm;
	const double llrcoefficient;
public:
	explicit AWGN(double sigmasq, std::int64_t seed);
	explicit AWGN(double sigmasq): AWGN(sigmasq, 0){}

	template<std::size_t L>
	auto noise(const std::bitset<L> &in);
	template<std::size_t L>
	auto LLR(const std::array<double,L> &y) const;
	template<std::size_t L>
	auto estimate(const std::array<double,L> &y) const;
};

template<>
class AWGN<float> {
	std::mt19937_64 mt;
	std::normal_distribution<double> norm;
	const float llrcoefficient;
public:
	explicit AWGN(double sigmasq, std::int64_t seed);
	explicit AWGN(double sigmasq): AWGN(sigmasq, 0){}

	template<std::size_t L>
	auto noise(const std::bitset<L> &in);
	template<std::size_t L>
	auto LLR(const std::array<float,L> &y) const;
	template<std::size_t L>
	auto estimate(const std::array<float,L> &y) const;
};

AWGN<double>::AWGN(double sigmasq, std::int64_t seed):
	mt(seed),
	norm(0.0, std::abs(std::sqrt(sigmasq))),
	llrcoefficient(2.0/sigmasq)
{}

AWGN<float>::AWGN(double sigmasq, std::int64_t seed):
	mt(seed),
	norm(0.0, std::abs(std::sqrt(sigmasq))),
	llrcoefficient(static_cast<float>(2.0/sigmasq))
{}

template<std::size_t L>
auto AWGN<double>::noise(const std::bitset<L> &in){
	std::array<double,L> out;
	for(std::size_t i=0u, iend=in.size(); i<iend; ++i) out[i] = (in.test(i)?-1.0:1.0)+norm(mt);
	return out;
}

template<std::size_t L>
auto AWGN<float>::noise(const std::bitset<L> &in){
	std::array<float,L> out;
	for(std::size_t i=0u, iend=in.size(); i<iend; ++i) out[i] = static_cast<float>((in.test(i)?-1.0:1.0)+norm(mt));
	return out;
}

template<std::size_t L>
auto AWGN<double>::LLR(const std::array<double,L> &y) const{
	std::array<double,L> LLR;
	for(std::size_t i=0u, iend=y.size(); i<iend; ++i) LLR[i] = y[i]*llrcoefficient;
	return LLR;
}

template<std::size_t L>
auto AWGN<float>::LLR(const std::array<float,L> &y) const{
	std::array<float,L> LLR;
	for(std::size_t i=0u, iend=y.size(); i<iend; ++i) LLR[i] = y[i]*llrcoefficient;
	return LLR;
}

template<std::size_t L>
auto AWGN<double>::estimate(const std::array<double,L> &y) const{
	std::bitset<L> est;
	for(std::size_t i=0u, iend=y.size(); i<iend; ++i) est[i] = y[i]<0.0;
	return est;
}

template<std::size_t L>
auto AWGN<float>::estimate(const std::array<float,L> &y) const{
	std::bitset<L> est;
	for(std::size_t i=0u, iend=y.size(); i<iend; ++i) est[i] = y[i]<0.0f;
	return est;
}

}

#endif