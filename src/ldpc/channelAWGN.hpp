//Additive White Gaussian Noise
#ifndef INCLUDE_GUARD_ldpc_channelAWGN
#define INCLUDE_GUARD_ldpc_channelAWGN

#include <bitset>
#include <vector>
#include <concepts>
#include <cstdint>
#include <cmath>
#include <random>

namespace channel {

class AWGN {
	std::mt19937_64 mt;
	std::normal_distribution<double> norm;
	const double llrcoefficient;
public:
	explicit AWGN(double sigmasq, std::uint64_t seed);
	explicit AWGN(double sigmasq): AWGN(sigmasq, 0){}

	template<std::floating_point T, std::size_t L>
	auto noise(const std::bitset<L> &in);
	template<std::floating_point T, std::size_t L>
	auto LLR(const std::array<T,L> &y) const;
	template<std::floating_point T, std::size_t L>
	auto estimate(const std::array<T,L> &y) const;
};

AWGN::AWGN(double sigmasq, std::uint64_t seed):
	mt(seed),
	norm(0.0, std::sqrt(std::abs(sigmasq))),
	llrcoefficient(2.0/std::abs(sigmasq))
{}

template<std::floating_point T, std::size_t L>
auto AWGN::noise(const std::bitset<L> &in){
	std::array<T,L> out;
	for(std::size_t i=0; i<L; ++i) out[i] = static_cast<T>((in.test(i)?-1.0:1.0)+norm(mt));
	return out;
}

template<std::floating_point T, std::size_t L>
auto AWGN::LLR(const std::array<T,L> &y) const{
	std::array<T,L> LLR;
	for(std::size_t i=0; i<L; ++i) LLR[i] = y[i]*static_cast<T>(llrcoefficient);
	return LLR;
}

template<std::floating_point T, std::size_t L>
auto AWGN::estimate(const std::array<T,L> &y) const{
	std::bitset<L> est;
	for(std::size_t i=0; i<L; ++i) est[i] = y[i]<0.0;
	return est;
}

}

#endif