#ifndef INCLUDE_GUARD_ldpc_codeSystematicLDPC
#define INCLUDE_GUARD_ldpc_codeSystematicLDPC

#include "LDPCdecoding.hpp"
#include "LDPCencoding.hpp"

namespace code {

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class SystematicLDPC {
	static constexpr std::uint64_t DEFAULT_ITERATION_LIMIT = 100;

	LDPC::validCheckMatrixType_t<S,C> H;
	LDPC::GenerationMatrix_encoding<decltype(H)> encoder;
	LDPC::Iterative_decoding<decltype(H)> decoder;
	std::uint64_t iterationlimit;//反復回数上限
public:
	struct DecoderType{
		struct SumProduct {using type = typename decltype(decoder)::SumProduct;};
		struct MinSum {using type = typename decltype(decoder)::MinSum;};
	};

	explicit SystematicLDPC(std::uint64_t iterationlimit);
	SystematicLDPC(): SystematicLDPC(DEFAULT_ITERATION_LIMIT){}

	auto encode(const std::bitset<S> &information) const;//引数から符号語を求める
	auto encode_redundancy(const std::bitset<S> &information) const;//引数から冗長を求める
	template<LDPC::DecoderType T, std::floating_point U>
	auto decode(const std::array<U,C> &LLR);

	auto sourcesize() const{return S;}
	auto codesize() const{return C;}
};

////////////////////////////////////////////////////////////////
//                                                            //
//                    class SystematicLDPC                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C>
SystematicLDPC<S,C>::SystematicLDPC(std::uint64_t iterationlimit):
	encoder(H),
	decoder(H),
	iterationlimit(iterationlimit)
{}

template<std::size_t S, std::size_t C>
auto SystematicLDPC<S,C>::encode(const std::bitset<S> &information) const{
	return encoder.systematic_encode(information);
}

template<std::size_t S, std::size_t C>
auto SystematicLDPC<S,C>::encode_redundancy(const std::bitset<S> &information) const{
	return encoder.systematic_redundancy(information);
}

template<std::size_t S, std::size_t C>
template<LDPC::DecoderType T, std::floating_point U>
auto SystematicLDPC<S,C>::decode(const std::array<U,C> &LLR){
	auto QLLR = encoder.inverse_substitution(LLR);
	std::array<U,C> QLPR;//対数事後確率比：列ごとのalphaの和+QLLR

	decoder.decode_init();
	for(auto iter=0ui64; !decoder.iterate<T>(QLPR, QLLR) && iter<iterationlimit; ++iter);

	return decoder.estimate(encoder.substitution(QLPR));
}

}

#endif