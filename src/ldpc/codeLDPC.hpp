#ifndef INCLUDE_GUARD_ldpc_codeLDPC
#define INCLUDE_GUARD_ldpc_codeLDPC

#include "LDPCbase.hpp"

namespace code {

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class Systematic_LDPC {
	static constexpr std::uint64_t DEFAULT_ITERATION_LIMIT = 100;
	std::unique_ptr<const LDPC::I_LDPC_Encoding<LDPC::validCheckMatrixType_t<S,C>>> encoder;
	std::uint64_t iterationlimit;//反復回数上限
public:
	struct DecoderType{
		struct SumProduct {using type = LDPC::SumProduct_Decoding<LDPC::validCheckMatrixType_t<S,C>>;};
		struct MinSum {using type = LDPC::MinSum_Decoding<LDPC::validCheckMatrixType_t<S,C>>;};
	};

	explicit Systematic_LDPC(std::uint64_t iterationlimit);
	Systematic_LDPC(): Systematic_LDPC(DEFAULT_ITERATION_LIMIT){}

	auto encode(const std::bitset<S> &information) const;//引数から符号語を求める
	auto encode_redundancy(const std::bitset<S> &information) const;//引数から冗長を求める
	auto decode(const std::array<LDPC::fptype,C> &LLR, const std::unique_ptr<LDPC::I_LDPC_Decoding<LDPC::validCheckMatrixType_t<S,C>>> &decoder);

	template<std::derived_from<LDPC::I_LDPC_Decoding<LDPC::validCheckMatrixType_t<S,C>>> T>
	auto make_decoder(){return std::unique_ptr<LDPC::I_LDPC_Decoding<LDPC::validCheckMatrixType_t<S,C>>>(new T(LDPC::validCheckMatrixType_t<S,C>()));}

	auto sourcesize() const{return S;}
	auto codesize() const{return C;}
};

////////////////////////////////////////////////////////////////
//                                                            //
//                      Systematic_LDPC                       //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C>
Systematic_LDPC<S,C>::Systematic_LDPC(std::uint64_t iterationlimit):
	encoder(new LDPC::Generation_Matrix_Encoding(LDPC::validCheckMatrixType_t<S,C>())),
	iterationlimit(iterationlimit)
{}

template<std::size_t S, std::size_t C>
auto Systematic_LDPC<S,C>::encode(const std::bitset<S> &information) const{
	std::bitset<C> code;
	auto parity = encoder->systematic_redundancy(information);
	for(std::size_t i=0u, iend=S; i<iend; ++i) code[i]=information[i];
	for(std::size_t i=S, iend=C; i<iend; ++i) code[i]=parity[i-S];
	return code;
}

template<std::size_t S, std::size_t C>
auto Systematic_LDPC<S,C>::encode_redundancy(const std::bitset<S> &information) const{
	return encoder->systematic_redundancy(information);
}

template<std::size_t S, std::size_t C>
auto Systematic_LDPC<S,C>::decode(const std::array<LDPC::fptype,C> &LLR, const std::unique_ptr<LDPC::I_LDPC_Decoding<LDPC::validCheckMatrixType_t<S,C>>> &decoder){
	auto QLLR = encoder->inverse_substitution(LLR);
	std::array<LDPC::fptype,C> QLPR;//対数事後確率比：列ごとのalphaの和+QLLR

	decoder->decode_init();
	for(auto iter=0ui64; !decoder->iterate(QLPR, QLLR) && iter<iterationlimit; ++iter);

	return decoder->estimate(encoder->substitution(QLPR));
}

}

#endif