﻿#ifndef __code_CODELDPC__
#define __code_CODELDPC__

#include "LDPCbase.hpp"

namespace code {

template<std::size_t S, std::size_t C>//S:Source length, C:Code length
class Systematic_LDPC {
	static constexpr std::uint64_t DEFAULT_ITERATION_LIMIT = 100;

	std::shared_ptr<const LDPC::CheckMatrix<S,C>> H;
	std::unique_ptr<const LDPC::I_LDPC_Encoding<S,C>> encoder;
	std::uint64_t iterationlimit;//反復回数上限
public:
	explicit Systematic_LDPC(std::uint64_t iterationlimit);
	Systematic_LDPC(): Systematic_LDPC(DEFAULT_ITERATION_LIMIT){}

	auto encode(const std::bitset<S> &information) const;//引数から冗長を求める
	auto decode(const std::array<LDPC::fptype,C> &LLR, const std::unique_ptr<LDPC::I_LDPC_Decoding<S,C>> &decoder);

	template<typename T> const auto make_decoder(){return std::unique_ptr<LDPC::I_LDPC_Decoding<S,C>>(new T(H));}

	auto sourcesize() const{return H->sourcesize();}
	auto codesize() const{return H->codesize();}
};

////////////////////////////////////////////////////////////////
//                                                            //
//                      Systematic_LDPC                       //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C>
Systematic_LDPC<S,C>::Systematic_LDPC(std::uint64_t iterationlimit):
	H(std::make_shared<const LDPC::CheckMatrix<S,C>>(LDPC::CheckMatrix<S,C>())),
	encoder(new LDPC::Generation_Matrix_Encoding<S,C>(*H)),
	iterationlimit(iterationlimit)
{}

template<std::size_t S, std::size_t C>
auto Systematic_LDPC<S,C>::encode(const std::bitset<S> &information) const{
	std::bitset<C> code;
	auto parity = encoder->systematic_parity(information);
	for(std::size_t i=0u, iend=S; i<iend; ++i) code[i]=information[i];
	for(std::size_t i=S, iend=C; i<iend; ++i) code[i]=parity[i-S];
	return code;
}

template<std::size_t S, std::size_t C>
auto Systematic_LDPC<S,C>::decode(const std::array<LDPC::fptype,C> &LLR, const std::unique_ptr<LDPC::I_LDPC_Decoding<S,C>> &decoder){
	auto QLLR = encoder->inverse_substitution(LLR);
	std::array<LDPC::fptype,C> QLPR;//対数事後確率比：列ごとのalphaの和+QLLR

	decoder->decode_init();
	for(auto iter=0ui64; !decoder->iterate(QLPR, QLLR) && iter<iterationlimit; ++iter);

	return decoder->estimate(encoder->substitution(QLPR));
}

}

#endif