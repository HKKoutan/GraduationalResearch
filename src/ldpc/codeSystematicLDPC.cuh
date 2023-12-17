#ifndef INCLUDE_GUARD_ldpc_codeSystematicLDPC
#define INCLUDE_GUARD_ldpc_codeSystematicLDPC

#include <optional>
#include "LDPCdecoding.cuh"
#include "LDPCencoding.cuh"

namespace code {

template<LDPC::CheckMatrix T>//S:Source length, C:Code length
class SystematicLDPC {
	static constexpr std::uint64_t DEFAULT_ITERATION_LIMIT = 100;
	static constexpr std::size_t S = T::sourcesize();
	static constexpr std::size_t C = T::codesize();

	T H;
	LDPC::GenerationMatrix_encoding<decltype(H)> encoder;
	inline static thread_local std::optional<LDPC::Sumproduct_decoding<decltype(H)>> decoders;
	std::uint64_t iterationlimit;//反復回数上限
public:
	explicit SystematicLDPC(const T &H, std::uint64_t iterationlimit);
	explicit SystematicLDPC(const T &H): SystematicLDPC(H, DEFAULT_ITERATION_LIMIT){}

	auto encode(const std::bitset<S> &information) const;//引数から符号語を求める
	auto encode_redundancy(const std::bitset<S> &information) const;//引数から冗長を求める
	template<std::floating_point F, LDPC::boxplusclass P>
	auto decode(const std::array<F,C> &LLR, const P &bp);

	static constexpr auto sourcesize(){return S;}
	static constexpr auto codesize(){return C;}
};

template<std::size_t S, std::size_t C, std::size_t W>//S:Source length, C:Code length
class SystematicLDPC<LDPC::CheckMatrix_regular<S,C,W>> {
	static constexpr std::uint64_t DEFAULT_ITERATION_LIMIT = 100;

	LDPC::CheckMatrix_regular<S,C,W> H;
	LDPC::GenerationMatrix_encoding<decltype(H)> encoder;
	inline static thread_local std::optional<LDPC::Sumproduct_decoding<decltype(H)>> decoders;
	std::uint64_t iterationlimit;//反復回数上限
public:
	explicit SystematicLDPC(const LDPC::CheckMatrix_regular<S,C,W> &H, std::uint64_t iterationlimit);
	explicit SystematicLDPC(const LDPC::CheckMatrix_regular<S,C,W> &H): SystematicLDPC(H, DEFAULT_ITERATION_LIMIT){}

	auto encode(const std::bitset<S> &information) const;//引数から符号語を求める
	auto encode_redundancy(const std::bitset<S> &information) const;//引数から冗長を求める
	template<std::floating_point F, LDPC::boxplusclass P>
	auto decode(const std::array<F,C> &LLR, const P &bp);

	static constexpr auto sourcesize(){return S;}
	static constexpr auto codesize(){return C;}
};

//ヘルパ関数
template<std::size_t S, std::size_t C>
auto make_SystematicLDPC(std::uint64_t iterationlimit){
	return SystematicLDPC(LDPC::validCheckMatrixType_t<S,C>(), iterationlimit);
}

template<std::size_t S, std::size_t C>
auto make_SystematicLDPC(){
	return SystematicLDPC(LDPC::validCheckMatrixType_t<S,C>());
}

////////////////////////////////////////////////////////////////
//                                                            //
//                    class SystematicLDPC                    //
//                                                            //
////////////////////////////////////////////////////////////////

template<LDPC::CheckMatrix T>
SystematicLDPC<T>::SystematicLDPC(const T &H, std::uint64_t iterationlimit):
	H(H),
	encoder(H),
	iterationlimit(iterationlimit)
{}

template<LDPC::CheckMatrix T>
auto SystematicLDPC<T>::encode(const std::bitset<S> &information) const{
	return encoder.systematic_encode(information);
}

template<LDPC::CheckMatrix T>
auto SystematicLDPC<T>::encode_redundancy(const std::bitset<S> &information) const{
	return encoder.systematic_redundancy(information);
}

template<LDPC::CheckMatrix T>
template<std::floating_point F, LDPC::boxplusclass P>
auto SystematicLDPC<T>::decode(const std::array<F,C> &LLR, const P &bp){
	if(!decoders) decoders.emplace(H);

	auto QLLR = encoder.inverse_substitution(LLR);
	std::array<F,C> QLPR;//対数事後確率比：列ごとのalphaの和+QLLR

	decoders->decode_init();
	for(auto iter=0ui64; !decoders->iterate(QLPR, QLLR, bp) && iter<iterationlimit; ++iter);

	return encoder.substitution(QLPR);
}

////////////////////////////////////////////////////////////////
//                                                            //
//         class SystematicLDPC<CheckMatrix_regular>          //
//                                                            //
////////////////////////////////////////////////////////////////

template<std::size_t S, std::size_t C, std::size_t W>
SystematicLDPC<LDPC::CheckMatrix_regular<S,C,W>>::SystematicLDPC(const LDPC::CheckMatrix_regular<S,C,W> &H, std::uint64_t iterationlimit):
	H(H),
	encoder(H),
	iterationlimit(iterationlimit)
{}

template<std::size_t S, std::size_t C, std::size_t W>
auto SystematicLDPC<LDPC::CheckMatrix_regular<S,C,W>>::encode(const std::bitset<S> &information) const{
	return encoder.systematic_encode(information);
}

template<std::size_t S, std::size_t C, std::size_t W>
auto SystematicLDPC<LDPC::CheckMatrix_regular<S,C,W>>::encode_redundancy(const std::bitset<S> &information) const{
	return encoder.systematic_redundancy(information);
}

template<std::size_t S, std::size_t C, std::size_t W>
template<std::floating_point F, LDPC::boxplusclass P>
auto SystematicLDPC<LDPC::CheckMatrix_regular<S,C,W>>::decode(const std::array<F,C> &LLR, const P &bp){
	if(!decoders) decoders.emplace(H);

	auto QLLR = encoder.inverse_substitution(LLR);
	F *QLPR_device = nullptr;//対数事後確率比：列ごとのalphaの和+QLLR
	F *QLLR_device = nullptr;
	cudaMalloc(&QLPR_device, sizeof(F)*C);
	cudaMalloc(&QLLR_device, sizeof(F)*C);
	cudaMemcpy(QLLR_device,QLLR.data(),sizeof(F)*C,cudaMemcpyHostToDevice);

	// decoders->decode_init();
	// for(auto iter=0ui64; !decoders->iterate(QLPR_device, QLLR_device, bp) && iter<iterationlimit; ++iter);
	decoders->decode(QLPR_device, QLLR_device, bp, iterationlimit);

	std::array<F,C> QLPR;
	cudaMemcpy(QLPR.data(),QLPR_device,sizeof(F)*C,cudaMemcpyDeviceToHost);
	cudaFree(QLPR_device);
	cudaFree(QLLR_device);

	return encoder.substitution(QLPR);
}

}

#endif