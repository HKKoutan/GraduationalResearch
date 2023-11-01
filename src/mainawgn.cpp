#include <thread>
#include <tuple>
#include "common/timekeep.hpp"
#include "common/randombits.hpp"
#include "ldpc/codeLDPC.hpp"
#include "ldpc/AWGN.hpp"

using std::array, std::bitset, std::vector;
using std::size_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 1000ul;
constexpr size_t SOURCE_LENGTH = 252u;
constexpr size_t CODE_LENGTH = 504u;
constexpr size_t NUM_THREADS = 12u;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>1) temp = std::stoull(argv[1]);
	const auto repeat_per_thread = temp;

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<repeat_per_thread<<"*"<<NUM_THREADS<<endl;

	constexpr array noise_factor = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
	//constexpr array noise_factor = {0.0};
	constexpr size_t nsize = noise_factor.size();

	code::Systematic_LDPC<SOURCE_LENGTH,CODE_LENGTH> ldpc;
	array<array<uint64_t,nsize>,3> biterror = {0};
	array<array<uint64_t,nsize>,NUM_THREADS> biterrors = {0};

	//スレッドで動かす関数オブジェクトを定義
	//符号化なし
	auto plain = [&ldpc, repeat_per_thread](size_t t, array<uint64_t,nsize> *bt){
		for(size_t n=0u; n<nsize; ++n){
			channel::AWGN ch(pow(10,-noise_factor[n]*0.1),t);
			util::RandomBits rb(t);
			bitset<SOURCE_LENGTH> info;
			auto &bn = (*bt)[n];

			for(size_t r=0u; r<repeat_per_thread; ++r){
				rb.generate(info);

				auto y = ch.noise<float>(info);
				auto est = ch.estimate(y);

				bn += (est^info).count();
			}
		}
	};
	//符号化あり
	auto encoded = [&ldpc, repeat_per_thread](size_t t, array<uint64_t,nsize> *bt, auto decodertype){
		const auto decoder = ldpc.make_decoder<decltype(decodertype)::type>();
		for(size_t n=0u; n<nsize; ++n){
			channel::AWGN ch(pow(10,-noise_factor[n]*0.1),t);
			util::RandomBits rb(t);
			bitset<SOURCE_LENGTH> info;
			auto &bn = (*bt)[n];

			for(size_t r=0u; r<repeat_per_thread; ++r){
				rb.generate(info);

				auto code = ldpc.encode(info);
				auto y = ch.noise<float>(code);
				auto est = ldpc.decode(ch.LLR(y), decoder);

				bn += (est^info).count();
			}
		}
	};

	//関数オブジェクトをスレッドに与えて実行
	vector<std::thread> threads;
	tk.split();
	for(biterrors = {0}; auto &bt: biterrors) threads.emplace_back(plain, threads.size(), &bt);
	for(auto &t: threads) t.join();
	for(auto &bs = biterror[0]; auto &bt: biterrors) for(size_t n=0u, nend=nsize; n<nend; ++n) bs[n]+=bt[n];
	threads.clear();
	tk.split();
	for(biterrors = {0}; auto &bt: biterrors) threads.emplace_back(encoded, threads.size(), &bt, std::type_identity<code::LDPC::SumProduct_Decoding<SOURCE_LENGTH,CODE_LENGTH>>());
	for(auto &t: threads) t.join();
	for(auto &bs = biterror[1]; auto &bt: biterrors) for(size_t n=0u, nend=nsize; n<nend; ++n) bs[n]+=bt[n];
	threads.clear();
	tk.split();
	for(biterrors = {0}; auto &bt: biterrors) threads.emplace_back(encoded, threads.size(), &bt, std::type_identity<code::LDPC::MinSum_Decoding<SOURCE_LENGTH,CODE_LENGTH>>());
	for(auto &t: threads) t.join();
	for(auto &bs = biterror[2]; auto &bt: biterrors) for(size_t n=0u, nend=nsize; n<nend; ++n) bs[n]+=bt[n];
	tk.stop();

	cout<<"S/N"
	<<"\tPlain"
	<<"\tSum-Product"
	<<"\tMin-sum"
	<<endl;

	for(size_t n=0u; n<nsize; ++n){
		cout<<noise_factor[n]
		<<"\t"<<static_cast<double>(biterror[0][n])/(ldpc.sourcesize()*repeat_per_thread*NUM_THREADS)
		<<"\t"<<static_cast<double>(biterror[1][n])/(ldpc.sourcesize()*repeat_per_thread*NUM_THREADS)
		<<"\t"<<static_cast<double>(biterror[2][n])/(ldpc.sourcesize()*repeat_per_thread*NUM_THREADS)
		<<endl;
	}

	return 0;
}
