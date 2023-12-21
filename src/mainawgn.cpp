#include <thread>
#include <tuple>
#include "common/util.hpp"
#include "common/codecommon.hpp"
#include "ldpc/codeSystematicLDPC.hpp"
#include "ldpc/channelAWGN.hpp"

using std::array, std::bitset, std::vector;
using std::size_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 1000;
constexpr size_t SOURCE_LENGTH = 5001;
constexpr size_t CODE_LENGTH = 10000;
constexpr size_t NUM_THREADS = 12;

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
	code::LDPC::phi_table<2> decodertype;

	auto ldpc = code::make_SystematicLDPC<SOURCE_LENGTH,CODE_LENGTH>();
	array<array<uint64_t,nsize>,3> biterror = {0};
	array<array<uint64_t,nsize>,NUM_THREADS> biterrors = {0};

	//スレッドで動かす関数オブジェクトを定義
	//符号化なし
	auto plain = [&ldpc, repeat_per_thread](size_t t, array<uint64_t,nsize> *bt){
		for(size_t n=0; n<nsize; ++n){
			channel::AWGN ch(pow(10,-noise_factor[n]*0.1),t);
			util::RandomBits rb(t);
			bitset<SOURCE_LENGTH> info;
			auto &bn = (*bt)[n];

			for(size_t r=0u; r<repeat_per_thread; ++r){
				rb.generate(info);

				auto y = ch.noise<float>(info);
				auto est = code::estimate_crop(y);

				bn += (est^info).count();
			}
		}
	};
	//符号化あり
	auto encoded = [&ldpc, repeat_per_thread](size_t t, array<uint64_t,nsize> *bt, auto decodertype){
		for(size_t n=0; n<nsize; ++n){
			channel::AWGN ch(pow(10,-noise_factor[n]*0.1),t);
			util::RandomBits rb(t);
			bitset<SOURCE_LENGTH> info;
			auto &bn = (*bt)[n];

			for(size_t r=0u; r<repeat_per_thread; ++r){
				rb.generate(info);

				auto code = ldpc.encode(info);
				auto y = ch.noise<float>(code);
				auto LLR = ch.LLR(y);
				LLR = ldpc.decode(LLR, decodertype);
				auto est = code::estimate_crop<SOURCE_LENGTH>(LLR);

				bn += (est^info).count();
			}
		}
	};

	//関数オブジェクトをスレッドに与えて実行
	vector<std::thread> threads;
	tk.split();
	for(biterrors = {0}; auto &bt: biterrors) threads.emplace_back(plain, threads.size(), &bt);
	for(auto &t: threads) t.join();
	for(auto &bs = biterror[0]; auto &bt: biterrors) for(size_t n=0; n<nsize; ++n) bs[n]+=bt[n];
	threads.clear();
	tk.split();
	for(biterrors = {0}; auto &bt: biterrors) threads.emplace_back(encoded, threads.size(), &bt, decodertype);
	for(auto &t: threads) t.join();
	for(auto &bs = biterror[1]; auto &bt: biterrors) for(size_t n=0; n<nsize; ++n) bs[n]+=bt[n];
	// threads.clear();
	// tk.split();
	// for(biterrors = {0}; auto &bt: biterrors) threads.emplace_back(encoded, threads.size(), &bt, decltype(ldpc)::DecoderType::MinSum());
	// for(auto &t: threads) t.join();
	// for(auto &bs = biterror[2]; auto &bt: biterrors) for(size_t n=0; n<nsize; ++n) bs[n]+=bt[n];
	tk.stop();

	cout<<"S/N"
	<<"\tPlain"
	<<"\tSum-Product"
	// <<"\tMin-sum"
	<<endl;

	for(size_t n=0; n<nsize; ++n){
		cout<<noise_factor[n]
		<<"\t"<<static_cast<double>(biterror[0][n])/(ldpc.sourcesize()*repeat_per_thread*NUM_THREADS)
		<<"\t"<<static_cast<double>(biterror[1][n])/(ldpc.sourcesize()*repeat_per_thread*NUM_THREADS)
		// <<"\t"<<static_cast<double>(biterror[2][n])/(ldpc.sourcesize()*repeat_per_thread*NUM_THREADS)
		<<endl;
	}

	return 0;
}
