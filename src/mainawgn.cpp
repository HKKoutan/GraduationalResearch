#include <thread>
#include "common/timekeep.hpp"
#include "common/randombits.hpp"
#include "ldpc/codeLDPC.hpp"
#include "ldpc/AWGN.hpp"

using std::array, std::bitset, std::vector;
using std::size_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 1000ul;
constexpr size_t SOURCE_LENGTH = 256u;
constexpr size_t CODE_LENGTH = 512u;
constexpr size_t NUM_THREADS = 12u;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	// if(argc<2){
	// 	cerr<<"Usage: "<<argv[0]<<" FILENAME"<<endl;
	// 	return 1;
	// }

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>1) temp = std::stoull(argv[1]);
	const auto repeat_per_thread = temp;

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<repeat_per_thread<<"*"<<NUM_THREADS<<endl;

	constexpr array noise_factor = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
	//constexpr array noise_factor = {0.0};

	vector<std::thread> threads;
	array<array<uint64_t, noise_factor.size()>, NUM_THREADS> pl_biterrors = {0};
	array<array<uint64_t, noise_factor.size()>, NUM_THREADS> sp_biterrors = {0};
	array<array<uint64_t, noise_factor.size()>, NUM_THREADS> ms_biterrors = {0};

	code::Systematic_LDPC<SOURCE_LENGTH,CODE_LENGTH> ldpc;

	tk.split();
	for(auto &bt: pl_biterrors){
		auto t = threads.size();
		threads.emplace_back([&ldpc, &bt, &noise_factor, repeat_per_thread, t](){
			for(size_t n=0u; n<noise_factor.size(); ++n){
				bitset<SOURCE_LENGTH> info;
				auto &bn = bt[n];

				channel::AWGN<float,SOURCE_LENGTH> ch(pow(10,-noise_factor[n]*0.1),t);
				util::RandomBits<SOURCE_LENGTH> rb(t);

				for(size_t r=0u; r<repeat_per_thread; ++r){
					rb.generate(info);

					auto y = ch.noise(info);

					for(size_t i=0u, iend=info.size(); i<iend; ++i) bn += info[i]^(y[i]<0.0);
				}
			}
		});
	}
	for(auto &t: threads) t.join();

	tk.split();
	threads.clear();
	for(auto &bt: sp_biterrors){
		auto t = threads.size();
		threads.emplace_back([&ldpc, &bt, &noise_factor, repeat_per_thread, t](){
			const auto decoder = ldpc.create_decoder<code::LDPC::SumProduct_Decoding<SOURCE_LENGTH,CODE_LENGTH>>();

			for(size_t n=0u; n<noise_factor.size(); ++n){
				bitset<SOURCE_LENGTH> info;
				auto &bn = bt[n];

				channel::AWGN<float,CODE_LENGTH> ch(pow(10,-noise_factor[n]*0.1),t);
				util::RandomBits<SOURCE_LENGTH> rb(t);

				for(size_t r=0u; r<repeat_per_thread; ++r){
					rb.generate(info);

					auto code = ldpc.encode(info);

					auto y = ch.noise(code);

					auto LLR = ch.LLR(y);
					auto est = ldpc.decode(LLR, decoder);

					bn += (est^info).count();
				}
			}
		});
	}
	for(auto &t: threads) t.join();

	tk.split();
	threads.clear();
	for(auto &bt: ms_biterrors){
		auto t = threads.size();
		threads.emplace_back([&ldpc, &bt, &noise_factor, repeat_per_thread, t](){
			const auto decoder = ldpc.create_decoder<code::LDPC::MinSum_Decoding<SOURCE_LENGTH,CODE_LENGTH>>();

			for(size_t n=0u; n<noise_factor.size(); ++n){
				bitset<SOURCE_LENGTH> info;
				auto &bn = bt[n];

				channel::AWGN<float,CODE_LENGTH> ch(pow(10,-noise_factor[n]*0.1),t);
				util::RandomBits<SOURCE_LENGTH> rb(t);

				for(size_t r=0u; r<repeat_per_thread; ++r){
					rb.generate(info);

					auto code = ldpc.encode(info);

					auto y = ch.noise(code);

					auto LLR = ch.LLR(y);
					auto est = ldpc.decode(LLR, decoder);

					bn += (est^info).count();
				}
			}
		});
	}
	for(auto &t: threads) t.join();

	tk.split();

	array<uint64_t, noise_factor.size()> pl_biterror = {0};
	array<uint64_t, noise_factor.size()> sp_biterror = {0};
	array<uint64_t, noise_factor.size()> ms_biterror = {0};

	for(size_t t=0u; t<NUM_THREADS; ++t){
		auto &pt = pl_biterrors[t], &st = sp_biterrors[t], &mt = ms_biterrors[t];
		for(size_t n=0u; n<noise_factor.size(); ++n){
			pl_biterror[n] += pt[n];
			sp_biterror[n] += st[n];
			ms_biterror[n] += mt[n];
		}
	}

	tk.stop();

	cout<<"S/N"
	<<"\tPlain"
	<<"\tSum-Product"
	<<"\tMin-sum"
	<<endl;

	for(size_t n=0u; n<noise_factor.size(); ++n){
		cout<<noise_factor[n]
		<<"\t"<<static_cast<double>(pl_biterror[n])/(ldpc.sourcesize()*repeat_per_thread*NUM_THREADS)
		<<"\t"<<static_cast<double>(sp_biterror[n])/(ldpc.sourcesize()*repeat_per_thread*NUM_THREADS)
		<<"\t"<<static_cast<double>(ms_biterror[n])/(ldpc.sourcesize()*repeat_per_thread*NUM_THREADS)
		<<endl;
	}

	return 0;
}
