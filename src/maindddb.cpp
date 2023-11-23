#include <array>
#include <tuple>
#include <thread>
#include <numeric>
#include "dnas/codeDNASadapter.hpp"
#include "dnas/codeDNASbalancing.hpp"
#include "dnas/codeDNASstats.hpp"
#include "dnas/channelsequencer.hpp"
#include "ldpc/codeSystematicLDPC.hpp"
#include "common/util.hpp"
#include "common/codecommon.hpp"

using std::array, std::bitset, std::vector, std::tuple, std::pair;
using std::size_t, std::uint8_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;
using code::DNAS::nucleotide_t;

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 1000;
constexpr size_t SOURCE_LENGTH = 512;
constexpr size_t CODE_LENGTH = 1024;
constexpr size_t NUM_THREADS = 12;
constexpr size_t BLOCK_SIZE = 0;
constexpr uint8_t ATGC = 0x1B;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	cout<<"Title: DNA storage simulation on nanopore sequencing channel with differential encoding and division balancing."<<endl;
	cout<<"ATGC: "<<std::bitset<8>(ATGC)<<endl;

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>1) temp = std::stoull(argv[1]);
	const auto repeat_per_thread = temp;

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<repeat_per_thread<<"*"<<NUM_THREADS<<endl;

	// constexpr array noise_factor = {0.0};
	constexpr array noise_factor = {0.04,0.03,0.02,0.01,0.0};
	// constexpr array noise_factor = {0.04,0.035,0.03,0.025,0.02,0.015,0.01,0.005,0.0};
	constexpr size_t nsize = noise_factor.size();

	auto ldpc = code::make_SystematicLDPC<SOURCE_LENGTH,CODE_LENGTH>();
	// tuple: biterrors, nterrors, GCper(average,var)
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,pair<double,double>>,5> stat = {};
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>>,NUM_THREADS> stats = {};

	auto plain = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcper = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto cm = code::DNAS::differential<ATGC>::encode(m);

				auto cmbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,0>::balance(cm);

				auto qty_AT = code::DNAS::countAT(cmbar);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cmbar.size());

				auto rm = ch.noise(cmbar);
				// auto rm=cmbar;

				auto mest = code::DNAS::differential<ATGC>::decode(rm);

				nterror[n] += code::DNAS::countDifferentialError(cm,rm);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto encoded = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcper = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto c = ldpc.encode(m);
				auto cc = code::DNAS::differential<ATGC>::encode(c);

				auto ccbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,0>::balance(cc);

				auto qty_AT = code::DNAS::countAT(ccbar);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(ccbar.size());

				auto rc = ch.noise(ccbar);
				// auto rc=ccbar;

				auto Lrc = ch.likelihood<float>(rc);
				auto LLR = code::DNAS::differential<ATGC>::decode_p(Lrc);

				auto LLRest = ldpc.decode<decltype(ldpc)::DecoderType::SumProduct>(LLR);
				auto mest = code::estimate_crop<SOURCE_LENGTH>(LLRest);

				nterror[n] += code::DNAS::countDifferentialError(cc,rc);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto plain_minchange = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcper = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto cm = code::DNAS::differential<ATGC>::encode(m);

				auto cmbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,2>::balance(cm);

				auto qty_AT = code::DNAS::countAT(cmbar);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cmbar.size());

				auto rm = ch.noise(cmbar);
				// auto rm=cmbar;

				auto mest = code::DNAS::differential<ATGC>::decode(rm);

				nterror[n] += code::DNAS::countDifferentialError(cm,rm);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto encoded_minchange = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcper = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto c = ldpc.encode(m);
				auto cc = code::DNAS::differential<ATGC>::encode(c);

				auto ccbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,2>::balance(cc);

				auto qty_AT = code::DNAS::countAT(ccbar);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(ccbar.size());

				auto rc = ch.noise(ccbar);
				// auto rc=ccbar;

				auto Lrc = ch.likelihood<float>(rc);
				auto LLR = code::DNAS::differential<ATGC>::decode_p(Lrc);

				auto LLRest = ldpc.decode<decltype(ldpc)::DecoderType::SumProduct>(LLR);
				auto mest = code::estimate_crop<SOURCE_LENGTH>(LLRest);

				nterror[n] += code::DNAS::countDifferentialError(cc,rc);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto aggregate = [&stat, &stats, repeat_per_thread](std::size_t dest){
		for(auto &st: stats) for(size_t n=0, nend=nsize; n<nend; ++n){
			std::get<0>(stat[dest])[n] += std::get<0>(st)[n];
			std::get<1>(stat[dest])[n] += std::get<1>(st)[n];
		}
		auto sum = 0.0, sqsum = 0.0;
		for(auto &st: stats) for(auto &pn: std::get<2>(st)) for(auto &pr: pn){
			sum += pr;
			sqsum += pr*pr;
		}
		auto num = static_cast<double>(nsize*NUM_THREADS*repeat_per_thread);
		auto average = sum/num;
		std::get<2>(stat[dest]).first = average;
		std::get<2>(stat[dest]).second = sqsum/num - average*average;
	};

	auto result = [&stat, repeat_per_thread](std::size_t target, std::size_t channel_size){
		cout<<"GCper var: "<<std::get<2>(stat[target]).second<<", ave: "<<std::get<2>(stat[target]).first<<endl;
		cout<<"Noise factor"
		<<"\tBER"
		<<"\tNER"
		<<endl;
		for(size_t n=0; n<noise_factor.size(); n++){
			cout<<noise_factor[n]
			<<"\t"<<static_cast<double>(std::get<0>(stat[target])[n])/static_cast<double>(SOURCE_LENGTH*NUM_THREADS*repeat_per_thread)
			<<"\t"<<static_cast<double>(std::get<1>(stat[target])[n])/static_cast<double>(channel_size/2*NUM_THREADS*repeat_per_thread)
			<<endl;
		}
	};

	//実行
	vector<std::thread> threads;
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(0);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(1);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain_minchange, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(2);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_minchange, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(3);
	tk.stop();

	//結果表示
	cout<<"Block Size: "<<BLOCK_SIZE<<endl;
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain"<<endl;
	result(0, SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded"<<endl;
	result(1, CODE_LENGTH);
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain(lesschange)"<<endl;
	result(2, SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded(lesschange)"<<endl;
	result(3, CODE_LENGTH);

	return 0;
}
