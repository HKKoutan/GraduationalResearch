#include <array>
#include <tuple>
#include <thread>
#include <numeric>
#include "dnas/DNASnttype.hpp"
#include "dnas/codeDNAS.hpp"
#include "dnas/sequencer.hpp"
#include "ldpc/codeLDPC.hpp"
#include "common/randombits.hpp"
#include "common/timekeep.hpp"

using std::array, std::bitset, std::vector, std::tuple, std::pair;
using std::size_t, std::uint8_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;
using code::DNAS::nucleotide_t;

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 1000u;
constexpr size_t SOURCE_LENGTH = 512u;
constexpr size_t CODE_LENGTH = 1024u;
constexpr size_t NUM_THREADS = 12u;
constexpr size_t BLOCK_SIZE = 0u;
constexpr uint8_t AGTC = 0x27;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	cout<<"Title: DNA storage simulation on nanopore sequencing channel with differential encoding and division balancing."<<endl;

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>1) temp = std::stoull(argv[1]);
	const auto repeat_per_thread = temp;

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<repeat_per_thread<<"*"<<NUM_THREADS<<endl;

	// constexpr array noise_factor = {0.0};
	constexpr array noise_factor = {0.04,0.03,0.02,0.01,0.0};
	// constexpr array noise_factor = {0.04,0.035,0.03,0.025,0.02};
	constexpr size_t nsize = noise_factor.size();

	code::Systematic_LDPC<SOURCE_LENGTH,CODE_LENGTH> ldpc;
	// tuple: biterrors, nterrors, GCper(average,var)
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,pair<double,double>>,4> stat = {};
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>>,NUM_THREADS> stats = {};

	auto plain_ATGC = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcper = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::Nanopore_Sequencing ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto qm = code::DNAS::binary_to_quarternary(m);
				auto cm = code::DNAS::differential::encode(qm);
				// for(auto &ci: cm) ci = 3;

				auto qty_AT = code::DNAS::count_AT(cm);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cm.size());

				auto rm = ch.noise(cm);
				// auto rm=cm;

				auto qmest = code::DNAS::differential::decode(rm);
				auto mest = code::DNAS::quarternary_to_binary(qmest);
				{
					uint64_t acc = 0u;
					for(size_t i=0u, iend=cm.size(); i<iend; ++i) acc += (cm[i]!=rm[i]);
					nterror[n] += acc;
				}
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto encoded_ATGC = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcper = std::get<2>(*st);
		auto decoder = ldpc.make_decoder<code::LDPC::SumProduct_Decoding<SOURCE_LENGTH,CODE_LENGTH>>();
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::Nanopore_Sequencing ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto c = ldpc.encode(m);
				auto qc = code::DNAS::binary_to_quarternary(c);
				auto cc = code::DNAS::differential::encode(qc);

				auto qty_AT = code::DNAS::count_AT(cc);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cc.size());

				auto rc = ch.noise(cc);
				// auto rc=cc;

				auto LLR = ch.differential_LLR(rc);

				auto mest = ldpc.decode(LLR, decoder);
				{
					uint64_t acc = 0u;
					for(size_t i=0u, iend=cc.size(); i<iend; ++i) acc += (cc[i]!=rc[i]);
					nterror[n] += acc;
				}
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto plain_AGTC = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcper = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::Nanopore_Sequencing<AGTC> ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto qm = code::DNAS::binary_to_quarternary<AGTC>(m);
				auto cm = code::DNAS::differential::encode(qm);
				// for(auto &ci: cm) ci = 3;

				auto qty_AT = code::DNAS::count_AT(cm);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cm.size());

				auto rm = ch.noise(cm);
				// auto rm=cm;

				auto qmest = code::DNAS::differential::decode(rm);
				auto mest = code::DNAS::quarternary_to_binary(qmest);
				{
					uint64_t acc = 0u;
					for(size_t i=0u, iend=cm.size(); i<iend; ++i) acc += (cm[i]!=rm[i]);
					nterror[n] += acc;
				}
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto encoded_AGTC = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcper = std::get<2>(*st);
		auto decoder = ldpc.make_decoder<code::LDPC::SumProduct_Decoding<SOURCE_LENGTH,CODE_LENGTH>>();
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::Nanopore_Sequencing<AGTC> ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto c = ldpc.encode(m);
				auto qc = code::DNAS::binary_to_quarternary<AGTC>(c);
				auto cc = code::DNAS::differential::encode(qc);

				auto qty_AT = code::DNAS::count_AT(cc);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cc.size());

				auto rc = ch.noise(cc);
				// auto rc=cc;

				auto LLR = ch.differential_LLR(rc);

				auto mest = ldpc.decode(LLR, decoder);
				{
					uint64_t acc = 0u;
					for(size_t i=0u, iend=cc.size(); i<iend; ++i) acc += (cc[i]!=rc[i]);
					nterror[n] += acc;
				}
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
	for(stats = {}; auto &st: stats) threads.emplace_back(plain_ATGC, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(0);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_ATGC, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(1);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain_AGTC, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(2);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_AGTC, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(3);
	tk.stop();

	//結果表示
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain ATGC"<<endl;
	result(0, SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded ATGC(no balancing)"<<endl;
	result(1, CODE_LENGTH);
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain AGTC"<<endl;
	result(2, SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded AGTC(no balancing)"<<endl;
	result(3, CODE_LENGTH);

	return 0;
}
