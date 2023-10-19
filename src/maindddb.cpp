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
using std::size_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;
using code::DNAS::nucleotide_t;

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 1000u;
constexpr size_t SOURCE_LENGTH = 256u;
constexpr size_t CODE_LENGTH = 512u;
constexpr size_t NUM_THREADS = 12u;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	cout<<"Title: DNA storage simulation on nanopore sequencing channel with differential encoding and division balancing."<<endl;

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>1) temp = std::stoull(argv[1]);
	const auto repeat_per_thread = temp;

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<repeat_per_thread<<"*"<<NUM_THREADS<<endl;

	// constexpr array noise_factor = {0};
	constexpr array noise_factor = {0.04,0.035,0.03,0.025,0.02};
	constexpr size_t nsize = noise_factor.size();

	// code::Systematic_LDPC<SOURCE_LENGTH,CODE_LENGTH> ldpc;
	// tuple: biterrors, bitcounts, nterrors, ntcounts, GCper(average,var)
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,pair<double,double>>,2> stat = {};
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>>,NUM_THREADS> stats = {};

	auto plain = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
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

	// auto encoded = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>> *st){
	// 	auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st), &ntcount = std::get<3>(*st);
	// 	auto &gcper = std::get<4>(*st);
	// 	auto decoder = ldpc.make_decoder<code::LDPC::SumProduct_Decoding<SOURCE_LENGTH,CODE_LENGTH>>();
	// 	for(size_t n=0; n<nsize; ++n){
	// 		bitset<SOURCE_LENGTH> m;
	// 		channel::Nanopore_Sequencing ch(noise_factor[n],t);
	// 		util::RandomBits rb(t);
	// 		gcper[n].resize(repeat_per_thread);

	// 		for(size_t r=0u; r<repeat_per_thread; r++){
	// 			rb.generate(m);

	// 			auto [cm, mmask] = code::DNAS::VLRLL_encode(m);
	// 			auto tm = code::DNAS::interim_map(cm);
	// 			auto tr = ldpc.encode_redundancy(tm);
	// 			auto cr = code::DNAS::modified_VLRLL_encode(tr, cm[cm.size()-1]);

	// 			auto qty_AT = code::DNAS::nt_qty_count(cm, code::DNAS::nt_qty_count(cr));
	// 			gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cm.size()+cr.size());

	// 			auto rm = ch.noise(cm);
	// 			auto rr = ch.noise(cr);
	// 			// auto rm=cm;
	// 			// auto rr=cr;
	// 			// auto test = tm;

	// 			auto LLR = ch.VLRLL_LLR(rm, rr, rm[rm.size()-1]);

	// 			auto test = ldpc.decode(LLR, decoder);
	// 			auto cest = code::DNAS::interim_demap(test);
	// 			auto mest = code::DNAS::VLRLL_decode(cest);
	// 			{
	// 				uint64_t acc = 0u;
	// 				for(size_t i=0u, iend=cm.size(); i<iend; ++i) acc += (cm[i]!=cest[i]);
	// 				nterror[n] = acc;
	// 			}
	// 			ntcount[n] += cm.size();
	// 			biterror[n] += (mest^(m&mmask)).count();
	// 			bitcount[n] += mmask.count();
	// 		}
	// 	}
	// };

	vector<std::thread> threads;
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain, threads.size(), &st);
	for(auto &t: threads) t.join();
	for(auto &st: stats){
		for(size_t n=0, nend=nsize; n<nend; ++n){
			std::get<0>(stat[0])[n] += std::get<0>(st)[n];
			std::get<1>(stat[0])[n] += std::get<1>(st)[n];
		}
	}
	{
		auto sum = 0.0, sqsum = 0.0;
		for(auto &st: stats){
			for(auto &pn: std::get<2>(st)){
				for(auto &pr: pn){
					sum += pr;
					sqsum += pr*pr;
				}
			}
		}
		auto num = static_cast<double>(nsize*NUM_THREADS*repeat_per_thread);
		auto average = sum/num;
		std::get<2>(stat[0]).first = average;
		std::get<2>(stat[0]).second = sqsum/num - average*average;
	}
	// threads.clear();
	// tk.split();
	// for(stats = {}; auto &st: stats) threads.emplace_back(encoded, threads.size(), &st);
	// for(auto &t: threads) t.join();
	// for(auto &st: stats){
	// 	for(size_t n=0, nend=nsize; n<nend; ++n){
	// 		std::get<0>(stat[1])[n] += std::get<0>(st)[n];
	// 		std::get<1>(stat[1])[n] += std::get<1>(st)[n];
	// 		std::get<2>(stat[1])[n] += std::get<2>(st)[n];
	// 		std::get<3>(stat[1])[n] += std::get<3>(st)[n];
	// 	}
	// }
	// {
	// 	auto sum = 0.0, sqsum = 0.0;
	// 	for(auto &st: stats){
	// 		for(auto &pn: std::get<4>(st)){
	// 			for(auto &pr: pn){
	// 				sum += pr;
	// 				sqsum += pr*pr;
	// 			}
	// 		}
	// 	}
	// 	auto num = static_cast<double>(nsize*NUM_THREADS*repeat_per_thread);
	// 	auto average = sum/num;
	// 	std::get<4>(stat[1]).first = average;
	// 	std::get<4>(stat[1]).second = sqsum/num - average*average;
	// }
	tk.stop();

	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain"<<endl;
	cout<<"GCper var: "<<std::get<2>(stat[0]).second<<", ave: "<<std::get<2>(stat[0]).first<<endl;
	cout<<"Noise factor"
	<<"\tBER"
	<<"\tNER"
	<<endl;

	for(size_t n=0; n<noise_factor.size(); n++){
		cout<<noise_factor[n]
		<<"\t"<<static_cast<double>(std::get<0>(stat[0])[n])/static_cast<double>(SOURCE_LENGTH*NUM_THREADS*repeat_per_thread)
		<<"\t"<<static_cast<double>(std::get<1>(stat[0])[n])/static_cast<double>(SOURCE_LENGTH/2*NUM_THREADS*repeat_per_thread)
		<<endl;
	}
	// cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	// cout<<"encoded(no balancing)"<<endl;
	// cout<<"GCper var: "<<std::get<4>(stat[1]).second<<", ave: "<<std::get<4>(stat[1]).first<<endl;
	// cout<<"Noise factor"
	// <<"\tBER"
	// <<"\tNER"
	// <<endl;

	// for(size_t n=0; n<noise_factor.size(); n++){
	// 	cout<<noise_factor[n]
	// 	<<"\t"<<static_cast<double>(std::get<0>(stat[1])[n])/static_cast<double>(SOURCE_LENGTH*NUM_THREADS*repeat_per_thread)
	// 	<<"\t"<<static_cast<double>(std::get<2>(stat[1])[n])/static_cast<double>(std::get<3>(stat[1])[n])
	// 	<<endl;
	// }
	return 0;
}
