#include <array>
#include <tuple>
#include <thread>
#include <numeric>
#include "dnas/DNASnttype.hpp"
#include "dnas/codeDNASadapter.hpp"
#include "dnas/codeDNASRLL.hpp"
#include "dnas/codeDNASstats.hpp"
#include "dnas/channelsequencer.hpp"
#include "ldpc/codeSystematicLDPC.hpp"
#include "common/util.hpp"
#include "common/codecommon.hpp"

using std::array, std::bitset, std::vector, std::tuple, std::pair;
using std::size_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;
using code::DNAS::nucleotide_t;

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 1000;
constexpr size_t SOURCE_LENGTH = 512;
constexpr size_t CODE_LENGTH = 1024;
constexpr size_t NUM_THREADS = 12;

constexpr std::uint8_t ATGC = 0x1B;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	cout<<"Title: DNA storage simulation on nanopore sequencing channel with VLRLL encoding and simple conversion or differential encoding."<<endl;
	cout<<"ATGC: "<<std::bitset<8>(ATGC)<<endl;

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>1) temp = std::stoull(argv[1]);
	const auto repeat_per_thread = temp;

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<repeat_per_thread<<"*"<<NUM_THREADS<<endl;

	// constexpr array noise_factor = {0};
	constexpr array noise_factor = {0.04,0.035,0.03,0.025,0.02};
	// constexpr array noise_factor = {0.04,0.035,0.03,0.025,0.02,0.015,0.01,0.005,0.0};
	constexpr size_t nsize = noise_factor.size();

	auto ldpc = code::make_SystematicLDPC<SOURCE_LENGTH,CODE_LENGTH>();
	// tuple: biterrors, bitcounts, nterrors, GCper(average,var), maxrunlength
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,pair<double,double>,std::size_t>,3> stat = {};
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>,std::size_t>,NUM_THREADS> stats = {};

	auto plain = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>,std::size_t> *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &gcper = std::get<3>(*st);
		auto &maxrunlength = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);

				auto qty_AT = code::DNAS::countAT(cm);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cm.size());
				auto runlength = code::DNAS::countRunlength(cm);
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cm);
				// auto rm=cm;

				auto mest = code::DNAS::VLRLL<ATGC>::decode(rm);

				nterror[n] += code::DNAS::countDifferentialError(cm,rm);
				biterror[n] += ((mest&mmask)^(m&mmask)).count();
				bitcount[n] += mmask.count();
			}
		}
	};

	auto encoded_conv = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>,std::size_t> *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &gcper = std::get<3>(*st);
		auto &maxrunlength = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);
				auto tm = code::DNAS::convert<ATGC>::nttype_to_binary(cm);
				auto tr = ldpc.encode_redundancy(tm);
				auto cr = code::DNAS::modifiedVLRLL<ATGC>::encode(tr, cm.back(), run);

				auto qty_AT = code::DNAS::countAT(cm) + code::DNAS::countAT(cr);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cm.size()+cr.size());
				auto runlength = code::DNAS::countRunlength(code::concatenate(cm,cr));
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cm);
				auto rr = ch.noise(cr);
				// auto rm=cm;
				// auto rr=cr;

				auto Lrm = ch.likelihood<float>(rm);
				auto Lrr = ch.likelihood<float>(rr);
				auto LLR = code::concatenate(code::DNAS::convert<ATGC>::nttype_to_binary_p(Lrm), code::DNAS::modifiedVLRLL<ATGC>::decode_p(Lrr, Lrm.back()));

				auto LLRest = ldpc.decode<decltype(ldpc)::DecoderType::SumProduct>(LLR);
				// auto LLRest = LLR;
				auto test = code::estimate_crop<SOURCE_LENGTH>(LLRest);
				auto cest = code::DNAS::convert<ATGC>::binary_to_nttype(test);
				auto mest = code::DNAS::VLRLL<ATGC>::decode(cest);

				nterror[n] += code::DNAS::countError(cm,cest);
				biterror[n] += ((mest&mmask)^(m&mmask)).count();
				bitcount[n] += mmask.count();
			}
		}
	};

	auto encoded_diff = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,array<vector<double>,nsize>,std::size_t> *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &gcper = std::get<3>(*st);
		auto &maxrunlength = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcper[n].resize(repeat_per_thread);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);
				auto tm = code::DNAS::differential<ATGC>::decode(cm);
				auto tr = ldpc.encode_redundancy(tm);
				auto cr = code::DNAS::modifiedVLRLL<ATGC>::encode(tr, cm.back(), run);

				auto qty_AT = code::DNAS::countAT(cm) + code::DNAS::countAT(cr);
				gcper[n][r] = 1.0 - static_cast<double>(qty_AT)/static_cast<double>(cm.size()+cr.size());
				auto runlength = code::DNAS::countRunlength(code::concatenate(cm,cr));
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cm);
				auto rr = ch.noise(cr);
				// auto rm=cm;
				// auto rr=cr;

				auto Lrm = ch.likelihood<float>(rm);
				auto Lrr = ch.likelihood<float>(rr);
				auto LLR = code::concatenate(code::DNAS::differential<ATGC>::decode_p(Lrm), code::DNAS::modifiedVLRLL<ATGC>::decode_p(Lrr, Lrm.back()));

				auto LLRest = ldpc.decode<decltype(ldpc)::DecoderType::SumProduct>(LLR);
				// auto LLRest = LLR;
				auto test = code::estimate_crop<SOURCE_LENGTH>(LLRest);
				auto cest = code::DNAS::differential<ATGC>::encode(test);
				auto mest = code::DNAS::VLRLL<ATGC>::decode(cest);

				nterror[n] += code::DNAS::countDifferentialError(cm,cest);
				biterror[n] += ((mest&mmask)^(m&mmask)).count();
				bitcount[n] += mmask.count();
			}
		}
	};

	auto aggregate = [&stat, &stats, repeat_per_thread](std::size_t dest){
		for(auto &st: stats) for(size_t n=0, nend=nsize; n<nend; ++n){
			std::get<0>(stat[dest])[n] += std::get<0>(st)[n];
			std::get<1>(stat[dest])[n] += std::get<1>(st)[n];
			std::get<2>(stat[dest])[n] += std::get<2>(st)[n];
			if(std::get<4>(stat[dest])<std::get<4>(st)) std::get<4>(stat[dest])=std::get<4>(st);
		}
		auto sum = 0.0, sqsum = 0.0;
		for(auto &st: stats) for(auto &pn: std::get<3>(st)) for(auto &pr: pn){
			sum += pr;
			sqsum += pr*pr;
		}
		auto num = static_cast<double>(nsize*NUM_THREADS*repeat_per_thread);
		auto average = sum/num;
		std::get<3>(stat[dest]).first = average;
		std::get<3>(stat[dest]).second = sqsum/num - average*average;
	};

	auto result = [&stat, repeat_per_thread](std::size_t target, std::size_t info_size){
		cout<<"GCper var: "<<std::get<3>(stat[target]).second<<", ave: "<<std::get<3>(stat[target]).first<<endl;
		cout<<"Run length: "<<std::get<4>(stat[target])<<endl;
		cout<<"Noise factor"
		<<"\tBER"
		<<"\tNER"
		<<endl;
		for(size_t n=0; n<noise_factor.size(); n++){
			cout<<noise_factor[n]
			<<"\t"<<static_cast<double>(std::get<0>(stat[target])[n])/static_cast<double>(std::get<1>(stat[target])[n])
			<<"\t"<<static_cast<double>(std::get<2>(stat[target])[n])/static_cast<double>(info_size/2*NUM_THREADS*repeat_per_thread)
			<<endl;
		}
	};

	vector<std::thread> threads;
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(0);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_conv, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(1);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_diff, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(2);
	tk.stop();

	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain"<<endl;
	result(0,SOURCE_LENGTH);

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded(conv)"<<endl;
	result(1,SOURCE_LENGTH);

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded(diff)"<<endl;
	result(2,SOURCE_LENGTH);

	return 0;
}
