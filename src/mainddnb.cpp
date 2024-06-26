﻿#include <array>
#include <tuple>
#include <thread>
#include <numeric>
#include "dnas/codeDNASadapter.hpp"
#include "dnas/codeDNASstats.hpp"
#include "dnas/channelsequencer.hpp"
#include "ldpc/codeSystematicLDPC.hpp"
#include "common/codecommon.hpp"
#include "common/util.hpp"

using std::array, std::bitset, std::vector, std::tuple, std::pair;
using std::size_t, std::uint8_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;
using code::DNAS::nucleotide_t;

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 5000;
constexpr size_t SOURCE_LENGTH = 512;
constexpr size_t CODE_LENGTH = 1024;
constexpr size_t NUM_THREADS = 20;
constexpr size_t BLOCK_SIZE = 16;//GCdevの集計にのみ影響
constexpr uint8_t ATGC = 0x1B;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	cout<<"Title: DNA storage simulation on nanopore sequencing channel with simple conversion or differential encoding."<<endl;
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
	code::LDPC::phi_table<> decodertype;

	auto ldpc = code::make_SystematicLDPC<SOURCE_LENGTH,CODE_LENGTH>();
	// tuple: biterrors, nterrors, maxGCdeviation
	using stattype = tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,vector<uint64_t>>;
	array<stattype,4> stat = {};
	array<stattype,NUM_THREADS> stats = {};

	auto plain_conv = [repeat_per_thread](size_t t, stattype *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcdist = std::get<2>(*st);
		gcdist.resize((BLOCK_SIZE==0?SOURCE_LENGTH/2:BLOCK_SIZE)+1);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			util::RandomBits rb(t);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto cm = code::DNAS::convert<ATGC>::binary_to_nttype(m);
				// for(auto &ci: cm) ci = 3;

				code::DNAS::countBlockGC<BLOCK_SIZE>(cm, gcdist);

				auto rm = ch.noise(cm);
				// auto rm=cm;

				auto mest = code::DNAS::convert<ATGC>::nttype_to_binary(rm);

				nterror[n] += code::DNAS::countError(cm,rm);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto encoded_conv = [&ldpc, &decodertype, repeat_per_thread](size_t t, stattype *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcdist = std::get<2>(*st);
		gcdist.resize((BLOCK_SIZE==0?CODE_LENGTH/2:BLOCK_SIZE)+1);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			util::RandomBits rb(t);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto c = ldpc.encode(m);
				auto cc = code::DNAS::convert<ATGC>::binary_to_nttype(c);

				code::DNAS::countBlockGC<BLOCK_SIZE>(cc, gcdist);

				auto rc = ch.noise(cc);
				// auto rc=cc;

				auto Lcc = ch.likelihood<float>(rc);
				auto Lc = code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcc);

				auto Lcest = ldpc.decode(Lc,decodertype);
				auto mest = code::estimate_crop<SOURCE_LENGTH>(Lcest);

				nterror[n] += code::DNAS::countError(cc,rc);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto plain_diff = [repeat_per_thread](size_t t, stattype *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcdist = std::get<2>(*st);
		gcdist.resize((BLOCK_SIZE==0?SOURCE_LENGTH/2:BLOCK_SIZE)+1);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			util::RandomBits rb(t);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto nm = code::DNAS::convert<ATGC>::binary_to_nttype(m);
				auto cm = code::DNAS::differential::encode(nm);
				// for(auto &ci: cm) ci = 3;

				code::DNAS::countBlockGC<BLOCK_SIZE>(cm, gcdist);

				auto rm = ch.noise(cm);
				// auto rm=cm;

				auto nmest = code::DNAS::differential::decode(rm);
				auto mest = code::DNAS::convert<ATGC>::nttype_to_binary(nmest);

				nterror[n] += code::DNAS::countDifferentialError(cm,rm);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto encoded_diff = [&ldpc, &decodertype, repeat_per_thread](size_t t, stattype *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &gcdist = std::get<2>(*st);
		gcdist.resize((BLOCK_SIZE==0?CODE_LENGTH/2:BLOCK_SIZE)+1);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			util::RandomBits rb(t);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto c = ldpc.encode(m);
				auto nc = code::DNAS::convert<ATGC>::binary_to_nttype(c);
				auto cc = code::DNAS::differential::encode(nc);

				code::DNAS::countBlockGC<BLOCK_SIZE>(cc, gcdist);

				auto rc = ch.noise(cc);
				// auto rc=cc;

				auto Lcc = ch.likelihood<float>(rc);
				auto Lnc = code::DNAS::differential::decode_p(Lcc);
				auto Lc = code::DNAS::convert<ATGC>::nttype_to_binary_p(Lnc);

				auto Lcest = ldpc.decode(Lc,decodertype);
				auto mest = code::estimate_crop<SOURCE_LENGTH>(Lcest);

				nterror[n] += code::DNAS::countDifferentialError(cc,rc);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto aggregate = [&stat, &stats, repeat_per_thread](std::size_t dest){
		for(auto &st: stats){
			for(size_t n=0, nend=nsize; n<nend; ++n){
				std::get<0>(stat[dest])[n] += std::get<0>(st)[n];
				std::get<1>(stat[dest])[n] += std::get<1>(st)[n];
			}
			std::get<2>(stat[dest]).resize(std::get<2>(st).size());
			for(size_t i=0, iend=std::get<2>(st).size(); i<iend; ++i) std::get<2>(stat[dest])[i] += std::get<2>(st)[i];
		}
	};

	auto result = [&stat, repeat_per_thread](std::size_t target, std::size_t channel_size){
		double num_blocks = static_cast<double>((BLOCK_SIZE==0?1:channel_size/2/BLOCK_SIZE)*nsize*NUM_THREADS*repeat_per_thread);
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
		cout<<"GCcontent distribution: "<<endl;
		for(auto i:std::get<2>(stat[target])) cout<<static_cast<double>(i)/num_blocks<<"\t"<<flush;
		cout<<endl;
	};

	//実行
	vector<std::thread> threads;
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain_conv, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(0);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_conv, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(1);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain_diff, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(2);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_diff, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(3);
	tk.stop();

	//結果表示
	cout<<"Block Size: "<<BLOCK_SIZE<<endl;
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain(conv)"<<endl;
	result(0, SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded(conv)"<<endl;
	result(1, CODE_LENGTH);
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain(diff)"<<endl;
	result(2, SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded(diff)"<<endl;
	result(3, CODE_LENGTH);

	return 0;
}
