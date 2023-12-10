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
constexpr size_t BLOCK_SIZE = 8;
constexpr uint8_t ATGC = 0x1B;
constexpr double TOLERANCE = 0.125;

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
	// constexpr array noise_factor = {0.04,0.03,0.02,0.01,0.0};
	constexpr array noise_factor = {0.04,0.035,0.03,0.025,0.02,0.015,0.01,0.005,0.0};
	constexpr size_t nsize = noise_factor.size();

	auto ldpc = code::make_SystematicLDPC<SOURCE_LENGTH,CODE_LENGTH>();
	// tuple: biterrors, nterrors, maxGCdeviation
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t>,6> stat = {};
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t>,NUM_THREADS> stats = {};

	auto plain = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &maxgcdev = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto nm = code::DNAS::convert<ATGC>::binary_to_nttype(m);
				auto dm = code::DNAS::differential::encode(nm);

				auto dmbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,0>::balance(dm);

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(dmbar);
				if(dev>maxgcdev) maxgcdev=dev;

				auto rm = ch.noise(dmbar);
				// auto rm=cmbar;

				auto nmest = code::DNAS::differential::decode(rm);
				auto mest = code::DNAS::convert<ATGC>::nttype_to_binary(nmest);

				nterror[n] += code::DNAS::countDifferentialError(dm,rm);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto encoded = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &maxgcdev = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto c = ldpc.encode(m);
				auto nc = code::DNAS::convert<ATGC>::binary_to_nttype(c);
				auto dc = code::DNAS::differential::encode(nc);

				auto dcbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,0>::balance(dc);

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(dcbar);
				if(dev>maxgcdev) maxgcdev=dev;

				auto rc = ch.noise(dcbar);
				// auto rc=dcbar;

				auto Ldcbar = ch.likelihood(rc);
				auto Lncbar = code::DNAS::differential::decode_p(Ldcbar);
				auto [Lnc,prev] = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,0>::restore_p(Lncbar);
				auto Lc = code::DNAS::convert<ATGC>::nttype_to_binary_p(Lnc);

				auto Lcest = ldpc.decode<decltype(ldpc)::DecoderType::SumProduct>(Lc);
				auto mest = code::estimate_crop<SOURCE_LENGTH>(Lcest);

				nterror[n] += code::DNAS::countDifferentialError(dc,rc);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto plain_lesschange = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &maxgcdev = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,2> bl(TOLERANCE);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto nm = code::DNAS::convert<ATGC>::binary_to_nttype(m);
				auto dm = code::DNAS::differential::encode(nm);

				auto dmbar = bl.balance(dm);

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(dmbar);
				if(dev>maxgcdev) maxgcdev=dev;

				auto rm = ch.noise(dmbar);
				// auto rm=cmbar;

				auto nmest = code::DNAS::differential::decode(rm);
				auto mest = code::DNAS::convert<ATGC>::nttype_to_binary(nmest);

				nterror[n] += code::DNAS::countDifferentialError(dm,rm);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto encoded_lesschange = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &maxgcdev = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,2> bl(TOLERANCE);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto c = ldpc.encode(m);
				auto nc = code::DNAS::convert<ATGC>::binary_to_nttype(c);
				auto dc = code::DNAS::differential::encode(nc);

				auto dcbar = bl.balance(dc);

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(dcbar);
				if(dev>maxgcdev) maxgcdev=dev;

				auto rc = ch.noise(dcbar);
				// auto rc=dcbar;

				auto Ldcbar = ch.likelihood(rc);
				auto Lncbar = code::DNAS::differential::decode_p(Ldcbar);
				auto [Lnc,prev] = bl.restore_p(Lncbar);
				auto Lc = code::DNAS::convert<ATGC>::nttype_to_binary_p(Lnc);

				auto Lcest = ldpc.decode<decltype(ldpc)::DecoderType::SumProduct>(Lc);
				auto mest = code::estimate_crop<SOURCE_LENGTH>(Lcest);

				nterror[n] += code::DNAS::countDifferentialError(dc,rc);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto plain_pitch = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &maxgcdev = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,6> bl(TOLERANCE);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto nm = code::DNAS::convert<ATGC>::binary_to_nttype(m);
				auto dm = code::DNAS::differential::encode(nm);

				auto dmbar = bl.balance(dm);

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(dmbar);
				if(dev>maxgcdev) maxgcdev=dev;

				auto rm = ch.noise(dmbar);
				// auto rm=cmbar;

				auto nmest = code::DNAS::differential::decode(rm);
				auto mest = code::DNAS::convert<ATGC>::nttype_to_binary(nmest);

				nterror[n] += code::DNAS::countDifferentialError(dm,rm);
				biterror[n] += (mest^m).count();
			}
		}
	};

	auto encoded_pitch = [&ldpc, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t> *st){
		auto &biterror = std::get<0>(*st), &nterror = std::get<1>(*st);
		auto &maxgcdev = std::get<2>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,6> bl(TOLERANCE);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto c = ldpc.encode(m);
				auto nc = code::DNAS::convert<ATGC>::binary_to_nttype(c);
				auto dc = code::DNAS::differential::encode(nc);

				auto dcbar = bl.balance(dc);

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(dcbar);
				if(dev>maxgcdev) maxgcdev=dev;

				auto rc = ch.noise(dcbar);
				// auto rc=dcbar;

				auto Ldcbar = ch.likelihood(rc);
				auto Lncbar = code::DNAS::differential::decode_p(Ldcbar);
				auto [Lnc,prev] = bl.restore_p(Lncbar);
				auto Lc = code::DNAS::convert<ATGC>::nttype_to_binary_p(Lnc);

				auto Lcest = ldpc.decode<decltype(ldpc)::DecoderType::SumProduct>(Lc);
				auto mest = code::estimate_crop<SOURCE_LENGTH>(Lcest);

				nterror[n] += code::DNAS::countDifferentialError(dc,rc);
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
			if(std::get<2>(st)>std::get<2>(stat[dest])) std::get<2>(stat[dest]) = std::get<2>(st);
		}
	};

	auto result = [&stat, repeat_per_thread](std::size_t target, std::size_t channel_size){
		const std::size_t block_size = BLOCK_SIZE==0?channel_size:BLOCK_SIZE;
		cout<<"max GCcontent deviation: "<<static_cast<double>(std::get<2>(stat[target]))/static_cast<double>(block_size)<<endl;
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
	for(stats = {}; auto &st: stats) threads.emplace_back(plain_lesschange, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(2);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_lesschange, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(3);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain_pitch, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(4);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_pitch, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(5);
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
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain(pitch)"<<endl;
	result(4, SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded(pitch)"<<endl;
	result(5, CODE_LENGTH);

	return 0;
}
