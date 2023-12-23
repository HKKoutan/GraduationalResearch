#include <array>
#include <tuple>
#include <thread>
#include <numeric>
#include "dnas/DNASnttype.hpp"
#include "dnas/codeDNASadapter.hpp"
#include "dnas/codeDNASbalancing.hpp"
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
constexpr size_t BLOCK_SIZE = 8;
constexpr std::uint8_t ATGC = 0x27;
constexpr double TOLERANCE = 0.125;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	cout<<"Title: DNA storage simulation on nanopore sequencing channel with VLRLL encoding, differential encoding and division balancing."<<endl;
	cout<<"ATGC: "<<std::bitset<8>(ATGC)<<endl;

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>1) temp = std::stoull(argv[1]);
	const auto repeat_per_thread = temp;

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<repeat_per_thread<<"*"<<NUM_THREADS<<endl;

	// constexpr array noise_factor = {0};
	// constexpr array noise_factor = {0.04,0.03,0.02,0.01,0.0};
	constexpr array noise_factor = {0.04,0.035,0.03,0.025,0.02,0.015,0.01,0.005,0.0};
	constexpr size_t nsize = noise_factor.size();
	code::LDPC::phi_table<> decodertype;

	auto ldpc = code::make_SystematicLDPC<SOURCE_LENGTH,CODE_LENGTH>();
	// tuple: biterrors, bitcounts, nterrors, maxGCdeviation, maxrunlength
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t,size_t>,4> stat = {};
	array<tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t,size_t>,NUM_THREADS> stats = {};

	auto plain = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t,size_t> *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &maxgcdev = std::get<3>(*st);
		auto &maxrunlength = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);

				auto cmbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,1>::balance(cm);
				// auto bm = cm;

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(cmbar);
				if(dev>maxgcdev) maxgcdev=dev;
				auto runlength = code::DNAS::countRunlength(cmbar);
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cmbar);
				// auto rm=cm;

				auto mest = code::DNAS::VLRLL<ATGC>::decode(rm);
				nterror[n] += code::DNAS::countDifferentialError(cm,rm);
				biterror[n] += ((mest&mmask)^(m&mmask)).count();
				bitcount[n] += mmask.count();
			}
		}
	};

	auto encoded = [&ldpc, &decodertype, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t,size_t> *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &maxgcdev = std::get<3>(*st);
		auto &maxrunlength = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);
				auto cmd = code::DNAS::differential::decode(cm);
				auto tm = code::DNAS::convert<ATGC>::nttype_to_binary(cmd);
				auto tr = ldpc.encode_redundancy(tm);
				auto crd = code::DNAS::convert<ATGC>::binary_to_nttype(tr);
				auto cr = code::DNAS::puncturedRLL<ATGC>::encode(crd, cm.back(), run);

				auto c = code::concatenate(cm, cr);
				auto cbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,1>::balance(c);
				auto [cmbar, crbar] = code::split<cmd.size()>(cbar);
				// auto cmbar = cm;
				// auto crbar = cr;

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(cmbar);
				if(dev>maxgcdev) maxgcdev=dev;
				dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(crbar);
				if(dev>maxgcdev) maxgcdev=dev;
				auto runlength = code::DNAS::countRunlength(code::concatenate(cmbar,crbar));
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cmbar);
				auto rr = ch.noise(crbar);
				// auto rm=cmbar;
				// auto rr=crbar;

				auto Lcmbar = ch.likelihood<float>(rm);
				auto Lcrbar = ch.likelihood<float>(rr);
				auto Lcmdbar = code::DNAS::differential::decode_p(Lcmbar);
				auto Lcrdbar = code::DNAS::puncturedRLL<ATGC>::decode_p(Lcrbar, Lcmbar.back());
				auto Lcbar = code::concatenate(Lcmdbar, Lcrdbar);
				auto [Lc, p1] = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,1>::restore_p(Lcbar);
				auto [Lcmd,Lcrd] = code::split<Lcmdbar.size()>(Lc);
				auto LLR = code::concatenate(code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcmd), code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcrd));

				auto LLRest = ldpc.decode(LLR,decodertype);
				// auto LLRest = LLR;
				auto tmest = code::estimate_crop<SOURCE_LENGTH>(LLRest);
				auto cmdest = code::DNAS::convert<ATGC>::binary_to_nttype(tmest);
				auto cmest = code::DNAS::differential::encode(cmdest);
				auto mest = code::DNAS::VLRLL<ATGC>::decode(cmest);

				nterror[n] += code::DNAS::countDifferentialError(cm,cmest);
				biterror[n] += ((mest&mmask)^(m&mmask)).count();
				bitcount[n] += mmask.count();
			}
		}
	};

	auto plain_pitch = [repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t,size_t> *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &maxgcdev = std::get<3>(*st);
		auto &maxrunlength = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,7> bl(TOLERANCE);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);

				auto cmbar = bl.balance(cm);
				// auto bm = cm;

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(cmbar);
				if(dev>maxgcdev) maxgcdev=dev;
				auto runlength = code::DNAS::countRunlength(cmbar);
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cmbar);
				// auto rm=cm;

				auto mest = code::DNAS::VLRLL<ATGC>::decode(rm);
				nterror[n] += code::DNAS::countDifferentialError(cm,rm);
				biterror[n] += ((mest&mmask)^(m&mmask)).count();
				bitcount[n] += mmask.count();
			}
		}
	};

	auto encoded_pitch = [&ldpc, &decodertype, repeat_per_thread](size_t t, tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,uint64_t,size_t> *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &maxgcdev = std::get<3>(*st);
		auto &maxrunlength = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			util::RandomBits rb(t);
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,7> bl(TOLERANCE);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);
				auto cmd = code::DNAS::differential::decode(cm);
				auto tm = code::DNAS::convert<ATGC>::nttype_to_binary(cmd);
				auto tr = ldpc.encode_redundancy(tm);
				auto crd = code::DNAS::convert<ATGC>::binary_to_nttype(tr);
				auto cr = code::DNAS::puncturedRLL<ATGC>::encode(crd, cm.back(), run);

				auto c = code::concatenate(cm, cr);
				auto cbar = bl.balance(c);
				auto [cmbar, crbar] = code::split<cmd.size()>(cbar);
				// auto cmbar = cm;
				// auto crbar = cr;

				auto dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(cmbar);
				if(dev>maxgcdev) maxgcdev=dev;
				dev = code::DNAS::countBlockGCmaxDeviation<BLOCK_SIZE>(crbar);
				if(dev>maxgcdev) maxgcdev=dev;
				auto runlength = code::DNAS::countRunlength(code::concatenate(cmbar,crbar));
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cmbar);
				auto rr = ch.noise(crbar);
				// auto rm=cmbar;
				// auto rr=crbar;

				auto Lcmbar = ch.likelihood<float>(rm);
				auto Lcrbar = ch.likelihood<float>(rr);
				auto Lcmdbar = code::DNAS::differential::decode_p(Lcmbar);
				auto Lcrdbar = code::DNAS::puncturedRLL<ATGC>::decode_p(Lcrbar, Lcmbar.back());
				auto [Lcmd, p1] = bl.restore_p(Lcmdbar);
				auto [Lcrd, p2] = bl.restore_p(Lcrdbar,p1);
				auto LLR = code::concatenate(code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcmd), code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcrd));

				auto LLRest = ldpc.decode(LLR,decodertype);
				// auto LLRest = LLR;
				auto tmest = code::estimate_crop<SOURCE_LENGTH>(LLRest);
				auto cmdest = code::DNAS::convert<ATGC>::binary_to_nttype(tmest);
				auto cmest = code::DNAS::differential::encode(cmdest);
				auto mest = code::DNAS::VLRLL<ATGC>::decode(cmest);

				nterror[n] += code::DNAS::countDifferentialError(cm,cmest);
				biterror[n] += ((mest&mmask)^(m&mmask)).count();
				bitcount[n] += mmask.count();
			}
		}
	};

	auto aggregate = [&stat, &stats, repeat_per_thread](std::size_t dest){
		for(auto &st: stats){
			for(size_t n=0, nend=nsize; n<nend; ++n){
				std::get<0>(stat[dest])[n] += std::get<0>(st)[n];
				std::get<1>(stat[dest])[n] += std::get<1>(st)[n];
				std::get<2>(stat[dest])[n] += std::get<2>(st)[n];
			}
			if(std::get<3>(st)>std::get<3>(stat[dest])) std::get<3>(stat[dest])=std::get<3>(st);
			if(std::get<4>(st)>std::get<4>(stat[dest])) std::get<4>(stat[dest])=std::get<4>(st);
		}
	};

	auto result = [&stat, repeat_per_thread](std::size_t target, std::size_t info_size){
		const std::size_t block_size = BLOCK_SIZE==0?info_size:BLOCK_SIZE;
		cout<<"max GCcontent deviation: "<<static_cast<double>(std::get<3>(stat[target]))/static_cast<double>(block_size)<<endl;
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
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(1);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain_pitch, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(2);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_pitch, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(3);
	tk.stop();

	cout<<"Block Size: "<<BLOCK_SIZE<<endl;
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain"<<endl;
	result(0,SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded"<<endl;
	result(1,SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain(pitch)"<<endl;
	result(2,SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded(pitch)"<<endl;
	result(3,SOURCE_LENGTH);

	return 0;
}
