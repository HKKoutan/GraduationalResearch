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
constexpr size_t NUM_THREADS = 20;
constexpr size_t BLOCK_SIZE = 0;
constexpr std::uint8_t ATGC = 0x1B;
constexpr std::uint8_t ATGC2 = 0x27;

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
	constexpr array noise_factor = {0.04,0.03,0.02,0.01,0.0};
	// constexpr array noise_factor = {0.04,0.035,0.03,0.025,0.02,0.015,0.01,0.005,0.0};
	constexpr size_t nsize = noise_factor.size();
	code::LDPC::phi_table<> decodertype;

	auto ldpc = code::make_SystematicLDPC<SOURCE_LENGTH,CODE_LENGTH>();
	// tuple: biterrors, bitcounts, nterrors, maxrunlength, GCdist
	using stattype = tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,size_t,vector<uint64_t>>;
	array<stattype,3> stat = {};
	array<stattype,NUM_THREADS> stats = {};

	auto plain = [repeat_per_thread](size_t t, stattype *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &maxrunlength = std::get<3>(*st);
		auto &gcdist = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcdist.resize((BLOCK_SIZE==0?SOURCE_LENGTH/2:BLOCK_SIZE)+1);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);

				code::DNAS::countBlockGC<BLOCK_SIZE>(cm, gcdist);
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

	auto encoded_conv = [&ldpc, &decodertype, repeat_per_thread](size_t t, stattype *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &maxrunlength = std::get<3>(*st);
		auto &gcdist = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcdist.resize((BLOCK_SIZE==0?CODE_LENGTH/2:BLOCK_SIZE)+1);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);
				auto tm = code::DNAS::convert<ATGC>::nttype_to_binary(cm);
				auto tr = ldpc.encode_redundancy(tm);
				auto crd = code::DNAS::convert<ATGC>::binary_to_nttype(tr);
				auto cr = code::DNAS::puncturedRLL<ATGC>::encode(crd, cm.back(), run);

				auto c = code::concatenate(cm,cr);
				code::DNAS::countBlockGC<BLOCK_SIZE>(c, gcdist);
				auto runlength = code::DNAS::countRunlength(c);
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cm);
				auto rr = ch.noise(cr);
				// auto rm=cm;
				// auto rr=cr;

				auto Lcm = ch.likelihood<float>(rm);
				auto Lcr = ch.likelihood<float>(rr);
				auto Lcrd = code::DNAS::puncturedRLL<ATGC>::decode_p(Lcr, Lcm.back());
				auto Ltr = code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcrd);
				auto Ltm = code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcm);
				auto Ltc = code::concatenate(Ltm, Ltr);

				auto Ltcest = ldpc.decode(Ltc,decodertype);
				// auto LLRest = LLR;
				auto tmest = code::estimate_crop<SOURCE_LENGTH>(Ltcest);
				auto cmest = code::DNAS::convert<ATGC>::binary_to_nttype(tmest);
				auto mest = code::DNAS::VLRLL<ATGC>::decode(cmest);

				nterror[n] += code::DNAS::countError(cm,cmest);
				biterror[n] += ((mest&mmask)^(m&mmask)).count();
				bitcount[n] += mmask.count();
			}
		}
	};

	auto encoded_diff = [&ldpc, &decodertype, repeat_per_thread](size_t t, stattype *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &maxrunlength = std::get<3>(*st);
		auto &gcdist = std::get<4>(*st);
		for(size_t n=0; n<nsize; ++n){
			bitset<SOURCE_LENGTH> m;
			channel::NanoporeSequencing<ATGC2> ch(noise_factor[n],t);
			util::RandomBits rb(t);
			gcdist.resize((BLOCK_SIZE==0?CODE_LENGTH/2:BLOCK_SIZE)+1);

			for(size_t r=0u; r<repeat_per_thread; r++){
				rb.generate(m);

				auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC2>::encode(m);
				auto cmd = code::DNAS::differential::decode(cm);
				auto tm = code::DNAS::convert<ATGC2>::nttype_to_binary(cmd);
				auto tr = ldpc.encode_redundancy(tm);
				auto crd = code::DNAS::convert<ATGC2>::binary_to_nttype(tr);
				auto cr = code::DNAS::puncturedRLL<ATGC2>::encode(crd, cm.back(), run);

				auto c = code::concatenate(cm,cr);
				code::DNAS::countBlockGC<BLOCK_SIZE>(c, gcdist);
				auto runlength = code::DNAS::countRunlength(c);
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cm);
				auto rr = ch.noise(cr);
				// auto rm=cm;
				// auto rr=cr;

				auto Lcm = ch.likelihood<float>(rm);
				auto Lcr = ch.likelihood<float>(rr);
				auto Lcmd = code::DNAS::differential::decode_p(Lcm);
				auto Ltm = code::DNAS::convert<ATGC2>::nttype_to_binary_p(Lcmd);
				auto Lcrd = code::DNAS::puncturedRLL<ATGC2>::decode_p(Lcr, Lcm.back());
				auto Ltr = code::DNAS::convert<ATGC2>::nttype_to_binary_p(Lcrd);
				auto Ltc = code::concatenate(Ltm, Ltr);

				auto Ltcest = ldpc.decode(Ltc,decodertype);
				// auto LLRest = LLR;
				auto tmest = code::estimate_crop<SOURCE_LENGTH>(Ltcest);
				auto cmdest = code::DNAS::convert<ATGC2>::binary_to_nttype(tmest);
				auto cmest = code::DNAS::differential::encode(cmdest);
				auto mest = code::DNAS::VLRLL<ATGC2>::decode(cmest);

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
			std::get<4>(stat[dest]).resize(std::get<4>(st).size());
			for(size_t i=0, iend=std::get<4>(st).size(); i<iend; ++i) std::get<4>(stat[dest])[i] += std::get<4>(st)[i];
		}
	};

	auto result = [&stat, repeat_per_thread](std::size_t target, std::size_t info_size){
		cout<<"Run length: "<<std::get<3>(stat[target])<<endl;
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
		cout<<"GCcontent distribution: "<<endl;
		for(auto i:std::get<4>(stat[target])) cout<<i<<"\t"<<flush;
		cout<<endl;
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

	cout<<"Block Size: "<<BLOCK_SIZE<<endl;
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
