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

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 5000;
constexpr size_t SOURCE_LENGTH = 512;
constexpr size_t CODE_LENGTH = 1024;
constexpr size_t NUM_THREADS = 20;
constexpr size_t BLOCK_SIZE = 32;
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
	using stattype = tuple<array<uint64_t,nsize>,array<uint64_t,nsize>,array<uint64_t,nsize>,size_t,vector<uint64_t>,vector<uint64_t>>;
	array<stattype,2> stat = {};
	array<stattype,NUM_THREADS> stats = {};

	// auto plain = [repeat_per_thread](size_t t, stattype *st){
	// 	auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
	// 	auto &maxrunlength = std::get<3>(*st);
	// 	auto &gcdistbef = std::get<4>(*st), &gcdistaft = std::get<5>(*st);
	// 	gcdistbef.resize((BLOCK_SIZE==0?SOURCE_LENGTH/2:BLOCK_SIZE)+1);
	// 	gcdistaft.resize((BLOCK_SIZE==0?SOURCE_LENGTH/2:BLOCK_SIZE)+1);
	// 	for(size_t n=0; n<nsize; ++n){
	// 		bitset<SOURCE_LENGTH> m;
	// 		util::RandomBits rb(t);
	// 		channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);

	// 		for(size_t r=0u; r<repeat_per_thread; r++){
	// 			rb.generate(m);

	// 			auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);

	// 			auto cmbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,1>::balance(cm);
	// 			// auto bm = cm;

	// 			code::DNAS::countBlockGC<BLOCK_SIZE>(cm, gcdistbef);
	// 			code::DNAS::countBlockGC<BLOCK_SIZE>(cmbar, gcdistaft);
	// 			auto runlength = code::DNAS::countRunlength(cmbar);
	// 			if(maxrunlength<runlength) maxrunlength = runlength;

	// 			auto rm = ch.noise(cmbar);
	// 			// auto rm=cm;

	// 			auto mest = code::DNAS::VLRLL<ATGC>::decode(rm);
	// 			nterror[n] += code::DNAS::countDifferentialError(cm,rm);
	// 			biterror[n] += ((mest&mmask)^(m&mmask)).count();
	// 			bitcount[n] += mmask.count();
	// 		}
	// 	}
	// };

	// auto encoded = [&ldpc, &decodertype, repeat_per_thread](size_t t, stattype *st){
	// 	auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
	// 	auto &maxrunlength = std::get<3>(*st);
	// 	auto &gcdistbef = std::get<4>(*st), &gcdistaft = std::get<5>(*st);
	// 	gcdistbef.resize((BLOCK_SIZE==0?CODE_LENGTH/2:BLOCK_SIZE)+1);
	// 	gcdistaft.resize((BLOCK_SIZE==0?CODE_LENGTH/2:BLOCK_SIZE)+1);
	// 	for(size_t n=0; n<nsize; ++n){
	// 		bitset<SOURCE_LENGTH> m;
	// 		util::RandomBits rb(t);
	// 		channel::NanoporeSequencing<ATGC> ch(noise_factor[n],t);

	// 		for(size_t r=0u; r<repeat_per_thread; r++){
	// 			rb.generate(m);

	// 			auto [cm, mmask, run] = code::DNAS::VLRLL<ATGC>::encode(m);
	// 			auto cmd = code::DNAS::differential::decode(cm);
	// 			auto tm = code::DNAS::convert<ATGC>::nttype_to_binary(cmd);
	// 			auto tr = ldpc.encode_redundancy(tm);
	// 			auto crd = code::DNAS::convert<ATGC>::binary_to_nttype(tr);
	// 			auto cr = code::DNAS::puncturedRLL<ATGC>::encode(crd, cm.back(), run);

	// 			auto c = code::concatenate(cm, cr);
	// 			auto cbar = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,1>::balance(c);
	// 			auto [cmbar, crbar] = code::split<cmd.size()>(cbar);
	// 			// auto cmbar = cm;
	// 			// auto crbar = cr;

	// 			code::DNAS::countBlockGC<BLOCK_SIZE>(c, gcdistbef);
	// 			code::DNAS::countBlockGC<BLOCK_SIZE>(cbar, gcdistaft);
	// 			auto runlength = code::DNAS::countRunlength(cbar);
	// 			if(maxrunlength<runlength) maxrunlength = runlength;

	// 			auto rm = ch.noise(cmbar);
	// 			auto rr = ch.noise(crbar);
	// 			// auto rm=cmbar;
	// 			// auto rr=crbar;

	// 			auto Lcmbar = ch.likelihood<float>(rm);
	// 			auto Lcrbar = ch.likelihood<float>(rr);
	// 			auto Lcmdbar = code::DNAS::differential::decode_p(Lcmbar);
	// 			auto Lcrdbar = code::DNAS::puncturedRLL<ATGC>::decode_p(Lcrbar, Lcmbar.back());
	// 			auto Lcbar = code::concatenate(Lcmdbar, Lcrdbar);
	// 			auto [Lc, p1] = code::DNAS::DivisionBalancing<ATGC,BLOCK_SIZE,1>::restore_p(Lcbar);
	// 			auto [Lcmd,Lcrd] = code::split<Lcmdbar.size()>(Lc);
	// 			auto LLR = code::concatenate(code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcmd), code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcrd));

	// 			auto LLRest = ldpc.decode(LLR,decodertype);
	// 			// auto LLRest = LLR;
	// 			auto tmest = code::estimate_crop<SOURCE_LENGTH>(LLRest);
	// 			auto cmdest = code::DNAS::convert<ATGC>::binary_to_nttype(tmest);
	// 			auto cmest = code::DNAS::differential::encode(cmdest);
	// 			auto mest = code::DNAS::VLRLL<ATGC>::decode(cmest);

	// 			nterror[n] += code::DNAS::countDifferentialError(cm,cmest);
	// 			biterror[n] += ((mest&mmask)^(m&mmask)).count();
	// 			bitcount[n] += mmask.count();
	// 		}
	// 	}
	// };

	auto plain_pitch = [repeat_per_thread](size_t t, stattype *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &maxrunlength = std::get<3>(*st);
		auto &gcdistbef = std::get<4>(*st), &gcdistaft = std::get<5>(*st);
		gcdistbef.resize((BLOCK_SIZE==0?SOURCE_LENGTH/2:BLOCK_SIZE)+1);
		gcdistaft.resize((BLOCK_SIZE==0?SOURCE_LENGTH/2:BLOCK_SIZE)+1);
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

				code::DNAS::countBlockGC<BLOCK_SIZE>(cm, gcdistbef);
				code::DNAS::countBlockGC<BLOCK_SIZE>(cmbar, gcdistaft);
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

	auto encoded_pitch = [&ldpc, &decodertype, repeat_per_thread](size_t t, stattype *st){
		auto &biterror = std::get<0>(*st), &bitcount = std::get<1>(*st), &nterror = std::get<2>(*st);
		auto &maxrunlength = std::get<3>(*st);
		auto &gcdistbef = std::get<4>(*st), &gcdistaft = std::get<5>(*st);
		gcdistbef.resize((BLOCK_SIZE==0?CODE_LENGTH/2:BLOCK_SIZE)+1);
		gcdistaft.resize((BLOCK_SIZE==0?CODE_LENGTH/2:BLOCK_SIZE)+1);
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
				auto crd = code::DNAS::puncturedRLL<ATGC>::encode(tr, run);
				auto cr = code::DNAS::differential::encode(crd, cm.back());
				// auto crd = code::DNAS::convert<ATGC>::binary_to_nttype(tr);
				// auto cr = code::DNAS::puncturedRLL<ATGC>::encode(crd, cm.back(), run);

				auto c = code::concatenate(cm, cr);
				auto cbar = bl.balance(c);
				auto [cmbar, crbar] = code::split<cmd.size()>(cbar);
				// auto cmbar = cm;
				// auto crbar = cr;

				code::DNAS::countBlockGC<BLOCK_SIZE>(c, gcdistbef);
				code::DNAS::countBlockGC<BLOCK_SIZE>(cbar, gcdistaft);
				auto runlength = code::DNAS::countRunlength(cbar);
				if(maxrunlength<runlength) maxrunlength = runlength;

				auto rm = ch.noise(cmbar);
				auto rr = ch.noise(crbar);
				// auto rm=cmbar;
				// auto rr=crbar;

				auto Lcmbar = ch.likelihood<float>(rm);
				auto Lcrbar = ch.likelihood<float>(rr);
				auto Lcmdbar = code::DNAS::differential::decode_p(Lcmbar);
				auto Lcrdbar = code::DNAS::differential::decode_p(Lcrbar, Lcmbar.back());
				// auto Lcrdbar = code::DNAS::puncturedRLL<ATGC>::decode_p(Lcrbar, Lcmbar.back());
				auto [Lcmd, p1] = bl.restore_p(Lcmdbar);
				auto [Lcrd, p2] = bl.restore_p(Lcrdbar,p1);
				auto LLR = code::concatenate(code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcmd), code::DNAS::puncturedRLL<ATGC>::decode_p(Lcrd));
				// auto LLR = code::concatenate(code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcmd), code::DNAS::convert<ATGC>::nttype_to_binary_p(Lcrd));

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
			std::get<4>(stat[dest]).resize(std::get<4>(st).size());
			for(size_t i=0, iend=std::get<4>(st).size(); i<iend; ++i) std::get<4>(stat[dest])[i] += std::get<4>(st)[i];
			std::get<5>(stat[dest]).resize(std::get<5>(st).size());
			for(size_t i=0, iend=std::get<5>(st).size(); i<iend; ++i) std::get<5>(stat[dest])[i] += std::get<5>(st)[i];
		}
	};

	auto result = [&stat, repeat_per_thread](std::size_t target, std::size_t channel_size){
		double num_blocks = static_cast<double>((BLOCK_SIZE==0?1:channel_size/2/BLOCK_SIZE)*nsize*NUM_THREADS*repeat_per_thread);
		cout<<"Run length: "<<std::get<3>(stat[target])<<endl;
		cout<<"Noise factor"
		<<"\tBER"
		<<"\tNER"
		<<endl;
		for(size_t n=0; n<noise_factor.size(); n++){
			cout<<noise_factor[n]
			<<"\t"<<static_cast<double>(std::get<0>(stat[target])[n])/static_cast<double>(std::get<1>(stat[target])[n])
			<<"\t"<<static_cast<double>(std::get<2>(stat[target])[n])/static_cast<double>(SOURCE_LENGTH/2*NUM_THREADS*repeat_per_thread)
			<<endl;
		}
		cout<<"GCcontent distribution(before): "<<endl;
		for(auto i:std::get<4>(stat[target])) cout<<static_cast<double>(i)/num_blocks<<"\t"<<flush;
		cout<<endl;
		cout<<"GCcontent distribution(after): "<<endl;
		for(auto i:std::get<5>(stat[target])) cout<<static_cast<double>(i)/num_blocks<<"\t"<<flush;
		cout<<endl;
	};

	vector<std::thread> threads;
	// tk.split();
	// for(stats = {}; auto &st: stats) threads.emplace_back(plain, threads.size(), &st);
	// for(auto &t: threads) t.join();
	// aggregate(0);
	// threads.clear();
	// tk.split();
	// for(stats = {}; auto &st: stats) threads.emplace_back(encoded, threads.size(), &st);
	// for(auto &t: threads) t.join();
	// aggregate(1);
	// threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(plain_pitch, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(0);
	threads.clear();
	tk.split();
	for(stats = {}; auto &st: stats) threads.emplace_back(encoded_pitch, threads.size(), &st);
	for(auto &t: threads) t.join();
	aggregate(1);
	tk.stop();

	cout<<"Block Size: "<<BLOCK_SIZE<<endl;
	// cout<<SOURCE_LENGTH<<endl;
	// cout<<"plain"<<endl;
	// result(0,SOURCE_LENGTH);
	// cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	// cout<<"encoded"<<endl;
	// result(1,SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<endl;
	cout<<"plain(pitch)"<<endl;
	result(0, SOURCE_LENGTH);
	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<"encoded(pitch)"<<endl;
	result(1, CODE_LENGTH);

	return 0;
}
