#include <array>
// #include <vector>
// #include <string>
// #include <cstdint>
#include "common/timekeep.hpp"
#include "common/randombits.hpp"
// #include "common/rangetools.hpp"
#include "ldpcidec/codeLDPCIDEC.hpp"
#include "idec/IDS.hpp"

using std::array, std::vector;
using std::uint32_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;

constexpr auto DEFAULT_REPEAT_PER_THREAD = 1000ul;
constexpr auto MARKER = {false,true};
constexpr auto INTERVAL = 10ul;
constexpr auto DRIFT_RANGE = 8;
constexpr auto FLIP_RATE = 0.05;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	if(argc<2){
		cerr<<"Usage: "<<argv[0]<<" FILENAME"<<endl;
		return 1;
	}

	cout<<argv[1]<<endl;

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>2){
		temp = std::stoul(argv[2]);
	}
	const auto repeat_per_thread = temp;
	cout<<"("<<INTERVAL<<")"<<"*"<<repeat_per_thread/*<<"*"<<NUM_THREADS*/<<endl;

	constexpr array indel_rate = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05};
	// const vector<double> indel_rate = {0,0.01,0.02,0.03,0.04,0.05};

	array<uint64_t, indel_rate.size()> sp_biterror = {0};
	array<uint64_t, indel_rate.size()> mp_biterror = {0};
	// vector<uint64_t> sp_biterror(indel_rate.size(), 0);
	// vector<uint64_t> mp_biterror(indel_rate.size(), 0);
	const vector<bool> marker = MARKER;

	code::Marker_LDPC_SPD spd(argv[1], MARKER, INTERVAL, DRIFT_RANGE);
	code::Marker_LDPC_MPD mpd(argv[1], MARKER, INTERVAL, DRIFT_RANGE);

	tk.split();

	auto sp_decoder = spd.create_decoder<code::LDPC::SumProduct_Decoding>();
	for(size_t n=0u; n<indel_rate.size(); ++n){
	// for(auto [ri, mai]: util::zip(indel_rate, sp_biterror)){
		vector source(spd.sourcesize(), false);
		auto &bn = sp_biterror[n];

		channel::IDS ch(indel_rate[n], DRIFT_RANGE, FLIP_RATE);
		util::RandomBits rb;

		for(uint32_t r=0; r<repeat_per_thread; ++r){
			rb.generate(source);

			auto code = spd.encode(source);
			auto y = ch.noise(code);
			auto est = spd.decode(y, ch.insertion_rate(), ch.deletion_rate(), ch.substitution_rate(), sp_decoder);

			for(size_t i=0u, iend=source.size(); i<iend; ++i) bn += source[i]^est[i];
			// for(auto [si, ei]: util::zip(source, est)) mai += si^ei;
		}
	}

	tk.split();

	auto mp_decoder = mpd.create_decoder<code::LDPC::SumProduct_Decoding>();
	for(size_t n=0u; n<indel_rate.size(); ++n){
	// for(auto [ri, mai]: util::zip(indel_rate, sp_biterror)){
		vector source(mpd.sourcesize(), false);
		auto &bn = mp_biterror[n];

		channel::IDS ch(indel_rate[n], DRIFT_RANGE, FLIP_RATE);
		util::RandomBits rb;

		for(uint32_t r=0; r<repeat_per_thread; ++r){
			rb.generate(source);

			auto code = mpd.encode(source);
			auto y = ch.noise(code);
			auto est = mpd.decode(y, ch.insertion_rate(), ch.deletion_rate(), ch.substitution_rate(), mp_decoder);

			for(size_t i=0u, iend=source.size(); i<iend; ++i) bn += source[i]^est[i];
			// for(auto [si, ei]: util::zip(source, est)) mai += si^ei;
		}
	}

	tk.stop();

	cout<<"Pid"
	<<"\tMarker+Sum-Product[SPD]"
	<<"\tMarker+Sum-Product[MPD]"
	<<endl;

	for(uint32_t n=0; n<indel_rate.size(); ++n){
		cout<<indel_rate[n]
		<<"\t"<<(double)sp_biterror[n]/(double)(spd.sourcesize()*repeat_per_thread)
		<<"\t"<<(double)mp_biterror[n]/(double)(mpd.sourcesize()*repeat_per_thread)
		<<endl;
	}
	return 0;
}
