#include <array>
// #include <vector>
// #include <string>
// #include <cstdint>
#include "common/timekeep.hpp"
#include "common/randombits.hpp"
// #include "main/rangetools.hpp"
#include "idec/codeIDEC.hpp"
#include "idec/IDS.hpp"

using std::array, std::vector;
using std::int32_t, std::uint64_t, std::size_t;
using std::cout, std::cerr, std::flush, std::endl;

constexpr size_t DEFAULT_REPEAT_PER_THREAD = 1000u;
constexpr auto MARKER = {false};
constexpr size_t SOURCE_LENGTH = 512u;
constexpr size_t INTERVAL = 10u;
constexpr int32_t DRIFT_RANGE = 8;
constexpr double FLIP_RATE = 0.0;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>1){
		temp = std::stoul(argv[1]);
	}
	const auto repeat_per_thread = temp;
	cout<<SOURCE_LENGTH<<"("<<INTERVAL<<")"<<"*"<<repeat_per_thread/*<<"*"<<NUM_THREADS*/<<endl;

	constexpr array indel_rate = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05};

	array<uint64_t, indel_rate.size()> pl_biterror = {0};
	array<uint64_t, indel_rate.size()> ma_biterror = {0};
	// vector<uint64_t> pl_biterror(indel_rate.size(), 0);
	// vector<uint64_t> ma_biterror(indel_rate.size(), 0);
	const vector<bool> marker = MARKER;

	code::Marker ma(SOURCE_LENGTH, MARKER, INTERVAL, DRIFT_RANGE);

	tk.split();

	for(size_t n=0u, nend=indel_rate.size(); n<nend; ++n){
		vector<bool> source(SOURCE_LENGTH);
		auto &bn = pl_biterror[n];

		channel::IDS ch(indel_rate[n], DRIFT_RANGE, FLIP_RATE);
		util::RandomBits rb;

		for(size_t r=0u, rend=repeat_per_thread; r<rend; ++r){
			rb.generate(source);

			auto y = ch.noise(source);

			y.resize(source.size());
			for(size_t i=0u, iend=source.size(); i<iend; ++i) bn += source[i]^y[i];
		}
	}

	tk.split();

	auto decoder = ma.create_decoder();
	for(size_t n=0u, nend=indel_rate.size(); n<nend; ++n){
		vector<bool> source(SOURCE_LENGTH);
		auto &bn = ma_biterror[n];

		channel::IDS ch(indel_rate[n], DRIFT_RANGE, FLIP_RATE);
		util::RandomBits rb;

		for(size_t r=0u, rend=repeat_per_thread; r<rend; ++r){
			rb.generate(source);

			auto code = ma.encode(source);
			auto y = ch.noise(code);
			auto est = ma.decode(y, ch.insertion_rate(), ch.deletion_rate(), ch.substitution_rate(), decoder);

			for(size_t i=0u, iend=source.size(); i<iend; ++i) bn += source[i]^est[i];
		}
	}

	tk.stop();

	cout<<"Pid"
	<<"\tPlain"
	<<"\tMarker"
	<<endl;

	for(size_t n=0, nend=indel_rate.size(); n<nend; ++n){
		cout<<indel_rate[n]
		<<"\t"<<static_cast<double>(pl_biterror[n])/(SOURCE_LENGTH*repeat_per_thread)
		<<"\t"<<static_cast<double>(ma_biterror[n])/(SOURCE_LENGTH*repeat_per_thread)
		<<endl;
	}
	return 0;
}
