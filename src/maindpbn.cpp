#include <array>
#include <thread>
#include <numeric>
#include "dnas/DNASnttype.hpp"
#include "dnas/codeDNAS.hpp"
#include "dnas/sequencer.hpp"
#include "ldpc/codeLDPC.hpp"
#include "common/randombits.hpp"
#include "common/timekeep.hpp"

using std::array, std::bitset, std::vector;
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

	auto temp = DEFAULT_REPEAT_PER_THREAD;
	if(argc>1) temp = std::stoull(argv[1]);
	const auto repeat_per_thread = temp;

	cout<<SOURCE_LENGTH<<"->"<<CODE_LENGTH<<endl;
	cout<<repeat_per_thread<<"*"<<NUM_THREADS<<endl;

	//const vector<double> noise_factor = {0.03};
	//const vector<double> noise_factor = {0.04,0.0375,0.035,0.0325,0.03};
	constexpr array noise_factor = {0.04,0.035,0.03,0.025,0.02};

	vector<std::thread> threads;
	array<array<uint64_t,noise_factor.size()>,NUM_THREADS> biterrors = {0};
	array<array<uint64_t,noise_factor.size()>,NUM_THREADS> bitcounts = {0};
	array<array<uint64_t,noise_factor.size()>,NUM_THREADS> nterrors = {0};
	array<array<uint64_t,noise_factor.size()>,NUM_THREADS> ntcounts = {0};
	array<vector<double>,NUM_THREADS> GCpercentages;
	for(auto &i: GCpercentages) i.resize(repeat_per_thread*noise_factor.size());

	code::Systematic_LDPC<SOURCE_LENGTH,CODE_LENGTH> ldpc;

	tk.split();

	for(uint32_t t=0; t<NUM_THREADS; t++){
		threads.emplace_back([&ldpc, &biterrors, &bitcounts, &nterrors, &ntcounts, &GCpercentages, &noise_factor, repeat_per_thread, t](){
			for(uint32_t n=0; n<noise_factor.size(); n++){
				bitset<SOURCE_LENGTH> m;
				vector<nucleotide_t> cm(ldpc.sourcesize()/2);
				vector<bool> tm;
				vector<bool> tr;
				vector<nucleotide_t> cr;
				uint32_t qty_AT, qty_GC, qty_ATs, qty_GCs;
				vector<nucleotide_t> crbar;
				vector<bool> info;
				vector<bool> ir;
				vector<nucleotide_t> ci;
				vector<double> LLRi;
				vector<double> LLR;
				vector<double> LLRr;
				vector<bool> iest;
				vector<bool> test;
				vector<nucleotide_t> cest;
				vector<bool> mest;

				channel::Nanopore_Sequencing ch(noise_factor[n],t);
				util::RandomBits<SOURCE_LENGTH> rb(t);

				for(uint32_t r=0; r<repeat_per_thread; r++){
					rb.generate(m);

					const auto mend = code::DNAS::VLRLL_encode(m.begin(),m.end(),cm.begin(),cm.end());
					// cout<<"m:"<<m.size()<<endl;
					// for(auto i:m) cout<<i<<flush;
					// cout<<endl;

					code::DNAS::interim_map(cm,tm);
					tr = ldpc.encode(tm);
					code::DNAS::modified_VLRLL_encode(tr,cr,cm[cm.size()-1]);

					code::DNAS::nt_qty_count(cm, qty_AT, qty_GC);
					code::DNAS::nt_addequalizing_encode(cr,crbar,info,qty_AT,qty_GC);

					code::DNAS::modified_VLRLL_encode(info,ci,crbar[crbar.size()-1]);

					code::DNAS::nt_qty_count(cm, qty_AT, qty_GC);
					qty_ATs = qty_AT;
					qty_GCs = qty_GC;
					code::DNAS::nt_qty_count(crbar, qty_AT, qty_GC);
					qty_ATs += qty_AT;
					qty_GCs += qty_GC;
					code::DNAS::nt_qty_count(ci, qty_AT, qty_GC);
					qty_ATs += qty_AT;
					qty_GCs += qty_GC;
					GCpercentages[t][n*repeat_per_thread + r] = double(qty_GCs)/double(qty_ATs+qty_GCs);

					vector<nucleotide_t> rm(cm);

					ch.noise(rm);
					ch.noise(crbar);
					ch.noise(ci);

					code::DNAS::modified_VLRLL_decode(ci,info,crbar[crbar.size()-1]);

					code::DNAS::nt_addequalizing_decode(crbar,info,cr);

					ch.message_LLR(rm,LLR);
					ch.redundancy_LLR(cr,LLRr,rm[rm.size()-1]);
					LLR.insert(LLR.end(),LLRr.begin(),LLRr.end());

					test = ldpc.sumproduct_decode(LLR, t);
					code::DNAS::interim_demap(test,cest);
					code::DNAS::VLRLL_decode(cest,mest);

					for(uint64_t i=0; i<cm.size(); i++) nterrors[t][n] += (cm[i]!=cest[i]);
					ntcounts[t][n] += cm.size();

					auto mi = m.begin();
					auto ei = mest.begin();
					const auto eend = mest.end();
					size_t i;
					for(i=0; mi!=mend&&ei!=eend; i++) biterrors[t][n] += (*mi++)^(*ei++);
					bitcounts[t][n] += i;
					//cout<<nterrors[t][n]<<","<<biterrors[t][n]<<","<<ntcounts[t][n]<<","<<bitcounts[t][n]<<endl;
				}
			}
		});
	}
	for(auto &t: threads) t.join();

	tk.split();

	double sum = 0;
	for(auto &p: GCpercentages) sum = std::reduce(p.begin(),p.end(),sum);
	double average = sum/(noise_factor.size()*NUM_THREADS*repeat_per_thread);
	sum = 0;
	for(auto &p: GCpercentages) sum = std::inner_product(p.begin(),p.end(),p.begin(),sum);
	double var = sum/(noise_factor.size()*NUM_THREADS*repeat_per_thread) - (average*average);

	array<uint64_t,noise_factor.size()> biterror = {0};
	array<uint64_t,noise_factor.size()> bitcount = {0};
	array<uint64_t,noise_factor.size()> nterror = {0};
	array<uint64_t,noise_factor.size()> ntcount = {0};

	for(size_t n=0; n<noise_factor.size(); n++){
		for(size_t t=0; t<NUM_THREADS; t++){
			biterror[n] += biterrors[t][n];
			bitcount[n] += bitcounts[t][n];
			nterror[n] += nterrors[t][n];
			ntcount[n] += ntcounts[t][n];
		}
	}

	tk.stop();

	cout<<"var GCper: "<<var<<endl;

	cout<<"Noise factor"
	<<"\tBER"
	<<"\tNER"
	<<endl;

	for(size_t n=0; n<noise_factor.size(); n++){
		cout<<noise_factor[n]
		<<"\t"<<static_cast<double>(biterror[n])/static_cast<double>(bitcount[n])
		<<"\t"<<static_cast<double>(nterror[n])/static_cast<double>(ntcount[n])
		<<endl;
	}
	return 0;
}
