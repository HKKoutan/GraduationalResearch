#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <thread>
#include <numeric>
#include "dnas/DNASmytype.hpp"
#include "dnas/codeDNAS.hpp"
#include "dnas/sequencer.hpp"
#include "ldpc/codeLDPC.hpp"
#include "main/randombits.hpp"
#include "main/timekeep.hpp"

using std::vector, std::string;
using std::int32_t, std::uint32_t, std::uint64_t;
using std::cout, std::cerr, std::flush, std::endl;
using code::DNAS::nucleotide_t;

constexpr uint32_t REPEAT_PER_THREAD = 5000;
constexpr uint32_t NUM_THREADS = 12;

int main(int argc, char* argv[]){
	util::Timekeep tk;
	tk.start();

	if(argc<2){
		cerr<<"Usage: "<<argv[0]<<" FILENAME"<<endl;
		return 1;
	}

	string path(argv[1]);
	cout<<path<<endl;

	uint32_t temp = REPEAT_PER_THREAD;
	if(argc>2){
		temp = std::stoi(string{argv[2]});
	}
	const uint32_t repeat_per_thread = temp;
	cout<<repeat_per_thread<<"*"<<NUM_THREADS<<endl;

	//const vector<double> noise_factor = {0.03};
	//const vector<double> noise_factor = {0.04,0.0375,0.035,0.0325,0.03};
	const vector<double> noise_factor = {0.04,0.035,0.03,0.025,0.02};

	vector<std::thread> threads;
	vector<vector<uint64_t>> biterrors(NUM_THREADS, vector<uint64_t>(noise_factor.size(), 0));
	vector<vector<uint64_t>> bitcounts(NUM_THREADS, vector<uint64_t>(noise_factor.size(), 0));
	vector<vector<uint64_t>> nterrors(NUM_THREADS, vector<uint64_t>(noise_factor.size(), 0));
	vector<vector<uint64_t>> ntcounts(NUM_THREADS, vector<uint64_t>(noise_factor.size(), 0));
	vector<vector<double>> GCpercentages(NUM_THREADS, vector<double>(repeat_per_thread*noise_factor.size()));

	code::Systematic_LDPC ldpc(path, NUM_THREADS);

	tk.split();

	for(uint32_t t=0; t<NUM_THREADS; t++){
		threads.emplace_back([&ldpc, &biterrors, &bitcounts, &nterrors, &ntcounts, &GCpercentages, &noise_factor, repeat_per_thread, t](){
			for(uint32_t n=0; n<noise_factor.size(); n++){
				vector<bool> m(ldpc.sourcesize());
				vector<nucleotide_t> cm(ldpc.sourcesize()/2);
				vector<bool> tm;
				vector<bool> tr;
				vector<nucleotide_t> cr;
				uint32_t qty_AT, qty_GC, qty_ATs, qty_GCs;
				vector<double> LLR;
				vector<double> LLRr;
				vector<bool> test;
				vector<nucleotide_t> cest;
				vector<bool> mest;

				channel::Nanopore_Sequencing ch(noise_factor[n],t);
				util::RandomBits rb(t);

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
					qty_ATs = qty_AT;
					qty_GCs = qty_GC;
					code::DNAS::nt_qty_count(cr, qty_AT, qty_GC);
					qty_ATs += qty_AT;
					qty_GCs += qty_GC;
					GCpercentages[t][n*repeat_per_thread + r] = double(qty_GCs)/double(qty_ATs+qty_GCs);

					vector<nucleotide_t> rm(cm);

					ch.noise(rm);
					ch.noise(cr);

					ch.message_LLR(rm,LLR);
					ch.redundancy_LLR(cr,LLRr,cm[cm.size()-1]);
					LLR.insert(LLR.end(),LLRr.begin(),LLRr.end());

					test = ldpc.sumproduct_decode(LLR,t);
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

	vector<uint64_t> biterror(noise_factor.size(), 0);
	vector<uint64_t> bitcount(noise_factor.size(), 0);
	vector<uint64_t> nterror(noise_factor.size(), 0);
	vector<uint64_t> ntcount(noise_factor.size(), 0);

	for(uint32_t n=0; n<noise_factor.size(); n++){
		for(uint32_t t=0; t<NUM_THREADS; t++){
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

	for(int32_t n=0; n<noise_factor.size(); n++){
		cout<<noise_factor[n]
		<<"\t"<<(double)biterror[n]/(double)(bitcount[n])
		<<"\t"<<(double)nterror[n]/(double)(ntcount[n])
		<<endl;
	}
	return 0;
}
