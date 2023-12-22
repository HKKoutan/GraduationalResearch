#include <iostream>
#include <array>
#include <bitset>
#include <vector>
#include <cstdint>
#include <cassert>
#include <thread>

constexpr std::size_t LENGTH = 16;
constexpr double TOLERANCE = 0.125;
constexpr std::size_t THREAD_PER_WORK_EXP = 2;
constexpr std::size_t THREAD_PER_WORK = 1<<THREAD_PER_WORK_EXP;

int main(){
	constexpr std::size_t divsize = LENGTH>>1;
	constexpr std::uint64_t qtytolerance = static_cast<std::uint64_t>(static_cast<double>(LENGTH)*TOLERANCE);
	constexpr std::uint64_t qtyhalflow = divsize-qtytolerance;
	constexpr std::uint64_t qtyhalfhigh = divsize+qtytolerance;

	std::array<std::uint64_t,LENGTH> count1 = {};
	std::array<std::uint64_t,LENGTH+1> dist1 = {};
	std::array<std::array<std::uint64_t,LENGTH>,THREAD_PER_WORK> counts1 = {};
	std::array<std::array<std::uint64_t,LENGTH+1>,THREAD_PER_WORK> dists1 = {};
	std::uint64_t noncount2 = 0;
	std::array<std::uint64_t,LENGTH> count2 = {};
	std::array<std::uint64_t,LENGTH+1> dist2 = {};
	std::array<std::uint64_t,THREAD_PER_WORK> noncounts2 = {};
	std::array<std::array<std::uint64_t,LENGTH>,THREAD_PER_WORK> counts2 = {};
	std::array<std::array<std::uint64_t,LENGTH+1>,THREAD_PER_WORK> dists2 = {};

	auto all = [&counts1, &dists1](std::size_t tid, std::bitset<LENGTH> begin, std::bitset<LENGTH> end){
		auto series = begin;
		std::bitset<divsize> half;
		do{
			series = series.to_ullong()+1;
			std::size_t qty1block = series.count();
			for(std::size_t i=0; i<divsize; ++i) half[i]=series[i];
			std::size_t pos = 0;
			std::uint64_t qty1div = half.count()<<1;
			std::uint64_t stride = qty1div>qty1block?qty1div-qty1block:qty1block-qty1div;
			while(stride!=0){
				for(std::size_t k=0, kend=qty1div>qty1block?qty1div-qty1block:qty1block-qty1div; k<kend; ++k){
					assert(pos<LENGTH);
					qty1div -= series.test(pos>>1);
					qty1div += series.test((pos+LENGTH)>>1);
					++pos;
				}
				std::uint64_t qty1divwide = qty1div + (series.test(pos>>1)+series.test((pos+LENGTH)>>1)-1)*(pos&1);
				stride = qty1divwide>qty1block?qty1divwide-qty1block:qty1block-qty1divwide;
			}
			++counts1[tid][pos];
			++dists1[tid][pos>>1];
			++dists1[tid][(pos+LENGTH+1)>>1];
		}while(series!=end);
	};

	auto pitch = [&noncounts2, &counts2, &dists2](std::size_t tid, std::bitset<LENGTH> begin, std::bitset<LENGTH> end){
		auto series = begin;
		std::bitset<divsize> half;
		do{
			series = series.to_ullong()+1;
			std::size_t qty1block = series.count();
			if(qty1block<qtyhalflow||qty1block>qtyhalfhigh){
				for(std::size_t i=0; i<divsize; ++i) half[i]=series[i];
				std::size_t pos = 0;
				std::uint64_t qty1div = half.count()<<1;
				std::uint64_t stride = qty1div>qty1block?qty1div-qty1block:qty1block-qty1div;
				while(stride>qtytolerance){
					for(std::size_t k=0, kend=qtytolerance*2+1; k<kend; ++k){
						assert(pos<LENGTH);
						qty1div -= series.test(pos>>1);
						qty1div += series.test((pos+LENGTH)>>1);
						++pos;
					}
					std::uint64_t qty1divwide = qty1div + (series.test(pos>>1)+series.test((pos+LENGTH)>>1)-1)*(pos&1);
					stride = qty1divwide>qty1block?qty1divwide-qty1block:qty1block-qty1divwide;
				}
				++counts2[tid][pos];
				++dists2[tid][pos>>1];
				++dists2[tid][(pos+LENGTH+1)>>1];
			}else{
				++noncounts2[tid];
			}
		}while(series!=end);
	};

	constexpr std::uint64_t blocksize = 1<<(LENGTH-THREAD_PER_WORK_EXP);
	constexpr std::uint64_t first = ~static_cast<std::uint64_t>(0);

	std::vector<std::thread> thread1;
	std::vector<std::thread> thread2;
	for(std::size_t i=0; i<THREAD_PER_WORK; ++i) thread1.emplace_back(all, thread1.size(), first+blocksize*i, first+blocksize*(i+1));
	for(std::size_t i=0; i<THREAD_PER_WORK; ++i) thread2.emplace_back(pitch, thread2.size(), first+blocksize*i, first+blocksize*(i+1));
	for(auto &t: thread1) t.join();
	for(auto &t: thread2) t.join();

	for(std::size_t i=0; i<LENGTH; ++i) for(auto &j: counts1) count1[i] += j[i];
	for(std::size_t i=0; i<LENGTH+1; ++i) for(auto &j: dists1) dist1[i] += j[i];
	for(auto &j: noncounts2) noncount2+=j;
	for(std::size_t i=0; i<LENGTH; ++i) for(auto &j: counts2) count2[i] += j[i];
	for(std::size_t i=0; i<LENGTH+1; ++i) for(auto &j: dists2) dist2[i] += j[i];

	for(auto i: count1) std::cout<<i<<"\t";
	std::cout<<std::endl;
	for(auto i: dist1) std::cout<<i<<"\t";
	std::cout<<std::endl;
	std::cout<<noncount2<<std::endl;
	for(auto i: count2) std::cout<<i<<"\t";
	std::cout<<std::endl;
	for(auto i: dist2) std::cout<<i<<"\t";
	std::cout<<std::endl;
}
