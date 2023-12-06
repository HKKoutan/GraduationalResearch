﻿#include <iostream>
#include <array>
#include <bitset>
#include <cstdint>
#include <cassert>

constexpr std::size_t LENGTH = 16;
constexpr double TOLERANCE = 0.125;

int main(){
	constexpr std::size_t divsize = LENGTH>>1;
	constexpr std::uint64_t qtytolerance = static_cast<std::uint64_t>(static_cast<double>(LENGTH)*TOLERANCE);
	constexpr std::uint64_t qtyhalflow = divsize-qtytolerance;
	constexpr std::uint64_t qtyhalfhigh = divsize+qtytolerance;

	std::bitset<LENGTH> series;
	series.flip();
	std::bitset<divsize> half;
	std::uint64_t noncount = 0;
	std::array<std::uint64_t,LENGTH> count = {};

	do{
		series = series.to_ullong()+1;
		std::size_t qty1block = series.count();
		if(TOLERANCE==0||qty1block<qtyhalflow||qty1block>qtyhalfhigh){
			for(std::size_t i=0; i<divsize; ++i) half[i]=series[i];
			std::size_t pos = 0;
			std::uint64_t qty1div = half.count()<<1;
			std::uint64_t stride = qty1div>qty1block?qty1div-qty1block:qty1block-qty1div;
			while(stride>qtytolerance){
				assert(pos<LENGTH);
				for(std::size_t k=0, kend=stride-qtytolerance; k<kend; ++k){
					qty1div -= series.test(pos>>1);
					qty1div += series.test((pos+LENGTH)>>1);
					++pos;
				}
				stride = qty1div>qty1block?qty1div-qty1block:qty1block-qty1div;
			}
			++count[pos];
		}else{
			++noncount;
		}
	}while(!series.all());
	std::cout<<noncount<<std::endl;
	for(auto i: count) std::cout<<i<<"\t";
	std::cout<<std::endl;
}
