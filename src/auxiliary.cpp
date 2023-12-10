#include <iostream>
#include <array>
#include <bitset>
#include <cstdint>
#include <cassert>

constexpr std::size_t LENGTH = 16;
constexpr double TOLERANCE = 0;
constexpr bool PITCH = false;

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
	std::array<std::uint64_t,LENGTH> dist = {};

	do{
		series = series.to_ullong()+1;
		std::size_t qty1block = series.count();
		if(TOLERANCE==0||qty1block<qtyhalflow||qty1block>qtyhalfhigh){
			for(std::size_t i=0; i<divsize; ++i) half[i]=series[i];
			std::size_t pos = 0;
			std::uint64_t qty1div = half.count()<<1;
			std::uint64_t stride = qty1div>qty1block?qty1div-qty1block:qty1block-qty1div;
			while(stride>qtytolerance){
				if constexpr(PITCH){
					for(std::size_t k=0, kend=qtytolerance*2+1; k<kend; ++k){
						assert(pos<LENGTH);
						qty1div -= series.test(pos>>1);
						qty1div += series.test((pos+LENGTH)>>1);
						++pos;
					}
				}else{
					for(std::size_t k=0, kend=1; k<kend; ++k){
						assert(pos<LENGTH);
						qty1div -= series.test(pos>>1);
						qty1div += series.test((pos+LENGTH)>>1);
						++pos;
					}
				}
				std::uint64_t qty1divwide = qty1div + (series.test(pos>>1)+series.test((pos+LENGTH)>>1)-1)*(pos&1);
				stride = qty1divwide>qty1block?qty1divwide-qty1block:qty1block-qty1divwide;
			}
			++count[pos];
			for(std::size_t j=pos>>1, jend=(pos+LENGTH+1)>>1; j<jend; ++j) ++dist[j];
		}else{
			++noncount;
		}
	}while(!series.all());
	std::cout<<noncount<<std::endl;
	for(auto i: count) std::cout<<i<<"\t";
	std::cout<<std::endl;
	for(auto i: dist) std::cout<<i<<"\t";
	std::cout<<std::endl;

	
}
