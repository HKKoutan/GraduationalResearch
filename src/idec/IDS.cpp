#include "IDS.hpp"

using std::vector;
using std::int32_t, std::int64_t;
using std::cerr, std::endl;
using channel::IDS;

IDS::IDS(double rate_ins, double rate_del, int32_t drift_max, double rate_sub, int64_t seed):
	mt(seed),
	uniform(0.0, 1.0),
	rate_ins(rate_ins),
	rate_indel(rate_ins+rate_del),
	rate_sub(rate_sub),
	drift_max(drift_max)
{
	if(rate_indel>1.0 || rate_sub>1.0){
		cerr<<"IDS: invalid error rate"<<endl;
		exit(12);
	}
}

vector<bool> IDS::noise(const vector<bool> &in){
	vector<bool> out;
	if(drift_max!=0){
		int32_t drift = 0;
		for(auto ii: in){
			double rand = uniform(mt);
			if(rand<rate_ins){
				if(drift<drift_max){//insert = duplicate
					out.push_back(ii);
					out.push_back(ii);
					++drift;
				}else{
					--drift;
				}
			}else if(rand<rate_indel){
				if(-drift<drift_max){
					--drift;
				}else{
					out.push_back(ii);
					out.push_back(ii);
					++drift;
				}
			}else{
				out.push_back((uniform(mt)<rate_sub)?!ii:ii);
			}
		}
	}else{
		for(auto ii: in){
			double rand = uniform(mt);
			if(rand<rate_ins){
				out.push_back(ii);
				out.push_back(ii);
			}else if(rand<rate_indel){//delete: outに何も入れない
			}else{
				out.push_back((uniform(mt)<rate_sub)?!ii:ii);
			}
		}
	}
	return out;
}
