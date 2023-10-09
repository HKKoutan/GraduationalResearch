#include "IDECbase.hpp"

using std::vector;
using std::int32_t, std::uint64_t, std::size_t;
using std::cout, std::cerr, std::flush, std::endl;
using code::IDEC::MarkerCode, code::IDEC::Marker_Encoding, code::IDEC::BCJR_IDEC;

////////////////////////////////////////////////////////////////
//                                                            //
//                           Marker                           //
//                                                            //
////////////////////////////////////////////////////////////////

MarkerCode::MarkerCode(const vector<bool> &marker, size_t interval, size_t sourcelength):
	marker(marker),
	interval(interval),
	sourcelength(sourcelength),
	codelength(sourcelength+(marker.size()*sourcelength/interval)-(sourcelength%interval==0?marker.size():0))
{}

////////////////////////////////////////////////////////////////
//                                                            //
//                      Marker_Encoding                       //
//                                                            //
////////////////////////////////////////////////////////////////

double Marker_Encoding::LVR_to_Pr0(double LVR){
	// constexpr auto LVR_BOUND = 11.0903548889591249;
	constexpr auto LVR_BOUND = 17.328679484196309;
	constexpr auto P_MAX = 0.99999997019767761;
	constexpr auto P_MIN = 2.98023223876953125e-8;
	// if(std::fabs(LVR)>=LVR_BOUND) LVR = (LVR<0.0?-LVR_BOUND:LVR_BOUND);
	if(std::fabs(LVR)>=LVR_BOUND) return LVR<0.0?P_MIN:P_MAX;
	// auto e = std::exp(LVR);
	// return e/(1+e);
	return 1.0/(1.0+std::exp(-LVR));
}

vector<double> Marker_Encoding::Pr0_init() const{
	vector Pr0(marker->codesize(), 0.0);
	auto pri = Pr0.begin(), prend = Pr0.end();
	while(pri!=prend){
		for(uint64_t i=0; i<marker->intervalsize() && pri!=prend; ++i) *pri++ = 0.5;
		for(auto mi=marker->begin(), mend=marker->end(); mi!=mend&&pri!=prend; ++mi) *pri++ = !*mi;
	}
	return Pr0;
}

vector<double> Marker_Encoding::Pr0_init(const vector<double> &LEVR) const{
	vector Pr0(marker->codesize(), 0.0);
	auto pri = Pr0.begin(), prend = Pr0.end();
	auto lei = LEVR.cbegin();
	while(pri!=prend){
		for(uint64_t i=0; i<marker->intervalsize() && pri!=prend; ++i) *pri++ = LVR_to_Pr0(*lei++);
		for(auto mi=marker->begin(), mend=marker->end(); mi!=mend&&pri!=prend; ++mi) *pri++ = !*mi;
	}
	return Pr0;
}

vector<bool> Marker_Encoding::insert(const vector<bool> &source) const{
	vector<bool> code;
	auto si = source.begin();
	const auto send = source.end();
	const auto interval = marker->intervalsize();
	while(si!=send){
		for(size_t i=0; i<interval && si!=send; ++i) code.push_back(*si++);
		if(si!=send) for(const auto &m: *marker) code.push_back(m);
	}
	return code;
}

vector<bool> Marker_Encoding::extract(const vector<bool> &code) const{
	vector<bool> source(marker->sourcesize());
	auto ci = code.cbegin(), cend = code.cend();
	auto si = source.begin();
	const auto interval = marker->intervalsize(), markersize = marker->size();
	while(ci!=cend){
		for(size_t i=0; i<interval && ci!=cend; ++i) *si++ = *ci++;
		for(size_t i=0; i<markersize && ci!=cend; i++) ++ci;
	}
	return source;
}

vector<double> Marker_Encoding::extract(const vector<double> &code) const{
	vector<double> source(marker->sourcesize());
	auto ci = code.cbegin(), cend = code.cend();
	auto si = source.begin();
	const auto interval = marker->intervalsize(), markersize = marker->size();
	while(ci!=cend){
		for(size_t i=0; i<interval && ci!=cend; ++i) *si++ = *ci++;
		for(size_t i=0; i<markersize && ci!=cend; i++) ++ci;
	}
	return source;
}

////////////////////////////////////////////////////////////////
//                                                            //
//                         BCJR_IDEC                          //
//                                                            //
////////////////////////////////////////////////////////////////

BCJR_IDEC::BCJR_IDEC(size_t length, const int32_t drift_ceil, const int32_t drift_floor):
	length(length),
	drift_ceil(std::abs(drift_ceil)),
	drift_floor(std::abs(drift_floor)),
	smax(this->drift_ceil+this->drift_floor+1),
	alpha(length+1, vector<double>(smax+2)),//第二軸を両端1ずつパディング
	beta(length+1, vector<double>(smax+2))
{}

void BCJR_IDEC::decode_init(const vector<bool> &y){
	for(auto ai: alpha) for(auto aij: ai) aij = 0.0;
	for(auto bi: beta) for(auto bij: bi) bij = 0.0;
	alpha.front()[drift_floor+1] = 1.0;
	beta.back()[drift_floor+1+y.size()-length] = 1.0;
}

//symbol_prob(y, x, i, s(i), s(i-1)): 時点iで状態s(i-1)から状態s(i)に遷移し、そのときの送信シンボルがxだったとしたときに受信系列がyになる確率

void BCJR_IDEC::iterate(vector<double> &LLR, const vector<bool> &y, const vector<double> &Pr0, const double rate_ins, const double rate_del, const double rate_sub){
	// constexpr auto LVR_BOUND = 17.328679484196309;
	// constexpr auto LVR_BOUND = 11.0903548889591249;
	const double rate_trans = 1.0-rate_ins-rate_del;

	//forward
	{
		size_t i=0u;
		for(size_t iend=length-1; i<iend; ++i){
			size_t s=(i<=drift_floor?drift_floor+1-i:1u);
				alpha[i+1][s]=
					(s==1?2.0*rate_ins:rate_ins)*((y[i+s-1-drift_floor]||y[i+s-drift_floor]?0.0:Pr0[i])+(y[i+s-1-drift_floor]&&y[i+s-drift_floor]?(1.0-Pr0[i]):0.0))*alpha[i][s-1]+
					rate_trans*(y[i+s-1-drift_floor]?Pr0[i]*rate_sub+(1.0-Pr0[i])*(1.0-rate_sub):Pr0[i]*(1.0-rate_sub)+(1.0-Pr0[i])*rate_sub)*alpha[i][s];
				++s;
			for(size_t send=(drift_ceil+i<y.size()? smax: drift_floor+y.size()-i); s<send; ++s){
				alpha[i+1][s]=
					rate_ins*((y[i+s-1-drift_floor]||y[i+s-drift_floor]?0.0:Pr0[i])+(y[i+s-1-drift_floor]&&y[i+s-drift_floor]?(1.0-Pr0[i]):0.0))*alpha[i][s-1]+
					rate_del*alpha[i][s+1]+
					rate_trans*(y[i+s-1-drift_floor]?Pr0[i]*rate_sub+(1.0-Pr0[i])*(1.0-rate_sub):Pr0[i]*(1.0-rate_sub)+(1.0-Pr0[i])*rate_sub)*alpha[i][s];
			}
			alpha[i+1][s]=
				(s==smax?2.0*rate_del:rate_del)*alpha[i][s+1]+
				(y[i+s-1-drift_floor]?Pr0[i]*rate_sub+(1.0-Pr0[i])*(1.0-rate_sub):Pr0[i]*(1.0-rate_sub)+(1.0-Pr0[i])*rate_sub)*rate_trans*alpha[i][s];
		}
		if(y.size()>length-drift_floor){
			size_t s=(i<=drift_floor?drift_floor+1-i:1u);
				alpha[i+1][s]=
					(s==1?2.0*rate_ins:rate_ins)*((y[i+s-1-drift_floor]||y[i+s-drift_floor]?0.0:Pr0[i])+(y[i+s-1-drift_floor]&&y[i+s-drift_floor]?(1.0-Pr0[i]):0.0))*alpha[i][s-1]+
					rate_trans*(y[i+s-1-drift_floor]?Pr0[i]*rate_sub+(1.0-Pr0[i])*(1.0-rate_sub):Pr0[i]*(1.0-rate_sub)+(1.0-Pr0[i])*rate_sub)*alpha[i][s];
				++s;
			for(size_t send=(drift_ceil+i<y.size()?smax:drift_floor+y.size()-i); s<send; ++s){
				alpha[i+1][s]=
					rate_ins*((y[i+s-1-drift_floor]||y[i+s-drift_floor]?0.0:Pr0[i])+(y[i+s-1-drift_floor]&&y[i+s-drift_floor]?(1.0-Pr0[i]):0.0))*alpha[i][s-1]+
					rate_del*alpha[i][s+1]+
					rate_trans*(y[i+s-1-drift_floor]?Pr0[i]*rate_sub+(1.0-Pr0[i])*(1.0-rate_sub):Pr0[i]*(1.0-rate_sub)+(1.0-Pr0[i])*rate_sub)*alpha[i][s];
			}
			alpha[i+1][s]=
				(s==smax?2.0*rate_del:rate_del)*alpha[i][s+1]+
				rate_trans*(y[i+s-1-drift_floor]?Pr0[i]*rate_sub+(1.0-Pr0[i])*(1.0-rate_sub):Pr0[i]*(1.0-rate_sub)+(1.0-Pr0[i])*rate_sub)*alpha[i][s];
		}else{
			size_t s=1u;
			alpha[i+1][s]=
				rate_del*alpha[i][s+1]+
				rate_trans*(y[i+s-1-drift_floor]?Pr0[i]*rate_sub+(1.0-Pr0[i])*(1.0-rate_sub):Pr0[i]*(1.0-rate_sub)+(1.0-Pr0[i])*rate_sub)*alpha[i][s];
		}
	}

	//backward
	{
		size_t i=length;
		if(y.size()>length-drift_floor){
			size_t s=(i<=drift_floor?drift_floor+1-i:1);
			beta[i-1][s]=
				(s==1?2.0*rate_ins:rate_ins)*((y[i+s-2-drift_floor]||y[i+s-1-drift_floor]?0.0:Pr0[i-1])+(y[i+s-2-drift_floor]&&y[i+s-1-drift_floor]?(1.0-Pr0[i-1]):0.0))*beta[i][s+1]+
				rate_trans*(y[i+s-1-drift_floor]?Pr0[i-1]*rate_sub+(1.0-Pr0[i-1])*(1.0-rate_sub):Pr0[i-1]*(1.0-rate_sub)+(1.0-Pr0[i-1])*rate_sub)*beta[i][s];
			++s;
			for(size_t send=(drift_ceil+i<y.size()?smax:drift_floor+y.size()-i); s<send; ++s){
				beta[i-1][s]=
					rate_ins*((y[i+s-2-drift_floor]||y[i+s-1-drift_floor]?0.0:Pr0[i-1])+(y[i+s-2-drift_floor]&&y[i+s-1-drift_floor]?(1.0-Pr0[i-1]):0.0))*beta[i][s+1]+
					rate_trans*(y[i+s-1-drift_floor]?Pr0[i-1]*rate_sub+(1.0-Pr0[i-1])*(1.0-rate_sub):Pr0[i-1]*(1.0-rate_sub)+(1.0-Pr0[i-1])*rate_sub)*beta[i][s]+
					rate_del*beta[i][s-1];
			}
			beta[i-1][s]=
				(s==smax?2.0*rate_del:rate_del)*beta[i][s-1]+
				(y[i+s-1-drift_floor]?Pr0[i-1]*rate_sub+(1.0-Pr0[i-1])*(1.0-rate_sub):Pr0[i-1]*(1.0-rate_sub)+(1.0-Pr0[i-1])*rate_sub)*rate_trans*beta[i][s];
		}else{
			size_t s=1u;
			beta[i-1][s]=
				rate_del*beta[i][s-1]+
				(y[i+s-1-drift_floor]?Pr0[i-1]*rate_sub+(1.0-Pr0[i-1])*(1.0-rate_sub):Pr0[i-1]*(1.0-rate_sub)+(1.0-Pr0[i-1])*rate_sub)*rate_trans*beta[i][s];
		}
		for(--i; i>0; --i){
			size_t s=(i<=drift_floor?drift_floor+1-i:1);
			beta[i-1][s]=
				(s==1?2.0*rate_ins:rate_ins)*((y[i+s-2-drift_floor]||y[i+s-1-drift_floor]?0.0:Pr0[i-1])+(y[i+s-2-drift_floor]&&y[i+s-1-drift_floor]?(1.0-Pr0[i-1]):0.0))*beta[i][s+1]+
				rate_trans*(y[i+s-1-drift_floor]?Pr0[i-1]*rate_sub+(1.0-Pr0[i-1])*(1.0-rate_sub):Pr0[i-1]*(1.0-rate_sub)+(1.0-Pr0[i-1])*rate_sub)*beta[i][s];
			++s;
			for(size_t send=(drift_ceil+i<y.size()?smax:drift_floor+y.size()-i); s<send; ++s){
				beta[i-1][s]=
					rate_ins*((y[i+s-2-drift_floor]||y[i+s-1-drift_floor]?0.0:Pr0[i-1])+(y[i+s-2-drift_floor]&&y[i+s-1-drift_floor]?(1.0-Pr0[i-1]):0.0))*beta[i][s+1]+
					rate_trans*(y[i+s-1-drift_floor]?Pr0[i-1]*rate_sub+(1.0-Pr0[i-1])*(1.0-rate_sub):Pr0[i-1]*(1.0-rate_sub)+(1.0-Pr0[i-1])*rate_sub)*beta[i][s]+
					rate_del*beta[i][s-1];
			}
			beta[i-1][s]=
				(s==smax?2.0*rate_del:rate_del)*beta[i][s-1]+
				(y[i+s-1-drift_floor]?Pr0[i-1]*rate_sub+(1.0-Pr0[i-1])*(1.0-rate_sub):Pr0[i-1]*(1.0-rate_sub)+(1.0-Pr0[i-1])*rate_sub)*rate_trans*beta[i][s];
		}
	}
	// if(drift_ceil+2<y.size()) for(size_t i=y.size()-drift_ceil-2; i>0; --i){
	// 	for(size_t s=(i<=drift_floor?drift_floor+2-i:1), send=drift_floor+drift_ceil+2; s<send; ++s){
	// 		beta[i-1][s]=
	// 			((y[i+s-2-drift_floor]||y[i+s-1-drift_floor]?0.0:Pr0[i-1])+(y[i+s-2-drift_floor]&&y[i+s-1-drift_floor]?(1.0-Pr0[i-1]):0.0))*rate_ins*beta[i][s+1]+
	// 			(y[i+s-2-drift_floor]?Pr0[i-1]*rate_sub+(1.0-Pr0[i-1])*(1.0-rate_sub):Pr0[i-1]*(1.0-rate_sub)+(1.0-Pr0[i-1])*rate_sub)*rate_trans*beta[i][s]+
	// 			rate_del*beta[i][s-1];
	// 	}
	// }

	//integrate
	for(size_t i=0; i<=drift_floor; ++i){
		double sum0 = 0.0;
		double sum1 = 0.0;
		{
			size_t s=drift_floor+1-i;
			sum0 += alpha[i][s]*(
				(y[i+s-1-drift_floor]?rate_sub:1.0-rate_sub)*rate_trans*beta[i+1][s]+
				rate_del*beta[i+1][s+1]
			);
			sum1 += alpha[i][s]*(
				(y[i+s-1-drift_floor]?1.0-rate_sub:rate_sub)*rate_trans*beta[i+1][s]+
				rate_del*beta[i+1][s+1]
			);
		}
		for(size_t s=drift_floor+2-i, send=(drift_ceil+i<y.size()? drift_floor+drift_ceil+2: drift_floor+y.size()+1-i); s<send; ++s){
			sum0 += alpha[i][s]*(
				(y[i+s-2-drift_floor]||y[i+s-1-drift_floor]?0.0:1.0)*rate_ins*beta[i+1][s-1]+
				(y[i+s-1-drift_floor]?rate_sub:1.0-rate_sub)*rate_trans*beta[i+1][s]+
				rate_del*beta[i+1][s+1]
			);
			sum1 += alpha[i][s]*(
				(y[i+s-2-drift_floor]&&y[i+s-1-drift_floor]?1.0:0.0)*rate_ins*beta[i+1][s-1]+
				(y[i+s-1-drift_floor]?1.0-rate_sub:rate_sub)*rate_trans*beta[i+1][s]+
				rate_del*beta[i+1][s+1]
			);
		}
		// auto llr = std::log(sum0*Pr0[i])-std::log(sum1*(1.0-Pr0[i]));
		// if(std::isnan(llr)) LLR[i] = 0.0;
		// else if(std::fabs(llr)>=LVR_BOUND) LLR[i] = llr<0.0?-LVR_BOUND:LVR_BOUND;
		// else LLR[i] = llr;
		// LLR[i] = std::log(sum0*Pr0[i])-std::log(sum1*(1.0-Pr0[i]));
		LLR[i] = std::log(sum0)-std::log(sum1);
	}
	for(size_t i=drift_floor+1u, iend=length; i<iend; ++i){
		double sum0 = 0.0;
		double sum1 = 0.0;
		for(size_t s=1u, send=(drift_ceil+i<y.size()? drift_floor+drift_ceil+2: drift_floor+y.size()+1-i); s<send; ++s){
			// sum0 +=
			// 	((y[i+s-2-drift_floor]||y[i+s-1-drift_floor]?0.0:1.0)*rate_ins*alpha[i-1][s-1]+
			// 	(y[i+s-1-drift_floor]?rate_sub:(1.0-rate_sub))*rate_trans*alpha[i-1][s]+
			// 	rate_del*alpha[i-1][s+1])*beta[i][s];
			// sum1 +=
			// 	((y[i+s-2-drift_floor]&&y[i+s-1-drift_floor]?1.0:0.0)*rate_ins*alpha[i-1][s-1]+
			// 	(y[i+s-1-drift_floor]?rate_sub:(1.0-rate_sub))*rate_trans*alpha[i-1][s]+
			// 	rate_del*alpha[i-1][s+1])*beta[i][s];
			sum0 += alpha[i][s]*(
				(y[i+s-2-drift_floor]||y[i+s-1-drift_floor]?0.0:1.0)*rate_ins*beta[i+1][s-1]+
				(y[i+s-1-drift_floor]?rate_sub:1.0-rate_sub)*rate_trans*beta[i+1][s]+
				rate_del*beta[i+1][s+1]
			);
			sum1 += alpha[i][s]*(
				(y[i+s-2-drift_floor]&&y[i+s-1-drift_floor]?1.0:0.0)*rate_ins*beta[i+1][s-1]+
				(y[i+s-1-drift_floor]?1.0-rate_sub:rate_sub)*rate_trans*beta[i+1][s]+
				rate_del*beta[i+1][s+1]
			);
		}
		// auto llr = std::log(sum0)-std::log(sum1);
		// if(std::isnan(llr)) LLR[i] = 0.0;
		// else if(std::fabs(llr)>=LVR_BOUND) LLR[i] = llr<0.0?-LVR_BOUND:LVR_BOUND;
		// else LLR[i] = llr;
		// LLR[i] = std::log(sum0*Pr0[i])-std::log(sum1*(1.0-Pr0[i]));
		LLR[i] = std::log(sum0)-std::log(sum1);
	}
}

vector<bool> BCJR_IDEC::estimate(const vector<double> &LLR){//Pr0=0.5のときのみ使うこと
	vector<bool> estimation(LLR.size());
	for(size_t i=0u, iend=LLR.size(); i<iend; ++i) estimation[i] = LLR[i]<0.0;
	// auto ei = estimation.begin();
	// for(const auto &li: LLR) *ei++ = li<0;
	return estimation;
}
