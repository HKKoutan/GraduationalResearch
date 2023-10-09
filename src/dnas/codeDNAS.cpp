#include "codeDNAS.hpp"

using std::vector;
using std::uint32_t;
using std::cerr, std::endl;

code::DNAS::nucleotide_t code::DNAS::VLRLL_encode(const vector<bool> &source, vector<nucleotide_t> &code, const nucleotide_t initial_state){
	vector<bool> previous_states;
	nucleotide_t current_state = initial_state;

	code.clear();
	for(const auto i: source){
		previous_states.push_back(i);
		switch(previous_states.size()){
		case 1:
		case 3:
			break;
		
		case 2:
		case 4:
			if(previous_states[previous_states.size()-2]){
				current_state = (current_state + (previous_states[previous_states.size()-1]?0:3)) & 3;//11とandを取る
				code.push_back(current_state);
				if(!previous_states[previous_states.size()-1]) previous_states.clear();
			}else{
				current_state = (current_state + (previous_states[previous_states.size()-1]?2:1)) & 3;//11とandを取る
				code.push_back(current_state);
				previous_states.clear();
			}
			break;

		case 5:
			if(previous_states[4]){
				current_state = (current_state + 3) & 3;//11とandを取る
				code.push_back(current_state);
				previous_states.clear();
			}
			break;
		
		case 6:
			current_state = (current_state + (previous_states[previous_states.size()-1]?2:1)) & 3;//11とandを取る
			code.push_back(current_state);
			previous_states.clear();
			break;
		
		default:
			cerr<<"error"<<endl;
			exit(9);
		}
	}
	//端数処理(0パディング)
	if(previous_states.size()%2 != 0){
		current_state = (current_state + (previous_states[previous_states.size()-1]?3:1)) & 3;//11とandを取る
		code.push_back(current_state);
	}

	return current_state;
}

vector<bool>::iterator code::DNAS::VLRLL_encode(const vector<bool>::iterator &source_begin, const vector<bool>::iterator &source_end, const vector<nucleotide_t>::iterator &code_begin, const vector<nucleotide_t>::iterator &code_end, const nucleotide_t initial_state){
	vector<bool> unprocessed;
	nucleotide_t current_state = initial_state;
	vector<bool>::iterator source_it = source_begin;
	vector<nucleotide_t>::iterator code_it = code_begin;

	while(code_it!=code_end){
		if(source_it!=source_end){
			unprocessed.push_back(*source_it);
			source_it++;
		}else{
			unprocessed.push_back(0);
		}
		switch(unprocessed.size()){
		case 1:
		case 3:
			break;
		
		case 2:
		case 4:
			if(unprocessed[unprocessed.size()-2]){
				current_state = (current_state + (unprocessed[unprocessed.size()-1]?0:3)) & 3;//11とandを取る
				*code_it = current_state;
				code_it++;
				if(!unprocessed[unprocessed.size()-1]) unprocessed.clear();
			}else{
				current_state = (current_state + (unprocessed[unprocessed.size()-1]?2:1)) & 3;//11とandを取る
				*code_it = current_state;
				code_it++;
				unprocessed.clear();
			}
			break;

		case 5:
			if(unprocessed[4]){
				current_state = (current_state + 3) & 3;//11とandを取る
				*code_it = current_state;
				code_it++;
				unprocessed.clear();
			}
			break;
		
		case 6:
			current_state = (current_state + (unprocessed[unprocessed.size()-1]?2:1)) & 3;//11とandを取る
			*code_it = current_state;
			code_it++;
			unprocessed.clear();
			break;
		
		default:
			cerr<<"error"<<endl;
			exit(9);
		}
	}
	return source_it;
}

code::DNAS::nucleotide_t code::DNAS::modified_VLRLL_encode(const vector<bool> &source, vector<nucleotide_t> &code, const nucleotide_t initial_state){
	vector<bool> unprocessed;
	nucleotide_t current_state = initial_state;

	code.clear();
	for(const auto i: source){
		unprocessed.push_back(i);
		switch(unprocessed.size()){
		case 0:
		case 1:
		case 3:
		case 5:
			break;
		
		case 2:
		case 4:
			if(unprocessed[unprocessed.size()-2]){
				current_state = (current_state + (unprocessed[unprocessed.size()-1]?0:3)) & 3;//11とandを取る
				code.push_back(current_state);
				if(!unprocessed[unprocessed.size()-1]) unprocessed.clear();
			}else{
				current_state = (current_state + (unprocessed[unprocessed.size()-1]?2:1)) & 3;//11とandを取る
				code.push_back(current_state);
				unprocessed.clear();
			}
			break;

		case 6:
			if(unprocessed[unprocessed.size()-2]){
				current_state = (current_state + 3) & 3;//11とandを取る
				code.push_back(current_state);
				unprocessed.clear();
			}else{
				current_state = (current_state + (unprocessed[unprocessed.size()-1]?2:1)) & 3;//11とandを取る
				code.push_back(current_state);
				unprocessed.clear();
			}
			break;
		
		default:
			cerr<<"error"<<endl;
			exit(9);
		}
	}
	if(unprocessed.size()%2 != 0){
		cerr<<"error"<<endl;
		exit(9);
	}

	return current_state;
}

void code::DNAS::interim_map(const vector<nucleotide_t> &source, vector<bool> &code){
	code.clear();
	for(const auto &i: source){
		code.push_back((i&2)>>1);
		code.push_back(i&1);
	}
}

void code::DNAS::interim_demap(const vector<bool> &source, vector<nucleotide_t> &code){
	if(source.size()%2){
		cerr<<"invalid source length."<<endl;
		exit(11);
	}
	code.clear();
	for(uint32_t i=0; i<source.size(); i+=2){
		code.push_back((source[i]<<1)+source[i+1]);
	}
}

void code::DNAS::VLRLL_decode(const vector<nucleotide_t> &source, vector<bool> &decode, const nucleotide_t initial_state){
	nucleotide_t previous = initial_state;
	uint32_t zeros = 0;

	decode.clear();
	for(const auto &i: source){
		switch((i-previous)&3){
		case 0:
			decode.push_back(1);
			decode.push_back(1);
			zeros++;
			break;
		
		case 1:
			decode.push_back(0);
			decode.push_back(0);
			zeros=0;
			break;
		
		case 2:
			decode.push_back(0);
			decode.push_back(1);
			zeros=0;
			break;

		case 3:
			decode.push_back(1);
			if(zeros<2) decode.push_back(0);
			zeros=0;
			break;

		default:
			cerr<<"error"<<endl;
			exit(9);
		}
		previous = i;
	}
}

void code::DNAS::modified_VLRLL_decode(const vector<nucleotide_t> &source, vector<bool> &decode, nucleotide_t initial_state){
	nucleotide_t previous = initial_state;

	decode.clear();
	for(const auto &i: source){
		switch((i-previous)&3){
		case 0:
			decode.push_back(1);
			decode.push_back(1);
			break;
		
		case 1:
			decode.push_back(0);
			decode.push_back(0);
			break;
		
		case 2:
			decode.push_back(0);
			decode.push_back(1);
			break;

		case 3:
			decode.push_back(1);
			decode.push_back(0);
			break;

		default:
			cerr<<"error"<<endl;
			exit(9);
		}
		previous = i;
	}
}
void code::DNAS::nt_addequalizing_encode(const vector<nucleotide_t> &cr, vector<nucleotide_t> &crbar, vector<bool> &info, uint32_t qty_AT, uint32_t qty_GC){
	nucleotide_t diff = 0;
	info.resize(cr.size());
	crbar.resize(cr.size());

	for(uint64_t i=0; i<cr.size(); i++){
		nucleotide_t prev_state = (cr[i]+diff)&3;
		if(qty_AT>qty_GC && prev_state>>1==0){
			crbar[i] = (~prev_state)&3;
			info[i] = 1;
		}else if(qty_AT<qty_GC && prev_state>>1==1){
			crbar[i] = (~prev_state)&3;
			info[i] = 1;
		}else{
			crbar[i] = prev_state;
			info[i] = 0;
		}
		diff = crbar[i]-cr[i];
		qty_AT += !(crbar[i]>>1);
		qty_GC += (crbar[i]>>1);
	}
}

void code::DNAS::nt_addequalizing_decode(const vector<nucleotide_t> &crbar, const vector<bool> &info, vector<nucleotide_t> &cr){
	nucleotide_t diff = 0;
	cr.resize(crbar.size());

	for(uint64_t i=0; i<crbar.size(); i++){
		if(info[i]) cr[i] = ((~crbar[i])-diff)&3;
		else cr[i] = (crbar[i]-diff)&3;
		diff = crbar[i]-cr[i];
	}
}

void code::DNAS::nt_qty_count(const vector<nucleotide_t> &c, uint32_t &qty_AT, uint32_t &qty_GC){
	qty_AT=0;
	qty_GC=0;
	for(const auto &i: c){
		qty_AT += i>>1;
		qty_GC += !(i>>1);
	}
}
