// #include "sequencer.hpp"

// using std::vector, std::array;
// using std::int64_t;
// using std::cerr, std::endl;
// using code::DNAS::nucleotide_t;
// using channel::Nanopore_Sequencing;

// Nanopore_Sequencing::Nanopore_Sequencing(double alpha) : Nanopore_Sequencing(alpha, 0){}

// Nanopore_Sequencing::Nanopore_Sequencing(double alpha, int64_t seed)
// 	: mt(seed), uniform(0,1), error_rate{4*alpha, alpha, 0.01, 0}, non_error_rate{1-error_rate[1]-error_rate[2]-error_rate[3], 1-error_rate[0]-error_rate[1]-error_rate[2], 1-error_rate[1]-error_rate[2]-error_rate[3], 1-error_rate[0]-error_rate[1]-error_rate[2]},
// 	condprob{
// 	non_error_rate[0],error_rate[1],error_rate[3],error_rate[2],
// 	error_rate[1],non_error_rate[1],error_rate[2],error_rate[0],
// 	error_rate[3],error_rate[2],non_error_rate[2],error_rate[1],
// 	error_rate[2],error_rate[0],error_rate[1],non_error_rate[3]
// 	}{
// 	if(alpha>0.198||alpha<0){
// 		cerr<<"Nanopore_Sequencing: invalid alpha"<<endl;
// 		exit(10);
// 	}
// }

// void Nanopore_Sequencing::noise(vector<nucleotide_t> &c){
// 	for(auto &i: c){
// 		//std::cout<<int(i)<<endl;
// 		nucleotide_t j=0;
// 		double rand = uniform(mt)-condprob[i][j];
// 		while(rand>=condprob[i][j]){
// 			j+=1;
// 			rand-=condprob[i][j];
// 		}
// 		i=j;
// 	}
// }

// void Nanopore_Sequencing::message_LLR(const vector<nucleotide_t> &cm, vector<double> &LLRm) const{
// 	LLRm.clear();
// 	for(const auto &ci: cm){
// 		switch(ci){
// 		case 0:
// 			LLRm.push_back(log((non_error_rate[0]+error_rate[1])/(error_rate[3]+error_rate[2])));
// 			LLRm.push_back(log((non_error_rate[0]+error_rate[3])/(error_rate[1]+error_rate[2])));
// 			break;
		
// 		case 1:
// 			LLRm.push_back(log((non_error_rate[1]+error_rate[1])/(error_rate[0]+error_rate[2])));
// 			LLRm.push_back(-log((non_error_rate[1]+error_rate[0])/(error_rate[1]+error_rate[2])));
// 			break;

// 		case 2:
// 			LLRm.push_back(-log((non_error_rate[2]+error_rate[1])/(error_rate[3]+error_rate[2])));
// 			LLRm.push_back(log((non_error_rate[2]+error_rate[3])/(error_rate[1]+error_rate[2])));
// 			break;
		
// 		case 3:
// 			LLRm.push_back(-log((non_error_rate[3]+error_rate[1])/(error_rate[0]+error_rate[2])));
// 			LLRm.push_back(-log((non_error_rate[3]+error_rate[0])/(error_rate[1]+error_rate[2])));
// 			break;
		
// 		default:
// 			cerr<<"Nanopore_Sequencing: invalid cm "<<int(ci)<<endl;
// 			exit(10);
// 		}
// 	}
// }

// void Nanopore_Sequencing::redundancy_LLR(const vector<nucleotide_t> &cm, vector<double> &LLRr, nucleotide_t initial_state) const{
// 	constexpr double p3_to_11 = 1.0/22.0;
// 	constexpr double p3_to_10 = 21.0/22.0;

// 	nucleotide_t previous = initial_state;

// 	LLRr.clear();

// 	for(const auto &current: cm){
// /*	condprob{
// 	non_error_rate[0],error_rate[1],error_rate[3],error_rate[2],
// 	error_rate[1],non_error_rate[1],error_rate[2],error_rate[0],
// 	error_rate[3],error_rate[2],non_error_rate[2],error_rate[1],
// 	error_rate[2],error_rate[0],error_rate[1],non_error_rate[3]}
// */
// 		double P0X =//遷移語が1 or 2になる組み合わせ
// 			condprob[previous][previous] * (condprob[(previous+1)&3][current] + condprob[(previous+2)&3][current]) +
// 			condprob[(previous-1)&3][previous] * (condprob[previous][current] + condprob[(previous+1)&3][current]) +
// 			condprob[(previous-2)&3][previous] * (condprob[(previous-1)&3][current] + condprob[previous][current]) +
// 			condprob[(previous-3)&3][previous] * (condprob[(previous-2)&3][current] + condprob[(previous-1)&3][current]);
// 		double P1X =//遷移語が0 or 3になる組み合わせ
// 			condprob[previous][previous] * (condprob[previous][current] + condprob[(previous+3)&3][current]) +
// 			condprob[(previous-1)&3][previous] * (condprob[(previous-1)&3][current] + condprob[(previous+2)&3][current]) +
// 			condprob[(previous-2)&3][previous] * (condprob[(previous-2)&3][current] + condprob[(previous+1)&3][current]) +
// 			condprob[(previous-3)&3][previous] * (condprob[(previous-3)&3][current] + condprob[previous][current]);
// 		double PX0 =//遷移語が1 or 3(*p3_to_10)になる組み合わせ
// 			(condprob[(current-1)&3][previous] + p3_to_10*condprob[(current-3)&3][previous]) * condprob[current][current] +
// 			(condprob[current][previous] + p3_to_10*condprob[(current-2)&3][previous]) * condprob[(current+1)&3][current] +
// 			(condprob[(current+1)&3][previous] + p3_to_10*condprob[(current-1)&3][previous]) * condprob[(current+2)&3][current] +
// 			(condprob[(current+2)&3][previous] + p3_to_10*condprob[current][previous]) * condprob[(current+3)&3][current];
// 		double PX1 =//遷移語が0 or 2 or 3(*p3_to_11)になる組み合わせ
// 			(condprob[current][previous] + condprob[(current-2)&3][previous] + p3_to_11*condprob[(current-3)&3][previous]) * condprob[current][current] +
// 			(condprob[(current+1)&3][previous] + condprob[(current-1)&3][previous] + p3_to_11*condprob[(current-2)&3][previous]) * condprob[(current+1)&3][current] +
// 			(condprob[(current+2)&3][previous] + condprob[current][previous] + p3_to_11*condprob[(current-1)&3][previous]) * condprob[(current+2)&3][current] +
// 			(condprob[(current+3)&3][previous] + condprob[(current+1)&3][previous] + p3_to_11*condprob[current][previous]) * condprob[(current+3)&3][current];

// 		previous = current;

// 		LLRr.push_back(log(P0X/P1X));
// 		LLRr.push_back(log(PX0/PX1));
// 	}
// }
