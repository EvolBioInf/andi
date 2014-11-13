/**
 * This program can create genome sequences with a specific distance.
 */

#include <iostream>
#include <random>
#include <functional>
#include <string>
#include <getopt.h>

using namespace std;

void usage();
void print_seq( unsigned, unsigned, int, int, double);

int main(int argc, char *argv[]){

	random_device rd{};
	auto seed = rd();
	int length = 1000;
	int line_length = 70;

	auto seqs = vector<double>({0,0.1});
	int num_seqs = 1;

	int check;
	while((check = getopt(argc, argv, "l:L:d:")) != -1){
		switch(check) {
			case 'l': length = stoi(optarg); break;
			case 'L': line_length = stoi(optarg); break;
			case 'd': seqs[num_seqs++] = stod(optarg); break;
			case '?':
			default: usage(); return 1;
		}
	}

	if( num_seqs < 2){
		seqs[num_seqs++] = 0.1;
	}

	for( int i=0; i< num_seqs; i++){
		cout << ">S" << i << endl;
		print_seq( seed, rd(), length, line_length, seqs[i]);
	}

	return 0;
}


static auto ACGT = "ACGT";
static auto NO_A = "CGT";
static auto NO_C = "AGT";
static auto NO_G = "ACT";
static auto NO_T = "ACG";

void print_seq( unsigned base_seed, unsigned mut_seed, int length, int line_length, double divergence){
	char line[line_length+1];
	line[line_length] = '\0';

	auto base_rand = default_random_engine{base_seed};
	auto base_dist = uniform_int_distribution<int>{0,3};
	auto base_acgt = [&]{return ACGT[base_dist(base_rand)];};

	auto mut_rand = default_random_engine{mut_seed};
	auto mut_dist = uniform_real_distribution<double>{0,1};
	auto mut = bind( mut_dist, mut_rand);
	auto mut_acgt = uniform_int_distribution<int>{0,2};
	auto mutate = [&](char c){
		int idx = mut_acgt(mut_rand);
		switch(c){
			case 'A': return NO_A[idx];
			case 'C': return NO_C[idx];
			case 'G': return NO_G[idx];
			case 'T': return NO_T[idx];
			default: return 'X';
		}
	};

	double nucleotides = (double)length;
	double mutations = nucleotides * divergence;

	for(int i= length, j; i > 0; i -= j){
		j = min(line_length, i);

		for(auto k=0; k<j; k++){
			char c = base_acgt();

			if( mut() < mutations / nucleotides ){
				c = mutate(c);
				mutations--;
			}

			line[k] = c;
			nucleotides--;
		}

		line[j] = '\0';
		cout << line << endl;
	}
}

void usage(){
	const static char *str = {
		"test_rand [-l length] [-d dist]\n"
	};
	cerr << str;
}
