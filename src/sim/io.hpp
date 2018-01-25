// Version 3.0
#include "defs.hpp"
#include <iostream>

typedef typename simulation::grid<std::string> char_grid_type;

void test(bool& pass, bool compar) {
	pass = pass * compar;
}

void print_test(bool& pass, std::string name) {
	std::cout << name <<"...";
	if(pass) std::cout << "pass" << std::endl;
	else std::cout << "fail" << std::endl;
	pass=true;
}

std::ostream& operator<<(std::ostream& os, simulation::point_type xy) {
	os << xy.x <<","<< xy.y;
	return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "[";
    for(auto it=v.begin(); it!=v.end()-1; ++it) os << *it <<",";
    os << *(v.end()-1) << "]";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::map<T,T>& m) {
    for(auto& pair : m) os << pair.first<<":"<< pair.second <<" ";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<T,T>& um) {
    os << "{";
    auto it=um.begin();
    if(it!=um.end()) {
	    os << it->first <<":" << it->second;
	    ++it;
	    for(;it!=um.end(); ++it) os <<","<< it->first <<":" << it->second;
    }
    os << "}";
    return os;
}

// print visits
std::ostream& operator<<(std::ostream& os, const std::unordered_map<int, std::vector<int>>& visits) {
    os << "visits-> ";
    for(auto& pair : visits) os << pair.first<<":"<< pair.second <<" ";
    return os;
}

// print brood
std::ostream& operator<<(std::ostream& os, const std::map<int, std::unordered_map<int,int>>& brood) {
    os << "brood-> ";
    for(auto& pair : brood) os << pair.first<<":"<< pair.second <<" ";
    return os;
}


void push_land(simulation& sim)
{ //SUCCESSFUL TEST 11-08-2016
	std::ofstream file("/home/colin/projects/bombus/init/land/gaussLand.comp" + sim.pars.rid);
	for(int m=0; m<sim.pars.nrow; ++m){
		for(int n=0; n<sim.pars.ncol; ++n) {
			file << sim.land[m][n].crop << " ";
		}
		file << std::endl;
	}
}


void print_grid(char_grid_type& g) {
	for(int n=0; n!=g.cols; ++n) std::cout <<"____";
	std::cout <<"_"<< std::endl;
	for(int m=g.rows-1; m>=0; --m) {
		for(int n=0; n!=g.cols; ++n) std::cout <<"|"<< g[m][n];
		std::cout <<"|"<< std::endl;
		for(int n=0; n!=g.cols; ++n) std::cout <<"|___";
		std::cout <<"|"<< std::endl;
	}	
}

void clear_grid(char_grid_type& g) {for(auto& i : g) i = "";}

void add_plant(simulation& S, char_grid_type& land) {
	for(int row=0; row!=S.pars.nrow; ++row) {
		for(int col=0; col!=S.pars.ncol; ++col) {
			if(S.land[row][col].crop) land[row][col]+="c";
			else if(S.land[row][col].plant) land[row][col]+="w";
			else land[row][col]+=" ";
		}
	}
}

void add_flower(simulation& S, char_grid_type& land) {
	for(int row=0; row!=S.pars.nrow; ++row) {
		for(int col=0; col!=S.pars.ncol; ++col) {
			if(S.is_flowering(S.land[row][col])) land[row][col]+="*";
			else land[row][col]+=" ";
		}
	}
}

void add_space(char_grid_type& land) {for(auto& i : land) i += " ";}

void add_nests(simulation& S, char_grid_type& land) {
	for(int row=0; row!=S.pars.nrow; ++row) {
		for(int col=0; col!=S.pars.ncol; ++col) {
			auto it = S.nests.find(S.land[row][col].id);
			if(it==S.nests.end()) land[row][col]+=" ";
			else land[row][col]+=std::to_string(it->second);
		}
	}	
}

/**void add_brood(simulation& S, char_grid_type& land) {
	for(int row=0; row!=S.pars.nrow; ++row) {
		for(int col=0; col!=S.pars.ncol; ++col) {
			auto it = S.brood.find(S.land[row][col].id);
			if(it==S.brood.end()) land[row][col]+=" ";
			else land[row][col]+=std::to_string(it->second);
		}
	}	
}**/

void add_poll(simulation& S, char_grid_type& land) {
	for(int row=0; row!=S.pars.nrow; ++row) {
		for(int col=0; col!=S.pars.ncol; ++col) {
			if(S.land[row][col].poll) land[row][col]+="P";
			else land[row][col]+=" ";
		}
	}
}

void add_bank(simulation& S, char_grid_type& land) {
	for(int row=0; row!=S.pars.nrow; ++row) {
		for(int col=0; col!=S.pars.ncol; ++col) {
			if(S.land[row][col].bank > 0) land[row][col]+=std::to_string(S.land[row][col].bank);
			else land[row][col]+=" ";
		}
	}
}

void add_flowDate(simulation& S, char_grid_type& land) {
	for(int row=0; row!=S.pars.nrow; ++row) {
		for(int col=0; col!=S.pars.ncol; ++col) {
			if(S.land[row][col].flowDate >= 0) land[row][col]+=std::to_string(S.land[row][col].flowDate);
			else land[row][col]+=" ";
		}
	}
}

void print_hist(std::vector<int> hist) {
	for(int k=0; k<hist.size(); ++k) std::cout << k << ": " << std::string(hist[k],'*') << std::endl;
}
