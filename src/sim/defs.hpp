// Version 3.x
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <fstream>

#define PI 3.14159265

struct simulation
{
	struct point_type
	{
		double x, y;
		point_type(double x, double y) : x(x), y(y) {}
	};

	template<class T>
	struct grid
	{
		typedef std::vector<T> storage_type;
		typedef typename storage_type::value_type value_type;
		typedef typename storage_type::iterator iterator;
		typedef typename storage_type::size_type size_type;

		grid(size_type nr, size_type nc) : rows(nr), cols(nc), storage(nr * nc) {}

		iterator operator[](size_type r) {return storage.begin() + r * cols;}

		value_type& operator[](const point_type& p) {return storage[int(p.y)*cols + int(p.x)];}

		iterator begin() {return storage.begin();}

		iterator end() {return storage.end();}

		size_type size() {return storage.size();}

		size_type rows, cols;
		storage_type storage;
	};

	struct cell_type
	{
		int id = -1;
		point_type center = point_type(-1,-1); // Center of cell
		bool crop = false; // Is this cell designated as cropland?
		bool plant = false; // Does this cell have a plant on it
		bool poll = false; // Has this plant been pollinated?
		int flowDate = -1; // First day of flowering [0...T)
		int pi = 0; // Number of flowering days for this plant
		int bank = 0; // seed bank
	};

	// Creates new point at random angle (uniform) and distance (negative exponential)
	point_type project(point_type xy, double r, double theta)
	{
		xy.x += r * std::cos(theta);
		xy.y += r * std::sin(theta);
		return xy;
	}

	// Enforces torus: points outside grid are wrapped around back into the grid
	point_type wrap(point_type& xy) 
	{
		if(xy.x < 0) {while(xy.x < 0) xy.x += pars.ncol;}
		else if(xy.x > pars.ncol) {while(xy.x > pars.ncol) xy.x -= pars.ncol;}
		if(xy.y < 0) {while(xy.y < 0) xy.y += pars.nrow;}
		else if(xy.y > pars.nrow) {while(xy.y > pars.nrow) xy.y -= pars.nrow;}
	}	

	int rdisp(int i_orig)
	{
		point_type dest_xy = project(land.storage[i_orig].center, disp_dist(engine), unif_dist(engine)*2*PI);
		wrap(dest_xy);
		return land[dest_xy].id;
	}

	int pick_one(std::vector<int>& V)
	{
		std::uniform_int_distribution<int> lottery(0, V.size()-1);
		return *std::next(V.begin(), lottery(engine)); // Choose the "first visit" randomly	
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// NEST DYNAMICS ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////// GENERATE VISITS : LOOP NESTS /////////////////////////////////////////
	bool is_late_bloomer(const cell_type& i) {return i.flowDate > pars.T-i.pi;} // TESTED 1-12-17

	bool is_flowering(const cell_type& i)
	{ // TESTED 1-12-17
		if(is_late_bloomer(i)) return i.plant  &&  (t >= i.flowDate || t < i.pi-(pars.T-i.flowDate));
		else return i.plant  &&  i.flowDate + i.pi > t  &&  t >= i.flowDate;
	}

	bool is_full(const cell_type& i) {return !i.poll;} // TESTED 1-12-17

	void visit(int id_orig, int id_dest)
	{ // TESTED 1-12-17
		cell_type& i_dest = land.storage[id_dest];
		// If plant is flowering and unpollinated...
		if(is_flowering(i_dest) && is_full(i_dest)) {
			// Add the nest to the lottery for that cell's nectar
			visits[id_dest].push_back(id_orig);
		} // Ignores visits with no potential for nectar collection or pollination		
	}

	void gen_visits() // Generate and record visits
	{ // TESTED 1-12-17
		for(auto& pair : nests) { // pair : {cell id, number of nests}
			for(int b = 0; b < pair.second; ++b) { // Each nest launches 1 foraging bout
				visit(pair.first, rdisp(pair.first)); // Disperse from nest to random cell
			}
		}
	}
	/////////////////////// PICK WINNERS : LOOP VISITS ///////////////////////////////////////////
	int get_epoch(int year, int day) {return year*pars.T + day;} // TESTED 1-12-17
	int get_year(int epoch) {return epoch/pars.T;} // TESTED 1-12-17
	int get_day(int epoch) {return epoch%pars.T;} // TESTED 1-12-17

	int get_eclosion_date();

	void pollinate(int id_nest, int id_visit)
	{ // TESTED 1-12-17
		assert(is_flowering(land.storage[id_visit]) && is_full(land.storage[id_visit]));
		land.storage[id_visit].poll = true; // Pollinate flower
		// Increment pollination event counter
		if(land.storage[id_visit].crop) {++Cc;} else {++Cw;}
		++brood[get_eclosion_date()][id_nest]; // Mass provision brood
	} 

	void pick_winners() // Pick winners to collect resources for home cell
	{
		// pair : {visited cell id, vector of visiting nest origins}
		for(auto& pair : visits) pollinate(pick_one(pair.second), pair.first);
		visits.clear();
	}
	//////////////////// NEST MORTALITY : LOOP NESTS ///////////////////////////////////////////
	int get_ndead(int n)
	{
		std::binomial_distribution<int> mort_dist(n, pars.delta_n);
		return mort_dist(engine);
	}

	bool kill_nests(int& n, int ndead)
	{ // TESTED 1-12-17
		assert(ndead <= n);
		n -= ndead; // Decrement nests in cell i
		N -= ndead; // Decrement landscape nest counter
		// if no nests left in cell i, tag cell i for removal from "nests" list
		return n == 0;
	}

	void nest_death() 
	{ // TESTED 1-12-17
		std::vector<int> to_be_erased;
		for(auto& pair : nests) {// pair :{cell id, number of nests}
			if(kill_nests(pair.second, get_ndead(pair.second))) to_be_erased.push_back(pair.first);
		}
		for(auto& i : to_be_erased) nests.erase(i); // Remove cells that no longer have nests
	}
	//////////////////// NEST DISPERSAL : LOOP BROOD[EPOCH] //////////////////////////////////////////
	bool brood_ready() {auto it=brood.begin(); return it != brood.end() && it->first == get_epoch(y,t);} // TESTED 1-12-17

	int get_nqueen(int b) 
	{
		std::binomial_distribution<int> eclosion(b, pars.lambda);
		return eclosion(engine);
	}	

	bool is_habitat(int id) {return !land.storage[id].crop;} // TESTED 1-12-17
	bool is_habitat(const cell_type& i) {return !i.crop;}

	void disperse_queen_to(int dest_id)
	{ // TESTED 1-12-17
		++nests[dest_id]; // Disperse there and create new nest
		++N; // Add to nest counter
	}	

	void nest_dispersal()
	{ // TESTED 1-13-17
		assert((brood.empty()) ? true : brood.begin()->first >= get_epoch(y,t));
		// If next brood is schedule to eclose today...
		if(brood_ready()) {
			for(auto& pair : brood.begin()->second) { // pair :{cell id, number of brood}
				for(int q = get_nqueen(pair.second); q > 0; --q) { // For each survivor...
					int dest_id = rdisp(pair.first); // Pick Random destination 
					// If you can nest in crops OR if destination is wild habitat...
					if(is_habitat(dest_id)) disperse_queen_to(dest_id);
				}
			}
			brood.erase(get_epoch(y,t));
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// PLANT DYNAMICS ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void assign_phen(int id)
	{
		cell_type& i = land.storage[id];
		if(i.plant) {
			i.pi = stats::random_round(pars.pi, engine);
			if(i.crop) i.flowDate = crop_phen_dist(engine);
			else i.flowDate = wild_phen_dist(engine);
		}
	}

	void plant_death()
	{
		for(cell_type& i : land) {
			if(is_habitat(i) && i.plant && plant_dies(engine)) { // only wild plants can die
				i.plant = false;
				i.poll = false;
				i.flowDate = -1;
				i.pi = 0;
			}
		}
	}
	
	void seed_dispersal()
	{
		for(cell_type& i_orig : land) {
			if(is_habitat(i_orig) && i_orig.poll) { // only pollinated wild plants produce seeds
				int id_dest = rdisp(i_orig.id);
				// If seed lands in habitat...
				if(is_habitat(id_dest)) ++land.storage[id_dest].bank; // Add to seed bank
			}
		}
	}

	void seed_death()
	{ // SUCCESSFUL VISUAL TEST 2016-09-09
		for(cell_type& i : land) {
			if(i.bank > 0) {
				std::binomial_distribution<int> mort_dist(i.bank, pars.delta_s);
				i.bank -= mort_dist(engine);
			}
		}
	}

	void recruitment(int id)
	{
		cell_type& i = land.storage[id];
		if(i.bank > 0) {
			std::binomial_distribution<int> germ_dist(i.bank, pars.gamma);
			int ngerm = germ_dist(engine); // Draw number of germinating seeds
			assert(is_habitat(i));
			// If 1+ seeds germinate in empty habitat...
			if(!i.plant && ngerm>0) {
				assert(i.flowDate==-1 && i.pi==0);
				i.plant = true; // Recruit 1 seedling to adult plant
				i.flowDate = wild_phen_dist(engine); // Assign random flowering date
				i.pi = stats::random_round(pars.pi, engine); // Assign discrete flowering period
			}
			i.bank -= ngerm; // Remove germinated seeds from seed bank
			assert(i.bank >= 0);
		}
	}

	void kill_plant(int id)
	{
		cell_type& i = land.storage[id];
		i.plant = false;
		i.poll = false;
		i.flowDate = -1;
		i.pi = 0;
	}

	void disperse_seed(int id_orig, int id_dest)
	{
		// If seed lands in habitat...
		if(is_habitat(id_dest)) ++land.storage[id_dest].bank; // Add to seed bank
	}

	void prune_bank(int id)
	{
		cell_type& i = land.storage[id];
		if(i.bank > 0) {
			std::binomial_distribution<int> mort_dist(i.bank, pars.delta_s);
			i.bank -= mort_dist(engine);
		}		
	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// INITIALIZE //////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void import_land()
	{ // SUCCESSFUL TEST 11-08-2016
		std::ifstream file("/home/colin/projects/bombus/init/land/gaussLand." + pars.rid);
		std::string line;
		int row=0;
		while(std::getline(file, line)) {
			std::istringstream iss(line);
			int col=0;
			for(std::string i; iss >> i; ++col) {
				land[row][col].crop = bool(stoi(i));
			}
			++row;
			assert(col==pars.ncol);
		}
		assert(row==pars.nrow);
	}

	void init_nests()
	{
		std::vector<int> habitat_ids;
		for(auto& i : land) {if(is_habitat(i)) habitat_ids.push_back(i.id);}
		int NI = pars.init_N * pars.nrow * pars.ncol; // Convert initial density to absolute number of nests
		for(int n=0; n<NI; ++n) {++brood[get_eclosion_date()][pick_one(habitat_ids)];}
	}
	
	void init()
	{
		import_land(); // Sets crop
		for(cell_type& i : land) i.plant = true; // Sets plant
		for(cell_type& i : land) assign_phen(i.id); // Sets pi, flowDate
		init_nests(); // Builds brood
	}

	simulation
	(	// parameters
		const pars_type& p
	) :	// initializations
		t(0), y(-1), N(0),
		over_I(1.0 / p.nrow / p.ncol),
		Rw(0), Cw(0), Fw(0),
		Rc(0), Cc(0), Fc(0),
		pars(p),
		engine(std::random_device{}()),
		unif_dist(0,1), // [0,1)
		land(p.nrow, p.ncol),
		plant_dies(p.delta_p),
		seed_germinates(p.gamma),
		disp_dist(p.vartheta),
		crop_phen_dist(p.alpha, p.beta, p.T),
		eclosion_dist(p.mu, p.nu, p.T),
		wild_phen_dist(0, p.T)
	{	// BEGIN CONSTRUCTOR BODY
		// assign cell IDs
		int ID = 0;
		for(cell_type& i : land) i.id = ID++;
		// fill cell center values
		for(int m=0; m!=pars.nrow; ++m){
			for(int n=0; n!=pars.ncol; ++n) {
				land[m][n].center = point_type(n+0.5, m+0.5);
			}
		}
	}

	void reset_CFR() {Cc=0; Fc=0; Rc=0; Cw=0; Fw=0; Rw=0;}

	void count_FR() 
	{	
		for(cell_type& i : land) {
			if(is_flowering(i)) {
				if(i.crop) {
					++Fc;
					if(is_full(i)) ++Rc;
				}
				else {
					++Fw;
					if(is_full(i)) ++Rw;
				}
			}
		}
	}

	void print_output()
	{
		std::cout << pars.rid<<" "<< y <<" "<< t <<" "<< N <<" "<< Cc <<" "<< Fc <<" "<< Rc <<" "<< Cw <<" "<< Fw <<" "<< Rw << std::endl;
	}

	void day()
	{
		reset_CFR();
		count_FR();
		gen_visits();
		pick_winners();
		nest_death();
		nest_dispersal();
		print_output();
	}

	bool runaway_growth()
	{
		return N > (pars.ncol*pars.nrow);
	}

	bool is_extinct()
	{
		return nests.empty() && brood.empty();
	}

	void winter();

	void year();

	void run()
	{
		y = 0;
		while(y < pars.Y && !is_extinct() && !runaway_growth()) {
			year(); // Year 0->Y-1
			++y; // Year 1->Y
		}
	}

	pars_type pars;
	int t, y;
	double over_I;
	grid<cell_type> land;
	std::unordered_map<int,int> nests;
	std::unordered_map<int,std::vector<int>> visits;
	std::map<int,std::unordered_map<int,int>> brood;
	std::default_random_engine engine;
	std::bernoulli_distribution plant_dies;
	std::bernoulli_distribution seed_germinates;
	std::uniform_real_distribution<double> unif_dist;
	std::exponential_distribution<double> disp_dist;
	stats::discrete_beta_distribution crop_phen_dist;
	stats::discrete_beta_distribution eclosion_dist;
	std::uniform_int_distribution<int> wild_phen_dist;
	int N, Fc, Rc, Cc, Fw, Rw, Cw;
};
