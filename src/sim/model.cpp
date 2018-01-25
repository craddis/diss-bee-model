#include <string>
#include "stats.hpp"

struct pars_type {
  std::string rid;
  std::string version;
  int T; // Number of day in a growing season
  int Y; // Number of years in simulation
  int ncol, nrow; // Dimensions of land
  double pi; // Number of days each flower blooms
  double phi; // Foraging trips per nests per day
  double alpha, beta; // Phenology Distribution Shape Parameters
  double init_P; // Initial density of plants in wild habitat
  double init_N; // Initial density of nests on whole landscape
  double gamma; // Germination fraction
  double delta_s; // Annual seed mortality rate
  double delta_p; // Annual wild plant mortality rate
  double delta_n; // Daily nest mortality
  double lambda; // Brood survival rate
  double vartheta; // Dispersal Distribution Parameter
  double mu, nu; // Eclosion Distribution Shape Parameters
	
  pars_type(std::string rid, double delta_n, double delta_p, 
            double alpha, int phpeak, double pi, double phi,
            double init_N, double lambda, double mu, int brpeak)
      : rid(rid),
        version("3.3"),
        T(180),
        Y(100),
        ncol(256),
        nrow(256),
        pi(pi),
        phi(phi),
        alpha(alpha),
        beta(stats::mode_to_shape(alpha, phpeak, 180)),
        init_P(1),
        init_N(init_N),
        gamma(0.5),
        delta_s(0.1),
        delta_p(delta_p),
        delta_n(delta_n),
        lambda(lambda),
        vartheta(0.0625),
        mu(mu),
        nu(stats::mode_to_shape(mu, brpeak, 180)) {}
    pars_type(std::string rid)
      : rid(rid),
        version("3.3"),
        T(10),
        Y(3),
        ncol(4),
        nrow(4),
        pi(2),
        phi(1),
        alpha(1),
        beta(1),
        init_N(0.25),
        init_P(1), 
        gamma(0.5),
        delta_s(0.5),
        delta_p(0.5),
        delta_n(0.5),
        lambda(0.5),
        vartheta(0.0625),
        mu(2),
        nu(2) {}
};

#include "defs.hpp"
int simulation::get_eclosion_date() {
  return get_epoch(y+1, eclosion_dist(engine));
}
void simulation::winter() {
  nests.clear(); N=0;
}
void simulation::year() {
  for (t=0; t<pars.T; ++t) {day();} // Pollinate
  winter(); // What happens to pollinators over winter?
  for (cell_type& i : land) {
    i.poll = false;
    assign_phen(i.id);
  }
}

int main(int argc, char* argv[]) {
  pars_type pars(argv[1], std::stod(argv[2]), std::stod(argv[3]),
                 std::stod(argv[4]), std::stoi(argv[5]),
                 std::stod(argv[6]), std::stod(argv[7]),
                 std::stod(argv[8]), std::stod(argv[9]),
                 std::stod(argv[10]), std::stoi(argv[11]));
  simulation sim(pars);
  sim.init();
  sim.run();
  return 0;
};

