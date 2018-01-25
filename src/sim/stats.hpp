#include <boost/math/distributions/beta.hpp>
#include <random>
#include <algorithm>

namespace stats {

int random_round(double x, std::default_random_engine& engine) 
{ // TESTED 1-17-17
	int base = std::floor(x);
	std::bernoulli_distribution roundup(x-base);
	return base + int(roundup(engine));
}

struct discrete_beta_distribution
{ // TESTED 1-17-17
	boost::math::beta_distribution<double> beta_dist;
	std::uniform_real_distribution<double> unif_dist;
	int n;
	
	discrete_beta_distribution(double a, double b, int n) : beta_dist(a,b), n(n) {}

	int operator() (std::default_random_engine& engine)
	{
		return std::trunc(n*boost::math::quantile(beta_dist, unif_dist(engine)));
	}
};

double mode_to_shape(double alpha, double mode, int T) {return T/mode*(alpha-1)-alpha+2;}

double mean(std::vector<int>::iterator begin, std::vector<int>::iterator end)
{
	double sum = std::accumulate(begin, end, 0);
	return sum / std::distance(begin, end);
}

double rmse(std::vector<int>::iterator begin, std::vector<int>::iterator end)
{
	double x_bar = mean(begin, end);
	double mse = std::accumulate(begin, end, 0, [&](double sum, double x){return sum + std::pow(x-x_bar, 2);}) / std::distance(begin, end);
	return std::sqrt(mse);
}

}

