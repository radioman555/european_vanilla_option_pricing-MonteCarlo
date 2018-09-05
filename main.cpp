#include <algorithm>    // Needed for the "max" function
#include <cmath>
#include <iostream>

using namespace std;

// A simple implementation of the Box-Muller algorithm, used to generate
// gaussian random numbers - necessary for the Monte Carlo method below
// Note that C++11 actually provides std::normal_distribution<> in
// the <random> library, which can be used instead of this function
double gaussian_box_muller() {
    double x = 0.0;
    double y = 0.0;
    double euclid_sq = 0.0;

    // Continue generating two uniform random variables
    // until the square of their "euclidean distance"
    // is less than unity
    do {
        x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
        y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
        euclid_sq = x*x + y*y;
    } while (euclid_sq >= 1.0);

    return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}

// Pricing a European vanilla call option with a Monte Carlo method
double monte_carlo_call_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T) {
    double S_adjust = S * exp(T*(r-0.5*v*v));
    double S_cur = 0.0;
    double payoff_sum = 0.0;

    for (int i=0; i<num_sims; i++) {
        double gauss_bm = gaussian_box_muller();
        S_cur = S_adjust * exp(sqrt(v*v*T)*gauss_bm);
        payoff_sum += std::max(S_cur - K, 0.0);
    }

    return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}

// Pricing a European vanilla put option with a Monte Carlo method
double monte_carlo_put_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T) {
    double S_adjust = S * exp(T*(r-0.5*v*v));
    double S_cur = 0.0;
    double payoff_sum = 0.0;

    for (int i=0; i<num_sims; i++) {
        double gauss_bm = gaussian_box_muller();
        S_cur = S_adjust * exp(sqrt(v*v*T)*gauss_bm);
        payoff_sum += std::max(K - S_cur, 0.0);
    }

    return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}

int main(int argc, char **argv) {
    // First we create the parameter list
    int num_sims = 10000000;   // Number of simulated asset paths

    // First we create the parameter list
    double S = 0;  // Option price
    double K = 0;  // Strike price
    double r = 0;   // Risk-free rate (5%)
    double v = 0;    // Volatility of the underlying (20%)
    double T = 0;    // One year until expiry


    // Then we ask the user for the vaules of the parameters

    std::cout << "Number of Paths: " << num_sims << std::endl;
    cout << "Option price?:    " <<  endl;
    cin >> S ;
    cout << "Strike price?:    " << endl;
    cin >> K ;
    cout << "Risk-Free Rate?:  " <<  endl;
    cin >> r ;
    cout << "Volatility?:     " <<  endl;
    cin >> v;
    cout << "Time to maturity?:   " << endl;
    cin >> T ;


// Then we calculate the call/put values via Monte Carlo
    double call = monte_carlo_call_price(num_sims, S, K, r, v, T);
    double put = monte_carlo_put_price(num_sims, S, K, r, v, T);

    // Finally we output the parameters and prices
    cout << "Underlying:      " << S << endl;
    cout << "Strike:          " << K << endl;
    cout << "Risk-Free Rate:  " << r << endl;
    cout << "Volatility:      " << v << endl;
    cout << "Maturity:        " << T << endl;

    cout << "Call Price:      " << call << endl;
    cout << "Put Price:       " << put << endl;

    return 0;
}