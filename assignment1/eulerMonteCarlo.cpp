#include <iostream> // cout 
#include <iomanip> // For formatting output
#include <vector> 
#include <cmath>
#include <iostream>
#include <matplot/matplot.h> // ploting
#include <boost/math/distributions/normal.hpp> // distributions
// for random generation
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace std;


double identity(double x) {
    return x;
}

double positivePart(double x) {
    return std::max(0.0, x);
}

double fabs(double x) {
    return fabs(x);
}

class HestonSimulator {
    public:
        // Parameters
        double interestRate, strike, initialStockPrice, initialVariance;
        double longTermVariance, corr, rateOfReversion, volOfVol, time;
        int steps;
    
        // Output
        std::unique_ptr<double[]> logPrices;
        std::unique_ptr<double[]> variances;
        std::unique_ptr<double[]> intermediateVars;
        std::unique_ptr<double[]> times;
        double executionTime;
    
        // Constructor
        HestonSimulator(double time, 
                        double interestRate, 
                        double strike, 
                        double initialStockPrice, 
                        double initialVariance, 
                        double longTermVariance, 
                        double corr, 
                        double rateOfReversion, 
                        double volOfVol, 
                        int steps,
                        unsigned int seed1 = 42,
                        unsigned int seed2 = 1337) : 
                        
                        rng1(seed1),
                        rng2(seed2),
                        nd(0.0, 1.0),
                        gen1(rng1, nd),
                        gen2(rng2, nd) 
        {
            this->interestRate = interestRate;
            this->time = time;
            this->strike = strike;
            this->initialStockPrice = initialStockPrice;
            this->longTermVariance = longTermVariance;
            this->corr = corr;
            this->rateOfReversion = rateOfReversion;
            this->volOfVol = volOfVol;
            this->initialVariance=initialVariance;
            this->steps = steps;
            // "Auto" deallocate off heap with smart pointers
            logPrices = std::make_unique<double[]>(steps);
            variances = std::make_unique<double[]>(steps);
            intermediateVars = std::make_unique<double[]>(steps);
            times = std::make_unique<double[]>(steps);
        }
    
        // Main simulation method
        void run(std::function<double(double)> f1,std::function<double(double)> f2,std::function<double(double)> f3) {
            logPrices[0] = log(initialStockPrice);
            variances[0] = initialVariance;
            intermediateVars[0] = initialVariance;
            times[0]=0;
        
            double delta = time / (double)steps;
            for (int i = 1; i<steps; i++) {
                // calc variance first
                double deltaW1 = sqrt(delta)*gen1();
                double deltaW2 = corr*deltaW1 + sqrt(1-pow(corr,2)) * sqrt(delta)*gen2();
                intermediateVars[i] = f1(intermediateVars[i-1]) - rateOfReversion*delta*(f2(intermediateVars[i-1])-longTermVariance) + volOfVol*sqrt(f3(intermediateVars[i-1]))*deltaW1;
                variances[i] = f3(intermediateVars[i]);
                logPrices[i] = logPrices[i-1] + (interestRate-0.5*variances[i])*delta + sqrt(variances[i])*deltaW2;
                times[i] = times[i-1]+delta;
            }
        }
    private:
        // random number generations -> alternatively we apply box muller transform to achieve same result
        boost::mt19937 rng1;
        boost::mt19937 rng2;
        boost::normal_distribution<> nd;
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > gen1;
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > gen2;
    };

// TODO 
// create a generic function 

int main() {
    double interestRate = 0.05;
    double strike = 100;
    double initialStockPrice = 100; 
    double initialVariance = 0.09;
    double longTermVariance = 0.09;
    double corr = -0.3;
    double rateOfReversion = 2; 
    double volOfVol = 1;
    double time = 5;
    int steps = 80;
    HestonSimulator absorptionSim = HestonSimulator(time,interestRate,strike,initialStockPrice,initialVariance,longTermVariance,corr,rateOfReversion,volOfVol,steps);
    int numPaths = 300;
    std::vector<std::vector<double>> allPaths(numPaths);
    std::vector<double> timeVec(steps);
    for (int i = 0; i < steps; i++) {
        timeVec[i] = i * (time / (double)steps);
    }

    for (int i = 0; i < numPaths; i++) {
        absorptionSim.run(positivePart, positivePart, positivePart);
        allPaths[i] = std::vector<double>(absorptionSim.logPrices.get(), absorptionSim.logPrices.get() + steps);
    }

    matplot::figure()->size(1920,1080);
    matplot::hold(true);  // keep all lines on the same figure
    for (const auto& path : allPaths) {
        matplot::plot(timeVec, path);
    };
    matplot::title("Monte Carlo sim Pricing");
    matplot::xlabel("M-time steps");
    matplot::ylabel("log option prices");
    matplot::show();
    matplot::save("plot.png");
    return -1;
}