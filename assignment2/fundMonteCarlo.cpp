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
#include <boost/random.hpp>
using namespace std;


float fundSimulation(float sigma, float years, float feeRate, float riskFreeRate, float startPrice, int paths, bool print) {

    // V_{i+1} = V_i * exp((r-f-0.5sigma^2)delta + sigma*sqrt(delta)*standardRandomNormal)

    // Singular array to store all of the path values
    // We could use multithreading here to boost performance
    std::unique_ptr<float[]> fundValue = std::make_unique<float[]>((int)(years*paths+paths));
    
    // Random number generation
    boost::random::mt19937 rng;
    boost::normal_distribution<> nd(0.0,1.0);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > generator(rng,nd);

    double delta = 1 / (double)years;
    float payOffSum = 0;
    float totalFeeSum = 0;
    for (int path = 0; path<paths ; path++) {
        int startIndex = path*years+path;
        fundValue[path] = startPrice;
        float totalFee = 0;
        float maxYearFundVal = 0;
        fundValue[startIndex]=startPrice;
        for (int year = 1; year<=years ; year++) {
            fundValue[startIndex+year] = fundValue[startIndex+year-1]*exp((riskFreeRate-feeRate-0.5*pow(sigma,2))*delta+sigma*sqrt(delta)*generator());
            // Add the fee charged and see if a new max has been reached
            totalFee += exp(-year*riskFreeRate*delta) * feeRate * delta * fundValue[startIndex+year];
            maxYearFundVal = max(maxYearFundVal, fundValue[startIndex+year]);

        }
        // Sum the total payoffs and fees to be later averaged
        payOffSum += exp(-riskFreeRate*years)*(maxYearFundVal-fundValue[startIndex+years]);
        totalFeeSum += totalFee;
    }
    float difference = totalFeeSum/paths - payOffSum/paths;
    if (print == true) {
        cout << "Average fee across paths " << totalFeeSum/paths << "\n";
        cout << "Average payoff in present value " << payOffSum/paths << "\n";
    }
    return difference;
}



int main() {

    float sigma = 0.2;
    float years = 5;
    float interestRate = 0.05;
    float startPrice = 1000;
    int paths = 1000000;
    float tolerance = 0.5;

    cout << "Parameters:\n"; 
    cout << "years: " << years << "\n";
    cout << "sigma: " << sigma*100 << "%\n";
    cout << "risk free iterest rate: " << interestRate*100 << "%\n";
    cout << "Start fund value: " << startPrice << "\n";
    cout << "Number of paths per simulation: " << paths << "\n";
    cout << "Tolerance: " << tolerance<< "\n";
    // Fee rates intial search bounds
    float low = 0.001;
    float high = 1;
    float mid = (low+high) / 2;


    float diffLow = fundSimulation(sigma,years,low,interestRate,startPrice,paths,false);
    float diffHi = fundSimulation(sigma,years,high,interestRate,startPrice,paths,false);
    float diffMid =fundSimulation(sigma,years,mid,interestRate,startPrice,paths,false);

    bool signLow = diffLow < 0;
    bool signHigh = diffHi < 0;
    bool signMid = diffMid < 0; 

    if (signLow == signHigh) throw("Error root doesnt exist in interval, check data");

    int maxIterations = 50;
    cout << "\nBegin Bisection\n\n";
    while (fabs(diffMid)>=tolerance || maxIterations == 0) {
        if (signLow != signMid) {
            high = mid;
            diffHi = fundSimulation(sigma,years,high,interestRate,startPrice,paths,true);
        } else if (signMid != signHigh) {
            low = mid;
            diffLow = fundSimulation(sigma,years,low,interestRate,startPrice,paths,true);
        } else {
            throw("Root doesnt exist in converging interval, check data");
        }
        cout << "Fee rate is between " << low << " and " << high << "\n";
        mid = (low+high) / 2;

        diffMid = fundSimulation(sigma,years,mid,interestRate,startPrice,paths,false);
    
        signLow = diffLow < 0;
        signHigh = diffHi < 0;
        signMid = diffMid < 0;
        maxIterations -= 1;
    }
    cout << "\nIterations used " << 50-maxIterations;
    cout << "\nFinal fee rate: " << mid*100 << "%\n";

    return 0;
}