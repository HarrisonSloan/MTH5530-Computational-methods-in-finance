#include <iostream> // cout 
#include <iomanip> // For formatting output
#include <vector> 
#include <cmath>
#include <iostream>
#include <matplot/matplot.h> // ploting
#include <boost/math/distributions/normal.hpp> // distributions

using namespace std;


float blackScholes(float sigma, float time, float interestRate, float strike, float initialStockPrice, int putOrCall) {
    // Here we simply tree time as T-t = Time
    boost::math::normal_distribution<> norm(0.0, 1.0);  // Standard normal

    float d1 = (std::log(initialStockPrice / strike) + (interestRate + 0.5f * sigma * sigma) * time) / (sigma * std::sqrt(time));
    float d2 = d1 - sigma * std::sqrt(time);

    float discountedStrike = strike * std::exp(-interestRate * time);

    float callPrice = initialStockPrice * boost::math::cdf(norm, d1) - discountedStrike * boost::math::cdf(norm, d2);

    if (putOrCall == 1) {
        // Call option
        return callPrice;
    } else {
        // Put option
        float putPrice = discountedStrike * boost::math::cdf(norm, -d2) - initialStockPrice * boost::math::cdf(norm, -d1);
        return putPrice;
    }
}

float blackScholesVega(float sigma, float time, float interestRate, float strike, float initialStockPrice) {
    // Note that putcall option is not needed as its the same for either type of option

    boost::math::normal_distribution<> norm(0.0, 1.0);  // Standard normal
    float d1 = (std::log(initialStockPrice / strike) + (interestRate + 0.5f * sigma * sigma) * time) / (sigma * std::sqrt(time));
    // need standard normal
    return boost::math::pdf(norm, d1) * std::sqrt(time);
}

float bisection(float marketPrice, float time, float interestRate, float strike, float initialStockPrice, int putOrCall, float tolerance) {

    // Improvements 
    // Dont need to substrct the market price off everything here, just add to the tolerance

    // volatility ranges
    float low = 0.001;
    float high = 1;
    float mid = (low+high) / 2;

    // BlackScholes 
    float diffLow = blackScholes(low, time, interestRate, strike, initialStockPrice, putOrCall)-marketPrice;
    float diffHi = blackScholes(high, time, interestRate, strike, initialStockPrice, putOrCall)-marketPrice;
    float diffMid = blackScholes(mid, time, interestRate, strike, initialStockPrice, putOrCall)-marketPrice;

    bool signLow = diffLow < 0;
    bool signHigh = diffHi < 0;
    bool signMid = diffMid < 0;

    if (signLow == signHigh) throw("Error root doesnt exist in interval, check data");

    // essentially binary search until you get tolerance.
    while (fabs(diffMid)>=tolerance) {
        if (signLow != signMid) {
            high = mid;
            diffHi = blackScholes(high, time, interestRate, strike, initialStockPrice, putOrCall)-marketPrice;
        } else if (signMid != signHigh) {
            low = mid;
            diffLow = blackScholes(low, time, interestRate, strike, initialStockPrice, putOrCall)-marketPrice;
        } else {
            throw("Root doesnt exist in converging interval, check data");
        }
        mid = (low+high) / 2;

        diffMid = blackScholes(mid, time, interestRate, strike, initialStockPrice, putOrCall)-marketPrice;
    
        signLow = diffLow < 0;
        signHigh = diffHi < 0;
        signMid = diffMid < 0;
    }

    return mid;
}

float newtonRaphson(float marketPrice, float time, float interestRate, float strike, float initialStockPrice, int putOrCall, float tolerance, float intialGuess) {
    const int maxIterations = 100;
    // intial guess of sigma reasonable amount is say 20 % ??? -> This could be better researched, likely better to have this be passed depending on the data we are working with
    for (int i = 0; i < maxIterations; i++) {
        float res = blackScholes(intialGuess, time, interestRate, strike, initialStockPrice, putOrCall)-marketPrice;
        // If close to evaluating to zero return the sigma
        if (fabs(res) < tolerance) return intialGuess;
        // revise sigma by finding tangent and its x intercept
        float dres = blackScholesVega(intialGuess, time, interestRate, strike, initialStockPrice);
        float dx = res/dres;
        intialGuess -= dx;
        if (intialGuess < 0) throw("sigma exceeded bounds (negative volatility)");
    }
    throw("Max iterations exceeded, did not converge to given tolerance");
}

float newtonBisectionHybrid(float marketPrice, float time, float interestRate, float strike, float initialStockPrice, int putOrCall, float tolerance, float intialGuess){
    // inspired by numerical recipes
}

int main() {

    return -1;
}