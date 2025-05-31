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
    return initialStockPrice * boost::math::pdf(norm, d1) * std::sqrt(time);
}

float blackScholesDelta(float sigma, float time, float interestRate, float strike, float initialStockPrice, int putOrCall) {

    boost::math::normal_distribution<> norm(0.0, 1.0);  // Standard normal
    float d1 = (std::log(initialStockPrice / strike) + (interestRate + 0.5f * sigma * sigma) * time) / (sigma * std::sqrt(time));
    // need standard normal
    if (putOrCall == 1) {
        // Call option
        return boost::math::cdf(norm, d1);
    } else {
        // Put option
        return boost::math::cdf(norm, d1) - 1;
    }
}

float blackScholesGamma(float sigma, float time, float interestRate, float strike, float initialStockPrice) {
    // Note that putcall option is not needed as its the same for either type of option

    boost::math::normal_distribution<> norm(0.0, 1.0);  // Standard normal
    float d1 = (std::log(initialStockPrice / strike) + (interestRate + 0.5f * sigma * sigma) * time) / (sigma * std::sqrt(time));
    // need standard normal
    return boost::math::pdf(norm, d1) / (initialStockPrice* sigma * std::sqrt(time));
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

int main() {
    // Create table 


    vector<vector<double>> volTable;
    vector<double> tenors = {1.0/12.0, 2.0/12.0, 3.0/12.0, 6.0/12.0, 9.0/12.0, 1.0, 18.0/12.0, 2.0, 3.0, 4.0, 5.0};
    vector<double> moneynesses = {.8,.9,.95,.975,1.0,1.025,1.05,1.1,1.2};
    // First row
    volTable.push_back({22.29,17.84,14.40,12.90,11.72,11.02,10.88,12.27,15.53});
    volTable.push_back({21.72,16.68,13.97,12.88,12.00,11.30,10.96,11.20,12.34});
    volTable.push_back({22.89,16.09,13.73,12.85,12.12,11.52,11.25,11.51,12.49});
    volTable.push_back({19.76,15.81,14.26,13.59,12.99,12.48,12.12,11.78,11.83});
    volTable.push_back({19.29,16.03,14.71,14.12,13.59,13.15,12.81,12.42,12.33});
    volTable.push_back({18.53,16.00,15.03,14.57,14.13,13.72,13.36,12.84,12.47});
    volTable.push_back({18.20,16.25,15.42,15.02,14.65,14.30,13.99,13.47,12.94});
    volTable.push_back({18.10,16.41,15.68,15.33,15.00,14.70,14.41,13.93,13.33});
    volTable.push_back({18.13,16.76,16.16,15.88,15.62,15.37,15.15,14.75,14.17});
    volTable.push_back({17.96,16.76,16.23,15.99,15.75,15.54,15.33,14.97,14.40});
    volTable.push_back({17.81,16.67,16.17,15.94,15.72,15.51,15.31,14.96,14.40});

    // Task 1 
    vector<vector<double>> blackSholesPrices;
    double intialVal = 1.0;
    double interestRate = 0.04;
    for (int i = 0; i<volTable.size(); i++) {
        double tenor = tenors[i];
        vector<double> temp;
        for (int j = 0; j<volTable[i].size();j++) {

            double impliedVol = volTable[i][j] / 100.00;
            double moneyness = moneynesses[j];
            double strike = intialVal*moneyness;
            temp.push_back(blackScholes(impliedVol,tenor,interestRate,strike,intialVal,1));
        }
        blackSholesPrices.push_back(temp);
    }

    for (auto row : blackSholesPrices) {
        for (auto price : row) {
            cout << price << " ";
        }
        cout << "\n";
    }

    // Task 2
    vector<vector<double>> bisectionImpliedVol;
    for (int i = 0; i<volTable.size(); i++) {
        vector<double> temp;
        double tenor = tenors[i];
        for (int j = 0; j<volTable[i].size();j++) {
            double moneyness = moneynesses[j];
            double strike = intialVal*moneyness;
            temp.push_back(bisection(blackSholesPrices[i][j],tenor,interestRate,strike,intialVal,1,0.0000001));
        }
        bisectionImpliedVol.push_back(temp);
    }
    cout << "Bisection method vol\n";
    for (auto row : bisectionImpliedVol) {
        for (auto vol : row) {
            cout << vol << " ";
        }
        cout << "\n";
    }

    vector<vector<double>> newtonRaphsonImpliedVol;
    vector<float> intialGuessesPerMoneyness = {0.2,0.17,0.15,0.14,0.13,0.13,0.13,0.13,0.13};
    for (int i = 0; i<volTable.size(); i++) {
        vector<double> temp;
        double tenor = tenors[i];
        for (int j = 0; j<volTable[i].size();j++) {
            double moneyness = moneynesses[j];
            double strike = intialVal*moneyness;
            temp.push_back(newtonRaphson(blackSholesPrices[i][j],tenor,interestRate,strike,intialVal,1,0.0000001,intialGuessesPerMoneyness[j]));
        }
        newtonRaphsonImpliedVol.push_back(temp);
    }
    cout << "newtonRaphson vol\n";
    for (auto row : newtonRaphsonImpliedVol) {
        for (auto vol : row) {
            cout << vol << " ";
        }
        cout << "\n";
    }

    // Task 3.3

    interestRate = 0.05;
    double strike = 100;
    double initialStockPrice = 100; 
    int call = 1;
    double time = 5;
    double sigma = 0.3;

    double delta = blackScholesDelta(sigma, time,interestRate,strike,initialStockPrice,call);
    double gamma = blackScholesGamma(sigma, time,interestRate,strike,initialStockPrice);
    cout << "Delta " << delta << " Gamma " << gamma << "\n";

    return 0;
}