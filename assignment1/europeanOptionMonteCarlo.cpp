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
#include <thread>
#include <mutex>
#include <boost/random.hpp>
#include <chrono>

using namespace std;

double identity(double x) {
    return x;
}

double dabs(double x) {
    return fabs(x);
}

double positivePart(double x) {
    return std::max(0.0, x);
}

class EuropeanOptionHeston {
    // Class for pricing European call/put (call implemented) options via Euler Discretisation
    public:
        // Parameters
        double interestRate, strike, initialStockPrice, initialVariance;
        double longTermVariance, corr, rateOfReversion, volOfVol, time;
        int paths;
        int stepsPerYear;
        bool simulationRan = false;
        std::function<double(double)> f1;
        std::function<double(double)> f2;
        std::function<double(double)> f3;
        size_t threadCount;
        std::chrono::duration<signed long, std::micro> simulationTime;
    
        // Output
        std::unique_ptr<double[]> logPrices;
        std::unique_ptr<double[]> variances;
        std::unique_ptr<double[]> intermediateVars;
        std::unique_ptr<double[]> times;
        std::unique_ptr<double[]> rng1Nums;
        std::unique_ptr<double[]> rng2Nums;
        double optionPricePayOffSum=0;
        double executionTime;
        double finalCallPrice;
        double finalPutPrice;
        double tempFinalPrice;

        // Constructor
        EuropeanOptionHeston(double time, 
                        double interestRate, 
                        double strike, 
                        double initialStockPrice, 
                        double initialVariance, 
                        double longTermVariance, 
                        double corr, 
                        double rateOfReversion, 
                        double volOfVol, 
                        int stepsPerYear,
                        int paths,
                        std::function<double(double)> f1,
                        std::function<double(double)> f2,
                        std::function<double(double)> f3,
                        unsigned int basedSeed = 1337,
                        size_t threadCount =  std::thread::hardware_concurrency()) : 
                        
                        nd(0.0, 1.0),
                        engines(threadCount*2)
        {
            // Intialise model parameters
            this->interestRate = interestRate;
            this->time = time;
            this->strike = strike;
            this->initialStockPrice = initialStockPrice;
            this->initialVariance=initialVariance;
            this->longTermVariance = longTermVariance;
            this->corr = corr;
            this->rateOfReversion = rateOfReversion;
            this->volOfVol = volOfVol;
            this->stepsPerYear = stepsPerYear;
            this->paths = paths;

            // Intialise functions for specific scheme
            this->f1 = f1;
            this->f2 = f2;
            this->f3 = f3;

            // Intialise thread count
            if (threadCount > std::thread::hardware_concurrency()) {
                throw("Computer doesnt have this many threads, use less");
            }
            this->threadCount=threadCount;

            // Create generator
            // Push the engines (the need for this is so threads dont produce the same paths, they need a different underlying engine)
            // Then use each engine to create a generator
            for (std::size_t i = 0; i < threadCount*2; ++i) {
                engines.emplace_back(boost::random::mt19937{basedSeed + static_cast<unsigned int>(i)});
                generators.emplace_back(engines.back(), nd);
            }

            // Create flat arrays for each paths data
            logPrices = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            variances = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            intermediateVars = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            // Store random variables for consistency of simulation if you want to alter other inputs (price, strike, ect)
            rng1Nums = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            rng2Nums = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            // Time for graphing the X axis
            times = std::make_unique<double[]>((int)(time*stepsPerYear+1));
            // intialise the time array here
        }
    
        // Main simulation method
        void computePath(int pathNum, double stockPriceIncrement, bool reuseRNG, int threadNum) {
            // since all arrays are flattened need to find the starting index 
            int startIndex = pathNum*time*stepsPerYear + pathNum;
            logPrices[startIndex] = log(initialStockPrice+stockPriceIncrement);
            variances[startIndex] = initialVariance;
            intermediateVars[startIndex] = initialVariance;
            double delta = 1 / (double)stepsPerYear;
            double sumTime = 0;
            for (int i = 1; i<=(int)(time*stepsPerYear); i++) {
                sumTime = i*delta;
                // generate and store random generated numbers for Greeks
                if (!reuseRNG) {
                    rng1Nums[startIndex+i] = generate(2*threadNum); // take from the threads 1st generator
                    rng2Nums[startIndex+i] = generate(2*threadNum + 1); // take from the threasd 2nd generator
                }
                // Calculate the Wiener increments for correlated procresses using chomsky decompistion
                double deltaW1 = sqrt(delta)*rng1Nums[startIndex+i];
                double deltaW2 = corr*deltaW1 + sqrt(1-corr*corr) * sqrt(delta)*rng2Nums[startIndex+i];
                // Apply formula 
                intermediateVars[startIndex+i] = f1(intermediateVars[startIndex+i-1]) - rateOfReversion*delta*(f2(intermediateVars[startIndex+i-1])-longTermVariance) + volOfVol*sqrt(f3(intermediateVars[startIndex+i-1]))*deltaW1;
                variances[startIndex+i] = f3(intermediateVars[startIndex+i]);
                logPrices[startIndex+i] = logPrices[startIndex+i-1] + (interestRate-0.5*variances[startIndex+i-1])*delta + sqrt(variances[startIndex+i-1])*deltaW2;
            }
            // Retrieve final prices add to running total and average it
            mtx.lock();
            optionPricePayOffSum+=max(0.0, exp(logPrices[startIndex+(int)(time*stepsPerYear)])-strike);
            mtx.unlock();       
        }

        // Note we may start to have errors where the time steps per year into how many "years" we have is no longer an integer
        // say time is 1.3 * 15 steps a year = 19.5
        // I guess we need to round down and account for this

        void run(double stockPriceIncrement,bool reuseRNG) {
            // Create variable for threads to sum too
            optionPricePayOffSum = 0.0;
            //std::cout << "Number of threads = " <<  threadCount << std::endl;
            std::vector<std::thread> threads;
            // Alter to allocate the thread space automatically based on the input
            threads.reserve(threadCount);

            // Evenly divide the work
            vector<int> pathsPerThread(threadCount, this->paths / threadCount);
            // Distribute remaining paths to the first threads
            for (int i = 0; i < this->paths % threadCount; i++) {
                pathsPerThread[i]++;
            } 
            int runningSum = 0;
            // Call thread function for each thread
            for (int i = 0; i < threadCount; i++) {
                int start = runningSum;
                int end = runningSum + pathsPerThread[i] - 1;
                threads.emplace_back([=]() {
                    for (int j = start ; j <= end ; j ++) {
                        computePath(j,stockPriceIncrement,reuseRNG,i);
                    };
                });
                runningSum += pathsPerThread[i];
            }
            // Join the threads
            for (auto& thread : threads) {
                thread.join();
            }
            // If reusing dont alter main sims price
            if (reuseRNG) {
                tempFinalPrice = exp(-interestRate*time) * (1.0/paths) * optionPricePayOffSum;
                return;
            }
            finalCallPrice = exp(-interestRate*time) * (1.0/paths) * optionPricePayOffSum;
            finalPutPrice = finalCallPrice - initialStockPrice + strike*exp(-interestRate*time);
        }

        void intialRun() {
            auto start = chrono::high_resolution_clock::now();
            run(0,false);
            simulationRan = true;
            auto end = chrono::high_resolution_clock::now();
            simulationTime = chrono::duration_cast<chrono::microseconds>(end - start);
        }
 

        void graph(int typeOfGraph) {
            if (!simulationRan) {
                throw("Run intial simulation");
            }
            if (paths>500) {
                throw("To many paths to graphs try a lower number");
            }
            // Create the time vector
            std::vector<double> timeVec(stepsPerYear*time+1);
            for (int i = 0; i <= (int)(time*stepsPerYear); i++) {
                timeVec[i] = i * (1 / (double)stepsPerYear);
            }
            // For log prices
            if (typeOfGraph==1) {
                std::vector<std::vector<double>> allPathsPrices(paths);
                for (int i = 0; i < paths; i++) {
                    int startIndex = i*time*stepsPerYear + i;
                    allPathsPrices[i] = std::vector<double>(logPrices.get() + startIndex, logPrices.get() + startIndex+ (int)(stepsPerYear*time)+1);
                }
                matplot::figure()->size(1920,1080);
                matplot::hold(true);  // keep all lines on the same figure
                for (const auto& path : allPathsPrices) {
                    matplot::plot(timeVec, path);
                };
                matplot::title("European Call option log price paths");
                matplot::xlabel("200 steps");
                matplot::ylabel("log option prices");
                matplot::show();
                matplot::save("logPricesAbsorption_2.png");
                return;
            }
            // For variances
            if (typeOfGraph==2) {
                std::vector<std::vector<double>> allPathsVariances(paths);
                for (int i = 0; i < paths; i++) {
                    int startIndex = i*time*stepsPerYear + i;
                    allPathsVariances[i] = std::vector<double>(variances.get() + startIndex, variances.get() + startIndex+ (int)(stepsPerYear*time)+1);
                }
                matplot::figure()->size(1920,1080);
                matplot::hold(true);  // keep all lines on the same figure
                for (const auto& path : allPathsVariances) {
                    matplot::plot(timeVec, path);
                };
                matplot::title("European Call option variance paths");
                matplot::xlabel("200 steps");
                matplot::ylabel("variances prices");
                matplot::show();
                matplot::save("variancesAbsorption_2.png");
                return;
            }
            // For intermediate variances
            if (typeOfGraph==3) {
                std::vector<std::vector<double>> allPathsIntermediateVariances(paths);
                for (int i = 0; i < paths; i++) {
                    int startIndex = i*time*stepsPerYear + i;
                    allPathsIntermediateVariances[i] = std::vector<double>(intermediateVars.get() + startIndex, intermediateVars.get() + startIndex+ (int)(stepsPerYear*time)+1);
                }
                matplot::figure()->size(1920,1080);
                matplot::hold(true);  // keep all lines on the same figure
                for (const auto& path : allPathsIntermediateVariances) {
                    matplot::plot(timeVec, path);
                };
                matplot::title("European Call option intermediate variance paths");
                matplot::xlabel("200 steps");
                matplot::ylabel("Intermediate variances prices");
                matplot::show();
                matplot::save("intermediateVariancesAbsorption_2.png");
                return;
            }
        }

        void deltaAndGamma() {
            // Delta = F(t, intialPrice + small amount) - F(t, intialPrice - small amount)  / 2h
            // Gama = F(t, intialPrice + small amount) - F(t, intialPrice - small amount) + 2 * oriingla Price / h^2
            if (!simulationRan) {
                throw("Run intial simulation");
            }
            double normal = finalCallPrice;
            double change = 0.001*initialStockPrice;
            run(change,true);
            double high = tempFinalPrice;
            run(-change,true);
            double low = tempFinalPrice;
            double delta = (high - low) / (2.0 * change);
            double gamma = (high + low - 2.0 * normal) / (change * change);
            cout << " delta is " << delta << " Gamma is " << gamma << "\n";
        }

        constexpr int64_t getRunTime() {
            // Note this is only for intial runs
            if (!simulationRan) {
                throw("Run intial simulation");
            }
            return simulationTime.count();
        }

        double getCallOptionPrice() {
            if (!simulationRan) {
                throw("Run intial simulation");
            }
            return finalCallPrice;
        }

        double getPutOptionPrice() {
            if (!simulationRan) {
                throw("Run intial simulation");
            }
            return finalPutPrice;
        }
    private:
        // random number generations -> alternatively we apply box muller transform to achieve same result less overhead, but trust robustness of this more
        mutex mtx; 
        // List of engines, normal distribution (0,1) and the generators list
        std::vector<boost::mt19937> engines;
        boost::normal_distribution<> nd;
        std::vector<boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >> generators;

        double generate(std::size_t i) {
            return generators[i]();
        }
        
    };

double averageSimPrice(std::vector<double> monteCarloPrices) {
    return accumulate(monteCarloPrices.begin(),monteCarloPrices.end(),0.0) / monteCarloPrices.size();
}

double standardError(std::vector<double> monteCarloPrices) {
    double avgSimPrice = averageSimPrice(monteCarloPrices);
    double sum = 0.0;
    for (double price: monteCarloPrices) {
        sum+=pow(price-avgSimPrice,2);
    }
    return sqrt(sum/monteCarloPrices.size());
}

double bias(std::vector<double> monteCarloPrices, double realisedPrice) {
    return fabs(averageSimPrice(monteCarloPrices)-realisedPrice);
}

double rootMeanSquareError(std::vector<double> monteCarloPrices, double realisedPrice) {
    return sqrt(pow(bias(monteCarloPrices,realisedPrice),2)+pow(standardError(monteCarloPrices),2));
}

int main() {
    // Parameters
    double interestRate = 0.05;
    double strike = 100;
    double initialStockPrice = 100; 
    double initialVariance = 0.09;
    double longTermVariance = 0.09;
    double corr = -0.3;
    double rateOfReversion = 2; 
    double volOfVol = 0.3;
    double time = 5;
    int stepsPerYear = 40;
    int numPaths = 100;
    double realisedPrice = 34.9998;

    int totalRuns = 100;
    // same set of sets across different methods
    srand(1337);
    std::vector<int> seeds;
    for (int i = 0; i<totalRuns; i++) {
        seeds.push_back(rand());
    }
    cout << "Generated the prices\n";
    std::vector<double> simPrices;
    int64_t CumlativeTime = 0;
    for (int i = 0; i<totalRuns; i++) {

        // Choose which method you want, comment out the rest

        // Absoprtion
        EuropeanOptionHeston sim = EuropeanOptionHeston(time,interestRate,strike,initialStockPrice,initialVariance,longTermVariance,corr,rateOfReversion,volOfVol,stepsPerYear,numPaths,positivePart, positivePart, positivePart,seeds[i]);
        
        // Reflection
        //EuropeanOptionHeston sim = EuropeanOptionHeston(time,interestRate,strike,initialStockPrice,initialVariance,longTermVariance,corr,rateOfReversion,volOfVol,stepsPerYear,numPaths,dabs, dabs, dabs,seeds[i]);
        
        // Partial truncation
        //EuropeanOptionHeston sim = EuropeanOptionHeston(time,interestRate,strike,initialStockPrice,initialVariance,longTermVariance,corr,rateOfReversion,volOfVol,stepsPerYear,numPaths,identity, identity, positivePart,seeds[i]);
        
        // Full truncation
        //EuropeanOptionHeston sim = EuropeanOptionHeston(time,interestRate,strike,initialStockPrice,initialVariance,longTermVariance,corr,rateOfReversion,volOfVol,stepsPerYear,numPaths,identity, positivePart, positivePart,seeds[i]);
        sim.intialRun();
        double price = sim.getCallOptionPrice();
        //cout << "Price : " << price << " Time : " << absorptionSim.getRunTime() << "\n";
        CumlativeTime+=sim.getRunTime();
        simPrices.push_back(price);
        //absorptionSim.deltaAndGamma();
    }
    cout << "Bias " << bias(simPrices,realisedPrice) << " Standard error : " << standardError(simPrices) << " RMSE : " << rootMeanSquareError(simPrices,realisedPrice) << " Average run time " << CumlativeTime/totalRuns << "\n";
    return 0;
}