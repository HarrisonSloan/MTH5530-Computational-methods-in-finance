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
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;

double identity(double x) {
    return x;
}

double positivePart(double x) {
    return std::max(0.0, x);
}

double fabs(double x) {
    return fabs(x);
}

class AmericanPutHeston {
    // Class for pricing American put options!
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
    
        // Output
        std::unique_ptr<double[]> logPrices;
        std::unique_ptr<double[]> variances;
        std::unique_ptr<double[]> intermediateVars;
        std::unique_ptr<double[]> times;
        std::unique_ptr<double[]> rng1Nums;
        std::unique_ptr<double[]> rng2Nums;
        double optionPricePayOffSum=0;
        std::chrono::duration<signed long, std::micro> simulationTime;
        double finalPutPrice;
        double tempFinalPrice;

        // Constructor
        // Constructor
        AmericanPutHeston(double time, 
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
                        size_t threadCount =  std::thread::hardware_concurrency(),
                        unsigned int basedSeed = 1337) : 
                        
                        nd(0.0, 1.0),
                        engines(threadCount*2)
    {       
        // Intialise model parameters
        this->interestRate = interestRate;
        this->time = time;
        this->strike = strike;
        this->initialStockPrice = initialStockPrice;
        this->longTermVariance = longTermVariance;
        this->corr = corr;
        this->rateOfReversion = rateOfReversion;
        this->volOfVol = volOfVol;
        this->initialVariance=initialVariance;
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
        // Push the engines
        // Then use each engine to create a generator
        for (std::size_t i = 0; i < threadCount*2; ++i) {
            engines.emplace_back(boost::random::mt19937{basedSeed + static_cast<unsigned int>(i)});
            generators.emplace_back(engines.back(), nd);
        }

        // Create flat arrays for each paths data
        logPrices = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
        variances = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
        intermediateVars = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
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
            for (int i = 1; i<=(int)(time*stepsPerYear); i++) {
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
        }

        // Note we may start to have errors where the time steps per year into how many "years" we have is no longer an integer
        // say time is 1.3 * 15 steps a year = 19.5
        // I guess we need to round down and account for this 

        void run(double stockPriceIncrement,bool reuseRNG) {
            optionPricePayOffSum = 0.0;
            std::cout << "Number of threads = " <<  threadCount << std::endl;
            std::vector<std::thread> threads;
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
            threads.clear();
            // Begin backwards induction

            // Payoff values
            auto terminalPayoffs = std::make_unique<double[]>(paths);
            auto terminalTimes = std::make_unique<double[]>(paths);
            for (int i = 0; i<paths; i++) {
                // Go to the end price
                terminalPayoffs[i] = putPayoff(this->logPrices[(i+1)*this->time*this->stepsPerYear + i]);
                // Time is measured here in just amount of intervals but corrected down the line
                terminalTimes[i]=this->time*this->stepsPerYear;
            }
        

            // Use to keep track of the paths that are in money
            auto regressionPointers = std::make_unique<int[]>(paths);
            const int features = 6;
            double delta = 1 / (double)this->stepsPerYear;
            // Look over each time step
            for (int i = (int)(time*stepsPerYear)-1; i>=1; i--) {

                // For each path determine if the new price is in money
                int numPathInMoney = 0;
                for (int pathNum = 0; pathNum<paths;pathNum++) {
                    int pathIndex = pathNum*time*stepsPerYear + pathNum + i;
                    // if an in money path
                    if (putPayoff(this->logPrices[pathIndex])>0) {
                        regressionPointers[pathNum] = numPathInMoney;
                        numPathInMoney++;
                    } else {
                        regressionPointers[pathNum] = -1;
                    }   
                }

                // Create the X matrix and Y matrix for the regression
                Eigen::MatrixXd X(numPathInMoney, features);
                Eigen::VectorXd Y(numPathInMoney);
                for (int pathNum = 0; pathNum<paths;pathNum++) {
                    // If path is in money
                    if (regressionPointers[pathNum] != -1) {
                            int pathIndex = pathNum*time*stepsPerYear + pathNum + i;
                            double price = exp(this->logPrices[pathIndex]);
                            double variance = variances[pathIndex];
                            X(regressionPointers[pathNum],0) = 1;
                            X(regressionPointers[pathNum],1) = price;
                            X(regressionPointers[pathNum],2) = price*price;
                            X(regressionPointers[pathNum],3) = variance;
                            X(regressionPointers[pathNum],4) = variance*variance;
                            X(regressionPointers[pathNum],5) = price*variance;
                            // Discounted g_x (between current time and previous time)
                            Y(regressionPointers[pathNum]) = exp(-this->interestRate*(terminalTimes[pathNum]-i)*delta)*terminalPayoffs[pathNum];
                    }
                }
                // Run the regression 
                // Recommend method from https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
                Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
                Eigen::VectorXd beta = svd.solve(Y);
                Eigen::VectorXd continuationValues = X * beta;

                // Regression details 
                // cout << " Betas 1 " << beta(0) << " Betas 2 " << beta(1) << " Betas 3 " << beta(2) << " Betas 4 " << beta(3) << " Betas 5 " << beta(4) << " Betas 6 " << beta(5) << "\n";
                // // Compute R²
                // double yMean = Y.mean();
                // double ssTotal = (Y.array() - yMean).square().sum();
                // double ssResidual = (Y - continuationValues).squaredNorm();  // same as ∑ (Y_i - Ŷ_i)^2

                // double rSquared = 1.0 - (ssResidual / ssTotal);
                // // Optional: Clamp to [0, 1] in case of tiny numerical noise
                // rSquared = std::max(0.0, std::min(1.0, rSquared));

                // std::cout << "R² at timestep " << i << " = " << rSquared << "\n";

                // Determine the set where the payoff is larger than the continuation value then update their payoffs and time of termination
                for (int pathNum = 0; pathNum<paths;pathNum++) {
                    if (regressionPointers[pathNum] != -1) {
                        int pathIndex = pathNum*time*stepsPerYear + pathNum + i;
                        double newPayOff = putPayoff(this->logPrices[pathIndex]);
                        if (continuationValues(regressionPointers[pathNum]) < newPayOff) {
                            terminalPayoffs[pathNum] = newPayOff;
                            terminalTimes[pathNum] = i;
                        }
                    }
                    // reset the regression pointers
                    regressionPointers[pathNum] = -1;
                }
            }
            // Sum and average optimal payoffs
            double continuationValuesSum = 0;
            for (int pathNum = 0; pathNum<paths;pathNum++) {
                continuationValuesSum += exp(-this->interestRate*(terminalTimes[pathNum])*delta)*terminalPayoffs[pathNum];
            }
            double continuationAverage = continuationValuesSum / this->paths;
            finalPutPrice = max(max(0.0,this->strike-initialStockPrice),continuationAverage);
            cout << "Final price is " << finalPutPrice << "\n";    
        }

        void intialRun() {
            auto start = chrono::high_resolution_clock::now();
            run(0,false);
            simulationRan = true;
            auto end = chrono::high_resolution_clock::now();
            simulationTime = chrono::duration_cast<chrono::microseconds>(end - start);
        }

        constexpr int64_t getRunTime() {
            // Note this is only for intial runs
            if (!simulationRan) {
                throw("Run intial simulation");
            }
            return simulationTime.count();
        }

        double getPutOptionPrice() {
            if (!simulationRan) {
                throw("Run intial simulation");
            }
            return finalPutPrice;
        }

    private:
        // random number generations -> alternatively we apply box muller transform to achieve same result
        mutex mtx; 
        // List of engines, normal distribution (0,1) and the generators list
        std::vector<boost::mt19937> engines;
        boost::normal_distribution<> nd;
        std::vector<boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >> generators;

        double generate(std::size_t i) {
            return generators[i]();
        }
        double putPayoff(double logPrice) {
            return max(0.0,this->strike-exp(logPrice));
        }
        double basisFunc(double a, double b, double c, double d, double e, double f, double price, double variance) {
            return a+b*price+c*pow(price,2)+d*variance+d*pow(variance,2)+f*price*variance;
        }
    };

int main() {
    // Parameters
    double interestRate = 0.05;
    double strike = 100;
    double initialStockPrice = 100; 
    double initialVariance = 0.09;
    double longTermVariance = 0.09;
    double corr = -0.3;
    double rateOfReversion = 2; 
    double volOfVol = 1;
    double time = 5;
    int stepsPerYear = 50;
    int numPaths = 10000;
    
    AmericanPutHeston FullTruncationSim = AmericanPutHeston(time,interestRate,strike,initialStockPrice,initialVariance,longTermVariance,corr,rateOfReversion,volOfVol,stepsPerYear,numPaths,identity, positivePart, positivePart);\
    FullTruncationSim.intialRun();
    cout << "Total time " << FullTruncationSim.getRunTime() << " Final price " << FullTruncationSim.getPutOptionPrice() << "\n";
    return 0;
}