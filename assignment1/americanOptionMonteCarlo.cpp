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

class HestonSimulator {
    // Class for pricing European call/put (call implemented) options via Euler Discretisation
    public:
        // Parameters
        const double EulerConstant = std::exp(1.0);
        double interestRate, strike, initialStockPrice, initialVariance;
        double longTermVariance, corr, rateOfReversion, volOfVol, time;
        int paths;
        int stepsPerYear;
        bool simulationRan = false;
        std::function<double(double)> f1;
        std::function<double(double)> f2;
        std::function<double(double)> f3;
    
        // Output
        std::unique_ptr<double[]> logPrices;
        std::unique_ptr<double[]> variances;
        std::unique_ptr<double[]> intermediateVars;
        std::unique_ptr<double[]> times;
        std::unique_ptr<double[]> rng1Nums;
        std::unique_ptr<double[]> rng2Nums;
        double optionPricePayOffSum=0;
        double executionTime;
        double finalPrice;
        double tempFinalPrice;

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
                        int stepsPerYear,
                        int paths,
                        std::function<double(double)> f1,
                        std::function<double(double)> f2,
                        std::function<double(double)> f3,
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
            this->stepsPerYear = stepsPerYear;
            this->paths = paths;
            this->f1 = f1;
            this->f2 = f2;
            this->f3 = f3;
            // "Auto" deallocate off heap with smart pointers
            logPrices = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            variances = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            intermediateVars = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            rng1Nums = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            rng2Nums = std::make_unique<double[]>((int)(time*stepsPerYear*paths+paths));
            times = std::make_unique<double[]>((int)(time*stepsPerYear+1));
        }
    
        // Main simulation method
        // currently a small error where im not actually getting the final day computed
        void computePath(int pathNum, double stockPriceIncrement, bool reuseRNG) {
            //cout << "Thread ID " << std::this_thread::get_id() << " Computing path " << pathNum << "\n";
            // since all arrays are flattened need to find the starting index 
            int startIndex = pathNum*time*stepsPerYear + pathNum;
            logPrices[startIndex] = log(initialStockPrice+stockPriceIncrement);
            variances[startIndex] = initialVariance;
            intermediateVars[startIndex] = initialVariance;
            //times[0]=0;
        
            double delta = 1 / (double)stepsPerYear;
            // review this indexing error here
            // might not be filling all the way
            for (int i = 1; i<=(int)(time*stepsPerYear); i++) {
                // generate and store random generated numbers for Greeks
                if (!reuseRNG) {
                    mtx.lock();
                    rng1Nums[startIndex+i] = gen1();
                    rng2Nums[startIndex+i] = gen2();
                    mtx.unlock();
                }
                // Calculate the Wiener increments for correlated procresses using chomsky decompistion
                double deltaW1 = sqrt(delta)*rng1Nums[startIndex+i];
                double deltaW2 = corr*deltaW1 + sqrt(1-pow(corr,2)) * sqrt(delta)*rng2Nums[startIndex+i];
                // Apply formula 
                intermediateVars[startIndex+i] = f1(intermediateVars[startIndex+i-1]) - rateOfReversion*delta*(f2(intermediateVars[startIndex+i-1])-longTermVariance) + volOfVol*sqrt(f3(intermediateVars[startIndex+i-1]))*deltaW1;
                variances[startIndex+i] = f3(intermediateVars[i]);
                logPrices[startIndex+i] = logPrices[startIndex+i-1] + (interestRate-0.5*variances[startIndex+i])*delta + sqrt(variances[startIndex+i])*deltaW2;
                //
                //times[i] = times[i-1]+delta;
            }
            // Retrieve final prices add to running total and average it    
        }

        // Note we may start to have errors where the time steps per year into how many "years" we have is no longer an integer
        // say time is 1.3 * 15 steps a year = 19.5
        // I guess we need to round down and account for this 

        void run(double stockPriceIncrement,bool reuseRNG) {
            optionPricePayOffSum = 0.0;
            const size_t threadCount =  std::thread::hardware_concurrency();
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
                        computePath(j,stockPriceIncrement,reuseRNG);
                    };
                });
                runningSum += pathsPerThread[i];
            }
            // Join the threads
            for (auto& thread : threads) {
                thread.join();
            }

            // Now need to backwards induct on the set of in path 
            // So keep of track of ones that are firstly being 
            // Continued, if continued and in money then include in regression


            // Logic!!!!
            // Check if its in path (K-discounted at S_t)^+, if this is more than 0, we will include in the regresison

            // what I need, latest g_k, the associated times this occured tao that sits wit hteh tao_k 
            // Each time find hte in money paths, retrieve their g_k and discount by the tao_k amount (the time difference between tao_j and where we are now)
            // complete a linear regression for hte subset
            // tricky bit is how do I correctly map back the values 
            // I then need to see what the predicted values are, if its better update g_k and tao_k

            // Payoff values
            auto terminalPayoffs = std::make_unique<double[]>(paths);
            auto terminalTimes = std::make_unique<double[]>(paths);
            for (int i = 0; i<paths; i++) {
                // Go to the end price
                terminalPayoffs[i] = putPayoff(this->logPrices[(i+1)*this->time*this->stepsPerYear + i]);
                // Time is measured here in just amount of intervals but corrected down the line
                terminalTimes[i]=this->time*this->stepsPerYear;
            }
            

            // intialise the above

            // In money subset 
            // regression pointers, this is to say what path is being used in the regression
            // For example if we have 10 paths but only 4 used in the regressions say (0,5,6,8)
            
            auto regressionPointers = std::make_unique<int[]>(paths);

            const int features = 6;
            // regression will be a 
            // # in money paths * 6 matrix 
            // int i = (int)(time*stepsPerYear) vs int i = (int)(time*stepsPerYear)-1
            double delta = 1 / (double)this->stepsPerYear;
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

                    // Create the X matrix and Y matrix for the regression
                }
                cout << "Total in money paths " << numPathInMoney << "\n";
                Eigen::MatrixXd X(numPathInMoney, features);
                Eigen::VectorXd Y(numPathInMoney);
                for (int pathNum = 0; pathNum<paths;pathNum++) {
                    if (regressionPointers[pathNum] != -1) {
                            int pathIndex = pathNum*time*stepsPerYear + pathNum + i;
                            double price = pow(this->EulerConstant,this->logPrices[pathIndex]);
                            double variance = variances[pathIndex];
                            X(regressionPointers[pathNum],0) = 1;
                            X(regressionPointers[pathNum],1) = price; // S
                            X(regressionPointers[pathNum],2) = pow(price,2); // S^2
                            X(regressionPointers[pathNum],3) = variance;
                            X(regressionPointers[pathNum],4) = pow(variance,2);
                            X(regressionPointers[pathNum],5) = price*variance;
                            // Discounted g_x
                            Y(regressionPointers[pathNum]) = pow(this->EulerConstant,-this->interestRate*(terminalTimes[pathNum]-i)*delta)*terminalPayoffs[pathNum];
                    }
                }
                // Run the regression 
                Eigen::VectorXd beta = X.colPivHouseholderQr().solve(Y);
                Eigen::VectorXd continuationValues = X * beta;
                cout << " Betas 1 " << beta(0) << " Betas 2 " << beta(1) << " Betas 3 " << beta(2) << " Betas 4 " << beta(3) << " Betas 5 " << beta(4) << " Betas 6 " << beta(5) << "\n";
                // Compute R²
                double yMean = Y.mean();
                double ssTotal = (Y.array() - yMean).square().sum();
                double ssResidual = (Y - continuationValues).squaredNorm();  // same as ∑ (Y_i - Ŷ_i)^2

                double rSquared = 1.0 - (ssResidual / ssTotal);

                // Optional: Clamp to [0, 1] in case of tiny numerical noise
                rSquared = std::max(0.0, std::min(1.0, rSquared));

                std::cout << "R² at timestep " << i << " = " << rSquared << "\n";
                // Determine the set where the payoff is larger than the continuation value
                for (int pathNum = 0; pathNum<paths;pathNum++) {
                    if (regressionPointers[pathNum] != -1) {
                        int pathIndex = pathNum*time*stepsPerYear + pathNum + i;
                        double newPayOff = putPayoff(this->logPrices[pathIndex]);
                        if (continuationValues(regressionPointers[pathNum]) < newPayOff) {
                            terminalPayoffs[pathNum] = newPayOff;
                            terminalTimes[pathNum] = i;
                        }
                    }
                    // reset the regression pointer
                    regressionPointers[pathNum] = -1;
                }
            }

            double continuationValuesSum = 0;
            for (int pathNum = 0; pathNum<paths;pathNum++) {
                continuationValuesSum += pow(this->EulerConstant,-this->interestRate*(terminalTimes[pathNum])*delta)*terminalPayoffs[pathNum];
            }
            double continuationAverage = continuationValuesSum / this->paths;
            finalPrice = max(max(0.0,this->strike-initialStockPrice),continuationAverage);
            cout << "Final price is " << finalPrice << "\n";
    
            
        }

        void intialRun() {
            run(0,false);
            simulationRan = true;
        }

        void graph() {
            // TODO 
        }

    private:
        // random number generations -> alternatively we apply box muller transform to achieve same result
        mutex mtx; // lock for random number generation
        boost::mt19937 rng1;
        boost::mt19937 rng2;
        boost::normal_distribution<> nd;
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > gen1;
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > gen2;

        double putPayoff(double logPrice) {
            return max(0.0,this->strike-pow(this->EulerConstant,logPrice));
        }
        double basisFunc(double a, double b, double c, double d, double e, double f, double price, double variance) {
            return a+b*price+c*pow(price,2)+d*variance+d*pow(variance,2)+f*price*variance;
        }
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
    int stepsPerYear = 50;
    int numPaths = 10000;
    HestonSimulator absorptionSim = HestonSimulator(time,interestRate,strike,initialStockPrice,initialVariance,longTermVariance,corr,rateOfReversion,volOfVol,stepsPerYear,numPaths,identity, positivePart, positivePart);\
    absorptionSim.intialRun();
    matplot::figure()->size(1920,1080);
    // matplot::hold(true);  // keep all lines on the same figure
    // for (const auto& path : allPaths) {
    //     matplot::plot(timeVec, path);
    // };
    // matplot::title("Monte Carlo sim Pricing");
    // matplot::xlabel("M-time steps");
    // matplot::ylabel("log option prices");
    matplot::show();
    // matplot::save("plot.png");
    return 0;
}