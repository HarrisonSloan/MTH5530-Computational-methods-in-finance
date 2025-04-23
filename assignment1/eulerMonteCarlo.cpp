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
    // Class for pricing European put options via Euler Discretisation
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
            logPrices = std::make_unique<double[]>((int)(time*stepsPerYear*paths));
            variances = std::make_unique<double[]>((int)(time*stepsPerYear*paths));
            intermediateVars = std::make_unique<double[]>((int)(time*stepsPerYear*paths));
            rng1Nums = std::make_unique<double[]>((int)(time*stepsPerYear*paths));
            rng2Nums = std::make_unique<double[]>((int)(time*stepsPerYear*paths));
            times = std::make_unique<double[]>((int)(time*stepsPerYear));
        }
    
        // Main simulation method
        // currently a small error where im not actually getting the final day computed
        void computePath(int pathNum, double stockPriceIncrement, bool reuseRNG) {
            //cout << "Thread ID " << std::this_thread::get_id() << " Computing path " << pathNum << "\n";
            // since all arrays are flattened need to find the starting index 
            int startIndex = pathNum*time*stepsPerYear;
            logPrices[startIndex] = log(initialStockPrice+stockPriceIncrement);
            variances[startIndex] = initialVariance;
            intermediateVars[startIndex] = initialVariance;
            //times[0]=0;
        
            double delta = 1 / (double)stepsPerYear;
            // review this indexing error here
            // might not be filling all the way
            for (int i = 1; i<(int)(time*stepsPerYear); i++) {
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
            mtx.lock();
            optionPricePayOffSum+=max((double)0, pow(EulerConstant, logPrices[(pathNum+1)*time*stepsPerYear-1])-strike);
            mtx.unlock();       
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

            if (reuseRNG) {
                tempFinalPrice = pow(EulerConstant,-interestRate*time) * (1.0/paths) * optionPricePayOffSum;
                cout << " Different Computed price is " << tempFinalPrice << "\n";
                return;
            }
            finalPrice = pow(EulerConstant,-interestRate*time) * (1.0/paths) * optionPricePayOffSum;
            cout << " Computed price is " << finalPrice << "\n";
            return;
        }

        void intialRun() {
            run(0,false);
            simulationRan = true;
        }

        void graph() {
            // TODO 
        }

        void deltaAndGamma() {
            // Delta = F(t, intialPrice + small amount) - F(t, intialPrice - small amount)  / 2h
            // Gama = F(t, intialPrice + small amount) - F(t, intialPrice - small amount) + 2 * oriingla Price / h^2
            if (!simulationRan) {
                throw("Run intial simulation");
            }
            double normal = finalPrice;
            double change = 0.01;
            run(change,true);
            double high = tempFinalPrice;
            run(-change,true);
            double low = tempFinalPrice;
            double delta = (high - low) / (2.0 * change);
            double gamma = (high + low - 2.0 * normal) / (change * change);
            cout << " delta is " << delta << " Gamma is " << gamma << "\n";
        }
    private:
        // random number generations -> alternatively we apply box muller transform to achieve same result
        mutex mtx; // lock for random number generation
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
    int stepsPerYear = 40;
    int numPaths = 40000;
    HestonSimulator absorptionSim = HestonSimulator(time,interestRate,strike,initialStockPrice,initialVariance,longTermVariance,corr,rateOfReversion,volOfVol,stepsPerYear,numPaths,identity, positivePart, positivePart);\
    absorptionSim.deltaAndGamma();
    absorptionSim.intialRun();
    // std::vector<std::vector<double>> allPaths(numPaths);
    // std::vector<double> timeVec(stepsPerYear*time);
    // for (int i = 0; i < (int)(time*stepsPerYear); i++) {
    //     timeVec[i] = i * (1 / (double)stepsPerYear);
    // }

    // for (int i = 0; i < numPaths; i++) {
    //     absorptionSim.run(positivePart, positivePart, positivePart);
    //     allPaths[i] = std::vector<double>(absorptionSim.logPrices.get(), absorptionSim.logPrices.get() + (int)(stepsPerYear*time));
    // }

    // matplot::figure()->size(1920,1080);
    // matplot::hold(true);  // keep all lines on the same figure
    // for (const auto& path : allPaths) {
    //     matplot::plot(timeVec, path);
    // };
    // matplot::title("Monte Carlo sim Pricing");
    // matplot::xlabel("M-time steps");
    // matplot::ylabel("log option prices");
    // matplot::show();
    // matplot::save("plot.png");
    return 0;
}