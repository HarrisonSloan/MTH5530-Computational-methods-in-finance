#include <iostream>
using namespace std;
#include <iomanip> // For formatting output
#include <vector>
#include <cmath>
#include <iostream>

// Templated class for Triangle (supports int, float, etc.)
template <typename T>
class Triangle {
private:
    int levels;
    std::vector<T> data;  // Underlining data is a single vector

    // Compute base index of level L
    int baseIndex(int L) const {
        return (L * (L + 1)) / 2;
    }

public:
    // Constructor (supports both int and float)
    Triangle(int m) {
        levels = m;
        data = std::vector<T>((m * (m + 1)) / 2, 0); // Default value 0
    }

    // Get value at (level, position)
    T get(int L, int P) const {
        return data[baseIndex(L) + P];
    }

    // Set value at (level, position)
    void set(int L, int P, T value) {
        data[baseIndex(L) + P] = value;
    }
    // Get the left child
    T getLeftChild(int L, int P) const {
        return data[baseIndex(L+1) + P];
    }
    // Get the right child
    T getRightChild(int L, int P) const {
        return data[baseIndex(L+1) + P+1];
    }

    void prettyPrint() const {
        int index = 0;
        int maxWidth = 6 * levels + 3* levels; // Controls spacing for alignment
    
        for (int L = 0; L < levels; L++) {
            int numElements = L + 1;
            
            // Compute proper padding to center-align the row
            int rowWidth = numElements * 5 + (numElements - 1) * 3; // Element width + spacing
            int padding = (maxWidth - rowWidth) / 2;
    
            // Print leading spaces for centering
            std::cout << std::string(padding, ' ');
    
            // Print row elements with spacing
            for (int P = 0; P <= L; P++) {
                std::cout << std::setw(5) << std::left << std::fixed << std::setprecision(2)
                          << data[index++] << "   ";
            }
            std::cout << "\n";
        }
    }
};

float binomialTreeEuropeOptionPricing(float sigma, float direction, int steps, int time, float interestRate, float strike, float intialStockPrice, int putOrCall) {
    const double EulerConstant = std::exp(1.0);
    // Calculate p, u and d
    float deltaTime = (float)time/steps;
    // Calculate p, u and d (general quadratic equation approach)
    // We dont assume a specific form of u and d, we decide this with direction (and therefore beta)
    // float beta = 0.5*(direction*pow(EulerConstant,-interestRate*deltaTime)+pow(EulerConstant,(interestRate+pow(sigma,2))*deltaTime));
    // float u = beta + sqrt(pow(beta,2)-direction);
    // float d = direction/u;

    // Alternative calculations of u and d using Cox-Ross-Rubinstein 
    // Ensures the volatility of the binomial tree model matches the variance of the continous Black-Schoels model
    // Gurantees recombininig trees -> the path in the tree doesnt inform the outcome
    float u = pow(EulerConstant,sigma*sqrt(deltaTime));
    float d = pow(EulerConstant,-sigma*sqrt(deltaTime));

    // calculate p
    float probability = (pow(EulerConstant,interestRate*deltaTime)-d)/(u-d);

    cout << "Delta time " << deltaTime << "\n";
    cout << "u " << u << "\n";
    cout << "d " << d << "\n";
    cout << "probability " << probability << "\n";
    
    // Create price tree
    Triangle<float> stockPrices(steps);
    for (int i=0; i<steps;i++) {
        for (int j=0; j<=i; j++) {
            float val = intialStockPrice*pow(u,j)*pow(d,i-j);
            stockPrices.set(i,j,val); 
        };
    }
    stockPrices.prettyPrint();

    // Solve option price values
    Triangle<float> optionValue(steps);
    // Solve the final "row"
    // For call -> V(j,m) = max(S(j,m)-strike,0)
    for (int i = 0; i<steps;i++) {
        cout << "Price of s j= " << i << " m = " << steps << " is " << stockPrices.get(steps-1,i) << "\n";
        float val = max(stockPrices.get(steps-1,i)-strike,0.0F);
        optionValue.set(steps-1,i,val); 
    }
    optionValue.prettyPrint();

    // backwards traverse and populate prices
    for (int i=steps-2; i>=0;i--) {
        for (int j=0; j<=i; j++) {
            // V(j,m) = e^(-interestRate * deltaTime) * ((1-p)V(j-1,m+1) + p * V(j,m+1)
            float val = pow(EulerConstant,-interestRate*deltaTime) * ((1-probability) * optionValue.getLeftChild(i,j) + probability*optionValue.getRightChild(i,j));
            optionValue.set(i,j,val); 
        };
    }
    optionValue.prettyPrint();
    return optionValue.get(0,0);
}

int main(int argc, char* argv[]) {
    cout << "Hello world\n";
    float price = binomialTreeEuropeOptionPricing(0.2,1,50,1,.05,3,5,1);
    cout << "call option price" << price;
    return 1;
}