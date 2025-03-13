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
    std::vector<T> data;  // Flattened 1D array

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

    // Pretty print the triangle
    void prettyPrint() const {
        int index = 0;
        int maxWidth = 4 * levels; // Controls spacing for alignment
        for (int L = 0; L < levels; L++) {
            // Print leading spaces for centering
            int numElements = L + 1;
            int padding = (maxWidth - numElements * 3) / 2;
            std::cout << std::string(padding, ' ');

            // Print the row elements
            for (int P = 0; P <= L; P++) {
                std::cout << std::setw(6) << std::fixed << std::setprecision(2) 
                          << data[index++] << " ";
            }
            std::cout << "\n";
        }
    }
};

float binomialTreeEuropeOptionPricing(float sigma, float direction, int steps, int time, float interestRate, float strike, float intialStockPrice, int putOrCall) {
    const double EulerConstant = std::exp(1.0);
    // Calculate p, u and d
    float deltaTime = (float)time/steps;
    float beta = 0.5*(direction*pow(EulerConstant,-interestRate*deltaTime)+pow(EulerConstant,(interestRate+pow(sigma,2))*deltaTime));
    float u = beta + sqrt(pow(beta,2)-direction);
    float d = direction/u;
    float probability = (pow(EulerConstant,interestRate*deltaTime)-d)/(u-d);

    cout << "Delta time " << deltaTime << "\n";
    cout << "beta " << beta << "\n";
    cout << "u " << u << "\n";
    cout << "d " << d << "\n";
    cout << "probability " << probability << "\n";

    Triangle<float> stockPrices(steps);  // Triangle with 4 levels
    for (int i=0; i<steps;i++) {
        for (int j=0; j<=i; j++) {
            float val = intialStockPrice*pow(u,j)*pow(d,i-j);
            stockPrices.set(i,j,val); 
        };
    }
    stockPrices.prettyPrint();
    
    return -1;
}

int main(int argc, char* argv[]) {
    cout << "Hello world\n";
    binomialTreeEuropeOptionPricing(0.2,1,256,1,.05,0,5,1);
    return 1;
}