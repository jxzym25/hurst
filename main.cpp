#include "hurst.hpp"
#include <iostream>

int main() {
    int n = 200;
    Hurst hurst(240, 20, 4);
    int x = 220;
    for (int i = 0; i < 120; i++){
        hurst.push_back(-log(2));
        hurst.push_back(log(2));
    }
    double h = hurst.compute();
    std::cout << h << std::endl;
    return 0;
}