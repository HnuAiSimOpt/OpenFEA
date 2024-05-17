#include <iostream>
#include <armadillo>
using namespace arma;

    int main()
    {
        std::cout << "Hello, World!" << std::endl;
        mat A = {{0.0013, 0.1741, 0.9885, 0.1662, 0.8760},
                 {0.1933, 0.7105, 0.1191, 0.4508, 0.9559},
                 {0.5850, 0.3040, 0.0089, 0.0571, 0.5393},
                 {0.3503, 0.0914, 0.5317, 0.7833, 0.4621},
                 {0.8228, 0.1473, 0.6018, 0.5199, 0.8622}};
        A.print("A: ");
//        mat B = inv(A);
//        B.print("inv(A): ");
//        mat I = A * B;
//        I.print("I: ");

        return 0;
    }



