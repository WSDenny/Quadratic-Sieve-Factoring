#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include "InfInt.h"
#include <omp.h>

using namespace std;

bool isPrime(long int x) {
    if (x == 0 || x == 1) {
        return false;
    } else {
        for (long int i = 2; i <= x / 2; ++i) {
            if (x % i == 0) {
                return false;
                i = x;
            }
        }
        return true;
    }
}

void MatElim(int *A, int rows, int cols) {
    int *temp = new int[cols];
    int piv = 0;
    for (int i = 0; i < cols; i++) {
        if (i < rows) {
            if (A[i * cols + i] == 0) {
                for (int j = i; j < cols; j++) {
                    for (int r = i; r < rows; r++) {
                        if (A[r * cols + j] != 0) {
                            for (int c = 0; c < cols; c++) {
                                temp[c] = A[i * cols + c];
                            }
                            for (int c = 0; c < cols; c++) {
                                A[i * cols + c] = A[r * cols + c];
                            }
                            for (int c = 0; c < cols; c++) {
                                A[r * cols + c] = temp[c];
                            }
                            r = rows;
                            piv = j;
                            j = cols;
                        }
                    }
                }
            }

            if (A[i * cols + piv] != 0) {

                for (int c = cols - 1; c >= 0; c--) {
                    A[i * cols + c] = A[i * cols + c] / A[i * cols + piv];
                }

                for (int c = 0; c < rows; c++) {
                    if (c != i) {
                        for (int k = cols - 1; k >= 0; k--) {
                            A[c * cols + k] = A[c * cols + k] - A[i * cols + k] * A[c * cols + piv];
                        }
                    }
                }
            }
        }
        piv++;
    }
    delete[] temp;
}

void resCheck(InfInt R[], InfInt L[], vector<int> loc, vector<int> primes, int size, InfInt x) {
    InfInt fac1 = 0;
    InfInt fac2 = 0;
    InfInt prodL = 1;
    InfInt prodR = 1;
    InfInt power;

    for (int i = 0; i < size; i++) {
        if (loc[i] > 0) {
            prodL *= L[i];
            for (int j = 0; j < size; j++) {
                if (R[j] > 0) {
                    power = R[j] / 2;
                    for (InfInt k = 0; k < power; k++)
                        prodR *= primes[j];
                }
            }
        }
    }

    fac1 = prodL + prodR;
    fac2 = prodL - prodR;

    if (fac2 < 0)
        fac2 *= -1;

    if (fac1 * fac2 == x) {
        cout << "IT WORKED?!?!?!" << endl;
        cout << "Factor 1: " << fac1 << endl;
        cout << "Factor 2: " << fac2 << endl;
    }
}

int main(void) {
    InfInt x;
    int bound = 100;
    vector<int> primes;
    int thr, tid;

    printf("Input Semiprime number to check for factors: ");
    cin >> x;
    InfInt start = x.intSqrt();

    printf("Input number of threads (or 0 for default): ");
    cin >> thr;

    if (thr <= 0) {
        thr = omp_get_max_threads();
        cout << "Using default number of threads: " << thr << endl;
    }

    for (int i = 2; i <= bound; i++) {
        if (isPrime(i))
            primes.push_back(i);
    }

    int size = primes.size();

    cout << size << endl;
    vector<int> expo;
    InfInt squares[size];

    int *expoUpdt = new int[(size) * size];
    int *expoOrig = new int[(size) * size];

    InfInt iVal = 1;
    InfInt jVal = 1;
    int check = 0;
    InfInt updt;  // Added declaration for updt variable

#pragma omp parallel for num_threads(thr) private(tid,updt)
    for (int i = 0; i < size; i++) {
        squares[i] = 0;
        expoUpdt[i] = 0;
        expoOrig[i] = 0;
    }

    InfInt prod;
    int k = 0;

// Parallelize the loop that generates exponents
#pragma omp parallel for num_threads(thr) private(updt, expo, prod, check, k, tid) 
    for (int loopVar = 0; loopVar < size; loopVar++) {
        tid = omp_get_thread_num();
        k = loopVar;

        printf("Thread %d on %d\n", tid, k);

        if (iVal == 100)
            continue;

        if ((start * start) - x < 0)
            updt = (start * start) % x;
        else
            updt = (start * start) - x;

        for (int j = 0; j < size; j++) {
            expo.push_back(0);

            while (updt % primes[j] == 0 && updt != 1) {
                expo[j] = expo[j] + 1;
                updt = updt / primes[j];
            }
        }

        if (updt == 1) {
            check = 0;
            for (int i = 0; i < size; i++) {
                if (squares[i] == start)
                    check = 1;
            }

            if (check == 0) {
                // Use critical section to update shared arrays safely
        #pragma omp critical
                {
                    for (int j = 0; j < size; j++) {
                        expoUpdt[k * size + j] = expo[j];
                        expoOrig[k * size + j] = expo[j];
                        if (expoUpdt[k * size + j] % 2 == 0) {
                            expoUpdt[k * size + j] = 0;
                        } else {
                            expoUpdt[k * size + j] = 1;
                        }
                    }

                    squares[k] = start;
                }
            }
        }

        if (jVal <= 10000) {
            prod = iVal * x;
            start = (prod.intSqrt()) + jVal;
        } else {
            prod = iVal * x;
            iVal++;
            jVal = 0;
            start = (prod.intSqrt()) + jVal;
        }
        jVal++;
        expo.clear();
    }

    for (int i = 0; i < size; i++)
        cout << primes[i] << " ";

    cout << endl;

    MatElim(expoUpdt, size, size);

    vector<int> freeVar;
    InfInt gcdL[size];
    InfInt gcdR[size];

    int i = 0;
    for (int j = 0; i < size; j++) {
        if (expoUpdt[i * size + j] == 0) {
            freeVar.push_back(1);
        } else {
            freeVar.push_back(0);
            i++;
        }
    }
    cout << endl;

    for (int i = 0; i < size; i++)
        cout << squares[i] << " ";

    cout << endl;

    vector<int> testLoc;

    cout << endl;

    for (int j = 0; j < size; j++) {
        if (freeVar[j] == 1) {
            for (int i = 0; i < size; i++) {
                if (expoUpdt[i * size + j] != 0) {
                    testLoc.push_back(i);
                } else {
                    testLoc.push_back(-1);
                }
            }

            testLoc[j] = j;

            for (int i = 0; i < size; i++) {
                gcdR[i] = 0;
                gcdL[i] = 0;
            }

            for (int i = 0; i < size; i++) {
                if (testLoc[i] >= 0) {
                    int r = testLoc[i];
                    for (int k = 0; k < size; k++) {
                        gcdR[k] += expoOrig[r * size + k];
                    }
                }
            }

            for (int i = 0; i < size; i++) {
                if (testLoc[i] >= 0) {
                    gcdL[i] = squares[i];
                }
            }

            resCheck(gcdR, gcdL, testLoc, primes, size, x);
        }

        testLoc.clear();
    }

    delete[] expoUpdt;
    delete[] expoOrig;

    return 0;
}
