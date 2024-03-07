#include <cmath>
#include <iostream>
#include <cstdlib>
#include <random>
#include <vector>
#include <chrono>

#include "InfInt.h"



using namespace std;

bool isPrime(long int x) {
    if(x == 0 || x == 1) {
        return false;
    } else {
        for(long int i = 2; i <= x/2; ++i) {
            if(x % i == 0) {
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
    for(int i = 0; i < rows; i++) {
        if(i < rows) {
            if(A[i * cols + piv] == 0) {
                for(int j = i; j < cols; j++) {
                    for(int r = i; r < rows; r++) {
                        if(A[r * cols + j] != 0) {
                            for(int c = 0; c < cols; c++) {
                                temp[c] = A[i * cols + c];
                            }
                            for(int c = 0; c < cols; c++) {
                                A[i * cols + c] = A[r * cols + c];
                            }
                            for(int c = 0; c < cols; c++) {
                                A[r * cols + c] = temp[c];
                            }
                            r = rows;
                            piv = j;
                            j = cols;
                        }
                    }
                }
            }
            
            
            if(A[i * cols + piv] != 0) {
                
                for(int c = cols - 1; c >= 0; c--){
                    A[i * cols + c] = A[i * cols + c] / A[i * cols + piv];
                }
                
                for(int c = 0; c < rows; c++) {
                    if(c != i) {
                        for(int k = cols - 1; k >= 0; k--) {
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
    
bool resCheck(int R[], InfInt L[], int loc[], vector<int> primes, int rows, int cols, InfInt x) {
    InfInt fac1 = 0;
    InfInt fac2 = 0;
    InfInt prodL = 1;
    InfInt prodR = 1;
    InfInt power;
    InfInt *pVal = new InfInt[rows];
    
    for(int i = 0; i < rows; i++)
        pVal[i] = 0;
    
    for(int i = 0; i < cols; i++) {
        if(loc[i] != -1 ) {
            prodL *= L[i];
            for(int j = 0; j < rows; j++) {
                if(R[j * cols + i] > 0) {
                    pVal[j] += R[j * cols + i];

                }
            }
        }
    }
    
    for(int i = 0; i < rows; i++) {
        if(pVal[i] > 0) {
            power = pVal[i] / 2;
            for(InfInt k = 0; k < power; k++)
                prodR *= primes[i];
        }
    }
        
    
    fac1 = prodL + prodR;
    
    fac2 = prodL - prodR;
    
    if(fac2 < 0)
        fac2 *= -1;
    
    delete[] pVal;
    
    if(fac1 * fac2 == x) {
        cout << "IT WORKED?!?!?!" << endl;
        cout << "Factor 1: " << fac1 << endl;
        cout << "Factor 2: " << fac2 << endl;
        return true;
    } else {
        return false;
    }
}

int main(void) {
    
    string tx;
    InfInt x;
    int bound;
    InfInt updt;
    vector<int> primes;
    
    std::chrono::time_point<std::chrono::system_clock> begin, end;

    
    cout << "Enter N: ";
    cin >> tx;
    
    x = tx;
    InfInt start = x.intSqrt();
    
    cout << "Enter B: ";
    cin >> bound;
        
    
    for(int i = 2; i <= bound; i++) {
        if(isPrime(i))
            primes.push_back(i);
    }
    
    int size = primes.size();
    
    vector <int> expo;
    //     vector<int> loc;
    
    int k = 0;
    
    int rows = size;
    int var;
    
    cout << "More relations to find > primes in B (0 recommended): ";
    cin >> var;
        
    int cols = size + var;
    
    InfInt iLim;
    InfInt jLim;
    
    cout << "Enter i limit: ";
    cin >> iLim;
    cout << "Enter j limit: ";
    cin >> jLim;
    
    int *expoUpdt = new int[rows * cols];
    int *expoOrig = new int[rows * cols];
    InfInt squares[cols];   

    
    InfInt iVal = 1;
    InfInt ii = 1;
    InfInt jVal = 0;
    int check = 0;
    
    for(int i = 0; i < rows * cols; i++) {
        expoUpdt[i] = 0;
        expoOrig[i] = 0;
    }
    
    for(int i = 0; i < cols; i++)
        squares[i] = 0;
    
    InfInt prod;
    
    begin = std::chrono::system_clock::now();
        
    while(k < cols && iVal < iLim) {
        //cout << start << endl;
        //cout << k << " ";
        if((start * start) - x < 0)
            updt = (start * start) % x;
        else
            updt = (start * start) - x;
       //    if (updt != 1) {
        
        for(int j = 0; j < rows; j++) {
            expo.push_back(0);
            
            while(updt % primes[j] == 0 && updt != 1) {
                //cout << updt << endl;
                expo[j] = expo[j] + 1;
                updt = updt / primes[j];
            }
        }
        
        if(updt == 1) {
            check = 0;
            for(int i = 0; i < rows; i++) {
                if(squares[i] == start) 
                    check = 1;
            }
                                
            
            if(check == 0) {
                for(int j = 0; j < rows; j++) {
                    expoUpdt[j * cols + k] = expo[j];
                    expoOrig[j * cols + k] = expo[j];
                    if(expoUpdt[j * cols + k] % 2 == 0) {
                        expoUpdt[j * cols + k] = 0;
                    } else {
                        expoUpdt[j * cols + k] = 1;
                    }
                }
                
                squares[k] = start;
                k++;
            }
            }

            
        if(jVal <= jLim) {
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
    
//     for(int i = 0; i < rows; i++) {
//         for(int j = 0; j < cols; j++) {
//             cout << expoOrig[i * cols + j] << " ";
//         }
//         cout << endl;
//     }
    
    MatElim(expoUpdt, rows, cols);
    
//     for(int i = 0; i < rows; i++) {
//         for(int j = 0; j < cols; j++) {
//             cout << expoUpdt[i * cols + j] << " ";
//         }
//         cout << endl;
//     }
    
    
    vector<int> freeVar;

    int i = 0;
        for(int j = 0; j < cols; j++) {
            if(expoUpdt[i * cols + j] == 0) {
                freeVar.push_back(1);
            } else {
                freeVar.push_back(0);
                i++;
            }
        }
        cout << endl;
    
    int testLoc[cols];
    bool display;
    int d = 0;
    
    for(int i = 0; i < cols; i++)
        testLoc[i] = -1;
    
    cout << endl;
    
    for(int j = 0; j < cols; j++) {
        
        for(int i = 0; i < cols; i++)
            testLoc[i] = -1;
        
        if(freeVar[j] == 1) {
            for(int i = 0; i < rows; i++) {
                if(expoUpdt[i * cols + j] != 0) {
                    testLoc[i] = i;
                }
            }
            
            testLoc[j] = j;
            
//             for(int k = 0; k < cols; k++)
//                 cout << testLoc[k] << " ";
//             
//             cout << endl;
                        
           display = resCheck(expoOrig, squares, testLoc, primes, rows, cols, x);
           
           if(display == true) {
               d++;
           }else{
               d = d;
           }
               
        }
    }

        end = std::chrono::system_clock::now();
        
        std::chrono::duration<double> elapsed_seconds = end - begin;
        
        if(d == 0)
            cout << "Could not find prime factors within bounds" << endl;
        
        cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    delete[] expoUpdt;
    delete[] expoOrig;

return 0;
}

