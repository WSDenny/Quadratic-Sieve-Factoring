#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>

#include <mpi.h>

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
    for(int i = 0; i < cols; i++) {
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
    
bool resCheck(int R[], InfInt *L, int loc[], int *primes, int rows, int cols, InfInt x) {
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
    InfInt x, start, updt, prod, iVal, iLim, jVal, jValScale, jLim;
    string tempx;
    int s, rows, cols, il, jl;
    int commsz, rank, size, padsize, bound;
    double begin, end;
    vector<int> primesVec;
    int *expoUpdt = NULL;
    int *expoOrig = NULL;
    int *primes = NULL, *expo = NULL, *locsquares2 = NULL, *squares2 = NULL;
    string *locsquares = NULL;
    InfInt *squares = NULL;
        
    MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &commsz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if(rank == 0) {
        cout << "Enter N: ";
        cin >> tempx;
    
        s = tempx.size();
    
        cout << "Enter B: ";
        cin >> bound;
                    
        for(int p = 2; p < bound; p++) {
            if(isPrime(p))
                primesVec.push_back(p);
        }
        
        size = primesVec.size();
    }
    
        MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
        primes = new int[size];
        
    if(rank == 0) {
        for(int i = 0; i < size; i++) {
            primes[i] = primesVec[i];
        }
        
        if (size % commsz == 0)
			padsize = size;
		else {
			int div = size / commsz;
			padsize = (div + 1) * (commsz * (size / 2));
		}
		
		rows = size;

        squares2 = new int[padsize];
        
        cout << "Enter i limit: ";
        cin >> il;
        cout << "Enter j limit: ";
        cin >> jl;
        
//         cout << padsize << endl;
//        cout << padsize / commsz << endl;
//         cout << cols << endl;
//         cout << rows << endl;
//         cout << il << endl;
//         cout << jl << endl;
        
        begin = MPI_Wtime();
    }
        
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    char brod[s];
    if(rank == 0) {
        for(int i = 0; i <= s; i++)
            brod[i] = tempx[i];
    }
    
    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&padsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(primes, rows, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&il, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&jl, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(brod, s, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    if(rank == 0)
        begin = MPI_Wtime();
    
    x = brod;
    iVal = 1;
    iLim = il;
    jValScale = jl / commsz;
    jVal = jValScale * rank;
    jLim = jValScale * (rank + 1);
    prod = iVal * x;
    start = (prod.intSqrt()) + jVal;
    string save = start.toString();
    int k = 0;
    int check = 0;
//     cout << " " << jVal << " " << jLim << " " << rank << endl;
    locsquares = new string[padsize / commsz];
    locsquares2 = new int[padsize / commsz];

    while(k < padsize / commsz && iVal < iLim) {
        if((start * start) - x < 0)
            updt = (start * start) % x;
        else
            updt = (start * start) - x;
        
        for(int j = 0; j < rows; j++) {            
            while(updt % primes[j] == 0 && updt != 1) {
                updt = updt / primes[j];
            }
        }
        
        if(updt == 1) {
            check = 0;
            for(int i = 0; i < rows; i++) {
                if(locsquares[i] == start.toString()) 
                    check = 1;
            }
                                
            if(check == 0) {
                
                locsquares[k] = start.toString();
                locsquares2[k] = stoi(locsquares[k]);
                k++;
            }
            }

            
        if(jVal < jLim) {
            prod = iVal * x;
            start = (prod.intSqrt()) + jVal;

        } else {
            prod = iVal * x;
            iVal++;
            jVal = jValScale * rank;
            start = (prod.intSqrt()) + jVal;
        }
        jVal++;           
    }
        
//     for(int i = 0; i < padsize / commsz; i++)
//         cout << locsquares[i] << " ";
    //cout << endl;
    
    MPI_Gather(locsquares2, padsize / commsz, MPI_INT, squares2, padsize / commsz,
		MPI_INT, 0, MPI_COMM_WORLD);
        
    if(rank == 0) {
    int count = 0;
    int save2 = stoi(save);
    
    for(int i = 0; i < padsize; i++) {
        if(squares2[i] < save2)
            count++;
    }
            
        cols = padsize - count;
        cout << padsize << endl;
        squares = new InfInt[cols];
        
        count = 0;
        
    for(int i = 0; i < padsize; i++) {
        if(squares2[i] >= save2) {
            squares[count] = squares2[i];
            count++;
        }
    }
        expoUpdt = new int[rows * cols];
        expoOrig = new int[rows * cols];
        expo = new int[rows];
        
        
        for(int i = 0; i < cols; i++) {
                updt = squares[i];
                
                for(int j = 0; j < rows; j++) {
                    expo[j] = 0;
                
                    while(updt % primes[j] == 0 && updt != 1) {
                        expo[j] = expo[j] + 1;
                        updt = updt / primes[j];
                    }
                }
                
                for(int j = 0; j < rows; j++) {
                    expoUpdt[j * cols + i] = expo[j];
                    expoOrig[j * cols + i] = expo[j];
                    if(expoUpdt[j * cols + i] % 2 == 0) {
                        expoUpdt[j * cols + i] = 0;
                    } else {
                        expoUpdt[j * cols + i] = 1;
                    }
                }
        }
        
//         for(int i = 0; i < rows; i++) {
//             for(int j = 0; j < cols; j++) {
//                 cout << expoOrig[i * cols + j] << " ";
//             }
//             cout << endl;
//         }


        MatElim(expoUpdt, rows, cols);
        
//         for(int i = 0; i < rows; i++) {
//             for(int j = 0; j < cols; j++) {
//                 cout << expoOrig[i * cols + j] << " ";
//             }
//             cout << endl;
//         }       
        
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
        
        if(d == 0)
            cout << "Could not find prime factors within bounds" << endl;
        
        end = MPI_Wtime();
        
        cout << "elapsed time: " << end - begin << "s\n";
        
    }
    

    delete[] expoUpdt;
    delete[] expoOrig;
    delete[] locsquares;
    delete[] locsquares2;
    delete[] primes;
    delete[] squares;
    delete[] squares2;
    delete[] expo;
    MPI_Finalize();

return 0;
}

