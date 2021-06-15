//THIS WORKS!!!!

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "Final2Dmat.cpp"
#include "FiveDmat.cpp"

//Instance Variables

//Points
const int tPoints = 4;
const int xPoints = 4;
const int yPoints = 4;
const int zPoints = 4;
const int cPoints = 4;//Cuboidal lattice requirement
const int points = tPoints*xPoints*yPoints*zPoints;//Number of Lattice Points

std::vector<double> avg;
double plaq = 0.0;

//Random Matrix size determiner
const int RowRandM = 8;
const int ColRandM = 8;
//Fixed Table of Random Matrices : These will contain pointers to my SU(3) matrices!
Final2Dmat < std::complex<double>** > TableRandom(RowRandM,ColRandM);
//Table full of coefficients corresponding to my SU(3) matrices I generate randomly.
//While I do not need these so far, I want to keep them
std::complex<double> CoeffTable[RowRandM*ColRandM][9]; 

//To count how many I accept
int accept = 0;
const int links = 4*points;//Number of links
const double beta = 9.0;//Defining Beta - Sort of like temperature, but really lattice spacing

//This is my five dimensional matrix which stores the links I generate and update
//Note that the FiveDMatrix stores data of type pointers to SU() matrices
//The first argument represents direction along which I travel, with all steps being positive 
//So for example, 0 represents t, 1 represents x etc
//The next 4 arguments represent which point in my spacetime lattice I am at
FiveDmat < std::complex<double>** > U(4,tPoints,xPoints,yPoints,zPoints);

//This is really for my debugging purposes; This just prints out a 3x3 matrix
void printmat(std::complex<double>** a){
    for(int i = 0; i<3; i++){
        for (int j = 0; j < 3; j++){
            std::cout<<a[i][j]<<" ";
        }    
        std::cout<<std::endl;
    }
}

//Simple power method for arbitrary complex numbers, with integer powers
std::complex<double> ComplexPow(std::complex<double> z, int p){
    if(p == 0){
        return std::complex<double>(1,0);
    }
    else return z*ComplexPow(z, p-1);

}

//Simple function which returns factorial of an integer number >= 0
int fact(int n){
    if(n == 0){
        return 1;
    }
    return n*fact(n-1);
}

//Returns a random matrix from my Random table fixed at the beginning
std::complex<double>** GetRandomMatrix(){
    int choose1 = rand()%RowRandM;//Random Row Number
    int choose2 = rand()%ColRandM;//Random Column Number
    return TableRandom.get(choose1,choose2);
}

//Deallocates memory from a 3x3 matrix. Use it for memory efficiency
void deallocate(std::complex<double>** a){
    for(int i = 0; i < 3; i++){
        delete[] a[i];//Deallocating memory from a to prevent memory leak
    }
    delete[] a;
}

//Gives Inverse of a SU(3) matrix
//I do this by returning conjugate transpose of a given matrix
//So assume that matrix is Uniatry SU(3) 
std::complex<double>** GetMatrixInverse(std::complex<double>** a){
    std::complex<double>** inv = new std::complex<double>* [3];
        for(int i = 0; i < 3; i++){
            inv[i] = new std::complex<double>[3];
            for(int j = 0; j<3; j++){//Transpose Conjugating
                inv[i][j] = conj(a[j][i]);
            }
        }
        return inv;
    }

//Fills coefficients table and also my Random Table
//Fills it to be within 10 percent of identity!!
//Change it according to your needs to get 60 - 80 % acceptance rate
//Half are matrices and other half are inverse
void FillRandomTable(){
    //Coefficient Matrix, where I write my SU(3) matrix as weighted by some coefficients of the Gellman matrices
    //So this linear combination of coefficients will give me my SU(3) matrix
    for(int i =0; i<RowRandM*ColRandM/2; i++){
        CoeffTable[i][0] = std::complex<double> (0,0);
        CoeffTable[i + RowRandM*ColRandM/2][0] = std::complex<double> (0,0);
        for (int j = 1; j < 9; j++){
            double p1 = (rand()%2)/1.0;//will real part be + or -?
            //Real part of coefficient being randomly generated
            //Imag part will be 0 for unitarity

            //Control percentage away from identity from this statement
            double a = pow(-1.0,p1)*(rand()%3000)/10000.0;
            CoeffTable[i][j] = std::complex<double> (a,0);//Tabulting them in the random table
            CoeffTable[i + RowRandM*ColRandM/2][j] = std::complex<double> (-a,0);//Putting inverse in the random table
        } 
    }

    const std::complex<double> i(0,1);//Defining complex i
    //Defining various Gellman matrices
    std::complex<double> I[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    std::complex<double> SU3_1[3][3] = {{0,1/2.0,0},{1/2.0,0,0},{0,0,0}};
    std::complex<double> SU3_2[3][3] = {{0,-i/2.0,0},{i/2.0,0,0},{0,0,0}};
    std::complex<double> SU3_3[3][3] = {{1/2.0,0,0},{0,-1/2.0,0},{0,0,0}};
    std::complex<double> SU3_4[3][3] = {{0,0,1/2.0},{0,0,0},{1/2.0,0,0}};
    std::complex<double> SU3_5[3][3] = {{0,0,-i/2.0},{0,0,0},{i/2.0,0,0}};
    std::complex<double> SU3_6[3][3] = {{0,0,0},{0,0,1/2.0},{0,1/2.0,0}};
    std::complex<double> SU3_7[3][3] = {{0,0,0},{0,0,-i/2.0},{0,i/2.0,0}};
    std::complex<double> SU3_8[3][3] = {{1/sqrt(12),0,0},{0,1/sqrt(12),0},{0,0,-2/sqrt(12)}};

    //Chose a random matrix to return - Choosing coefficients from the table
    for(int s = 0; s < RowRandM/2; s ++){
        for(int t = 0; t < ColRandM; t++){
            //choose essentially allows us to travel through the coefficient table to calculate all 64 matrices we want
            int choose = t + s*ColRandM;
            //randomMat is the SU(3) matrix we will store
            //randomMatInv is its inverse
            std::complex<double>** randomMat = new std::complex<double>* [3];
            std::complex<double>** randomMatInv = new std::complex<double>* [3];
            std::complex<double>** product = new std::complex<double>* [3];
            //generator is the random generator matrix which generates the SU(3) matrix 
            //combination of Gellmann matrices
            std::complex<double>** generator = new std::complex<double>* [3];
            //temp is the temp matrix in which values calculated are stored  
            std::complex<double>** temp = new std::complex<double>* [3];
            std::complex<double>** copy = new std::complex<double>* [3];
            for (int j = 0; j < 3; j++){
                randomMat[j] = new std::complex<double> [3];
                randomMatInv[j] = new std::complex<double> [3];
                generator[j] = new std::complex<double> [3];
                temp[j] = new std::complex<double> [3];
                product[j] = new std::complex<double> [3];
                copy[j] = new std::complex<double> [3];
                for (int k = 0; k < 3; k++){//Choosing generator matrix
                    generator[j][k] = ( CoeffTable[choose][1]*SU3_1[j][k] + CoeffTable[choose][2]*SU3_2[j][k] + CoeffTable[choose][3]*SU3_3[j][k] + CoeffTable[choose][4]*SU3_4[j][k] + CoeffTable[choose][5]*SU3_5[j][k] + CoeffTable[choose][6]*SU3_6[j][k] + CoeffTable[choose][7]*SU3_7[j][k] + CoeffTable[choose][8]*SU3_8[j][k]);
                    randomMat[j][k] = I[j][k];//Initialising values of random matrix to identity
                    randomMatInv[j][k] = I[j][k];//Initialising values of randomInv matrix to identity
                    temp[j][k] = generator[j][k];
                    copy[j][k] = 0;
                }
            }
            //Starting loop from 1, since 0^th order expansion = identity already done, going to 10^th order
            for(int l = 1; l<10; l++){
                for(int j = 0; j<3; j++){//Creating matrix from generator
                    for (int k = 0; k < 3; k++){
                        //Calculating based on exponential power series
                        randomMat[j][k] += ComplexPow(i, l)*temp[j][k]*std::complex<double>(1.0/fact(l), 0);
                        randomMatInv[j][k] += pow(-1.0,l)*ComplexPow(i, l)*temp[j][k]*std::complex<double>(1.0/fact(l), 0);
                    }
                }
                //Increasing power of generator matrix by 1 in temp variable by using a copy vairable
                for(int j = 0;j<3;j++){
                    for (int k = 0; k < 3; k++){
                        copy[j][k] = 0;
                        for (int l = 0; l < 3; l++){
                            copy[j][k] += temp[j][l]*generator[l][k];
                        }
                    }
                }
                //Rewriting copy vairable back to temp
                for(int j = 0;j<3;j++){
                    for (int k = 0; k < 3; k++){
                        temp[j][k] = copy[j][k];
                    }
                }
            }
            for(int j = 0;j<3;j++){
                for (int k = 0; k < 3; k++){
                    for (int l = 0; l < 3; l++){
                        product[j][k] += randomMat[j][l]*randomMatInv[l][k];
                    }
                }
            }

            //Use it for my own testing purpses to see that the matrices are Unitary and within 10 percent or not
            std::cout<<"\n";
            printmat(randomMat);
            std::cout<<"\n";
            printmat(randomMatInv);
            std::cout<<"\n";
            printmat(product);

            //deallocate memory from arrays
            deallocate(generator);
            deallocate(temp);
            deallocate(copy);
            //Store my randomly generated unitary matrices in the table
            TableRandom.set(randomMat, s,t);
            TableRandom.set(randomMatInv, s+(RowRandM/2), t);
        }
    }
}

//Don't really need this, just feels more natural to call it singleLink though so I have it
std::complex<double>** singleLink(int i, int j, int k, int l, int m){
    return U.get(i,j,k,l,m);
}

//Calculates action associated with a single plane specified by a single link
//Calculates changes to two plaquettes by going around a loop for both the plaquettes
double singlePlane1(std::complex<double>** a, int axes1, int axes2, int j, int k, int l, int m){
    if(axes1 == axes2){
        return 0;
    }
    //Going around first loop
    std::complex<double> trace(0,0);
    int x[] = {0,0,0,0,0};
    x[axes1 + 1] = 1;
    x[axes2 + 1] = 0;
    std::complex<double>** b1 = singleLink(axes2,(j + x[1])%tPoints,(k + x[2])%xPoints,(l + x[3])%yPoints,(m + x[4])%zPoints);
    x[axes1 + 1] = 0;
    x[axes2 + 1] = 1;
    std::complex<double>** c1 = singleLink(axes1,(j + x[1])%tPoints,(k + x[2])%xPoints,(l + x[3])%yPoints,(m + x[4])%zPoints);
    x[axes1 + 1] = 0;
    x[axes2 + 1] = 0;
    std::complex<double>** d1 = singleLink(axes2,(j + x[1])%tPoints,(k + x[2])%xPoints,(l + x[3])%yPoints,(m + x[4])%zPoints);
    std::complex<double>** cInv1 = GetMatrixInverse(c1);
    std::complex<double>** dInv1 = GetMatrixInverse(d1);
    std::complex<double>** product1 = new std::complex<double>*[3];
    for(int i = 0; i<3; i++){
        product1[i] = new std::complex<double>[3];
        for(int j = 0; j<3; j++){
            product1[i][j] = 0;
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    for(int m =0; m<3; m++){
                        product1[i][j] += a[i][k]*b1[k][l]*cInv1[l][m]*dInv1[m][j];
                    }
                }
            }
        }
    }
    trace+= product1[0][0] + product1[1][1] + product1[2][2];
    //memory-efficiency
    deallocate(cInv1);
    deallocate(dInv1);
    deallocate(product1);

    //Going around second loop
    //Note : I need it to be a cuboidal lattice, otherwise I need to check for different rows differently.
    x[axes1 + 1] = 1;
    x[axes2 + 1] = cPoints-1;
    std::complex<double>** b2 = singleLink(axes2,(j + x[1])%tPoints,(k + x[2])%xPoints,(l + x[3])%yPoints,(m + x[4])%zPoints);
    x[axes1 + 1] = 0;
    x[axes2 + 1] = cPoints-1;
    std::complex<double>** c2 = singleLink(axes1,(j + x[1])%tPoints,(k + x[2])%xPoints,(l + x[3])%yPoints,(m + x[4])%zPoints);
    x[axes1 + 1] = 0;
    x[axes2 + 1] = cPoints-1;
    std::complex<double>** d2 = singleLink(axes2,(j + x[1])%tPoints,(k + x[2])%xPoints,(l + x[3])%yPoints,(m + x[4])%zPoints);
    std::complex<double>** bInv2 = GetMatrixInverse(b2);
    std::complex<double>** cInv2 = GetMatrixInverse(c2);
    
    std::complex<double>** product2 = new std::complex<double>*[3];
    for(int i = 0; i<3; i++){
        product2[i] = new std::complex<double>[3];
        for(int j = 0; j<3; j++){
            product2[i][j] = 0;
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    for(int m =0; m<3; m++){
                        product2[i][j] += a[i][k]*bInv2[k][l]*cInv2[l][m]*d2[m][j];
                    }
                }
            }
        }
    }
    trace+= product2[0][0] + product2[1][1] + product2[2][2];

    //memory-efficiency
    deallocate(product2);
    deallocate(bInv2);
    deallocate(cInv2);

    //returning normalised action
    return std::real(trace/3.0);
}

//Calculates action due to all neighbouring planes induced by a given matrix at a single site
double NormalizedRealTrace(std::complex<double>** a, int i, int j, int k, int l, int m){
    double trace = 0.0;
    for(int x = 0; x<4; x++){
        trace+= singlePlane1(a,i,x,j,k,l,m);//Adding over all planes containing that edge at the points
    }
    return trace;
}

//This method is intended to initially configure the entire lattice links with random matrices
void SetLinks(){//Starting matrix
    std::complex<double> Start[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    for(int i = 0; i<4; i++){
        for(int j = 0; j<tPoints; j++){
            for(int k = 0; k<xPoints; k++){
                for(int l = 0; l<yPoints; l++){
                    for(int m = 0; m<zPoints; m++){
                        std::complex<double>** init = new std::complex<double>* [3];
                        for(int i = 0; i<3; i++){
                            init[i] = new std::complex<double>[3];
                            for (int j = 0; j < 3; j++){
                                init[i][j] = Start[i][j];
                            }    
                        }
                        U.set(init,i,j,k,l,m);
                    }
                }
            }
        }
    }
}

//This method goes over and updates the lattice link using Monte Carlo methods
void MonteCarlo(){
    for(int i = 0; i<4; i++){
        for(int j = 0; j<tPoints; j++){
            for(int k = 0; k<xPoints; k++){
                for(int l = 0; l<yPoints; l++){
                    for(int m = 0; m<zPoints; m++){
                        //Generate a random matrix for the MonetCarlo procedure
                        std::complex<double>** V = GetRandomMatrix();
                        std::complex<double>** currentU = U.get(i,j,k,l,m);

                        //Multiply random matrix V with current matrix to generate the Monte Carlo configuration
                        std::complex<double>** trialU = new std::complex<double>* [3];
                        for(int p = 0; p<3; p++){
                            trialU[p] = new std::complex<double>[3];
                            for (int q = 0; q < 3; q++){
                                for(int r = 0; r < 3; r++){
                                    trialU[p][q] += V[p][r]*currentU[r][q];
                                }    
                            }   
                        }
                        //Calculate action (weights)
                        double act_old = 6.0 - 1.0*NormalizedRealTrace(currentU, i,j,k,l,m);
                        double act_new = 6.0 - 1.0*NormalizedRealTrace(trialU,i,j,k,l,m) ;
                        double delta_act = act_new - act_old;

                        //Updating links according to the probabilities
                        double prob = (rand()%10000)/10000.0;
                        double compare = exp(-1.0*beta*(delta_act));
                        if(prob < compare){
                            U.set(trialU,i,j,k,l,m);
                            accept++;
                            deallocate(currentU);
                        }
                        else deallocate(trialU);                   
                    }
                }
            }
        }
    }
}

//This gives back the average Plaquette value
//The method right now is very inefficient, please change later if needed in complex calculations.
void avgPlaq(){
    plaq = 0.0;
    double avgaction = 0.0;
    for(int i = 0;i<4;i++){
        for(int j = 0;j<tPoints;j++){
            for(int k = 0;k<xPoints;k++){
                for(int l = 0;l<yPoints;l++){
                    for(int m = 0; m<zPoints; m++){//Summing over all 6 plquettes for each link
                        avgaction += NormalizedRealTrace(singleLink(i,j,k,l,m), i,j,k,l,m);
                    }
                }
            }
        }
    }
    //Realizing we have over-counted by 6 times so additional division by 6
    plaq = avgaction/(links*6.0);
    std::cout<<avgaction/(links*6.0)<<std::endl;
}

int main(){
    srand( time (NULL) );
    //First Fill table, and then hold that table constant for all loops
    FillRandomTable();
    //Set Links - In this Program it is a cold start
    SetLinks();
    for(int i = 0; i< 5000; i++){
        accept = 0;
        std::cout<<i<<std::endl;
        //Keep doing MonteCarlo
        MonteCarlo();
        std::cout<<"Acceptance Percentage : "<<accept/(links*1.0)<<std::endl;
        std::cout<<"Average plaquette : ";
        avgPlaq();
        avg.push_back(1 - plaq);
    }
    double tot = 0.0;
    double sd = 0.0;
    for(int i = 1000;i <avg.size(); i++){
        tot+= avg[i];
    }
    double mean = tot/(avg.size()-1000);
    for(int i = 1000;i<avg.size(); i++){
        sd+= (avg[i] - mean)*(avg[i] - mean);
    }
    sd = sqrt(sd)/sqrt(avg.size()-1000);
    std::cout<<"Mean is "<<mean<<"\n";
    std::cout<<"Standard Deviation is "<<sd<<"\n";
}
