//Hit Method 
//THIS WORKS!!!
//Calculating RxT Wilson Loops
// In this problem due to symmetry over all directions, note that 
// I am actually summing over all possible 2x1 loops regardless of the axes

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "Final2Dmat.cpp"
#include "FiveDmat.cpp"
#include "fstream"

//myfile for output purposes
std::ofstream myfile;
//Instance Variables

//Points
const int tPoints = 6;
const int xPoints = 6;
const int yPoints = 6;
const int zPoints = 6;
const int cPoints = 6;//Cuboidal lattice requirement
const int points = tPoints*xPoints*yPoints*zPoints;//Number of Lattice Points

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
double beta = 9.0;//Defining Beta - Sort of like temperature, but really lattice spacing

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
//Please deallocate this once used
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
//Fills it to be within 10 (easily changeable) percent of identity!!
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
            double a = pow(-1.0,p1)*(rand()%2000)/10000.0;
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
            //Just for checking that matrix multiplied by inverse does give me identity right
            for(int j = 0;j < 3;j++){
                for (int k = 0; k < 3; k++){
                    for (int l = 0; l < 3; l++){
                        product[j][k] += randomMat[j][l]*randomMatInv[l][k];
                    }
                }
            }

            //deallocate memory from arrays
            deallocate(generator);
            deallocate(temp);
            deallocate(copy);
            //Store my randomly generated unitary matrices in the table
            TableRandom.set(randomMat, s,t);
            TableRandom.set(randomMatInv, s+(RowRandM/2), t);
            deallocate(product);
        }
    }
}

//Don't really need this, just feels more natural to call it singleLink though so I have it
std::complex<double>** singleLink(int i, int j, int k, int l, int m){
    return U.get(i,j,k,l,m);
}

//Returns product matrix asocciated with a single plaquette in a single plane (out of 3 possible planes)
std::complex<double>** singlePlane(int axes1, int axes2, int j, int k, int l, int m){
    if(axes1 == axes2){
        return nullptr;
    }
    //Going around first loop
    std::complex<double> trace(0,0);
    int x[] = {0,0,0,0,0};//Changing x vector as shown helps me go around a plaquette
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
    std::complex<double>** product = new std::complex<double>*[3];

    //Multiplying 3 matrices along the first plaquette I consider
    for(int i = 0; i<3; i++){
        product[i] = new std::complex<double>[3];
        for(int j = 0; j<3; j++){
            product[i][j] = 0;
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    product[i][j] += b1[i][k]*cInv1[k][l]*dInv1[l][j];
                }
            }
        }
    }
    
    //memory-efficiency
    deallocate(cInv1);
    deallocate(dInv1);

    //Going around second loop
    //Note : I need it to be a cuboidal lattice, otherwise I need to check for different rows differently.
    //Have done so in Problem3 - Phase Transition: Check it out if required!
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
    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    product[i][j] += bInv2[i][k]*cInv2[k][l]*d2[l][j];
                }
            }
        }
    }
    //memory-efficiency
    deallocate(bInv2);
    deallocate(cInv2);

    //returning normalised action
    return product;
}

//returns all matrices in all plaquette associated with a link (6 plaquettes in 3 planes)
std::complex<double>** singlePoint(int i, int j, int k, int l, int m){
    std::complex<double>** mat = new std::complex<double>*[3];
    for (int x = 0; x < 3; x++){
        mat[x] = new std::complex<double>[3];
    }
    for (int x = 0; x < 4; x++){
        std::complex<double>** planar = singlePlane(i,x,j,k,l,m);
        if(planar == nullptr){//Happens when both axes are similar
            continue;
        }
        for(int a = 0; a < 3; a++){
            for (int b = 0; b < 3; b++){
                mat[a][b]+= planar[a][b];//Summing all plaquettes over all planes
            }    
        }
        deallocate(planar);
    }
    return mat;    
}

//Multiply the current link with the sum of all matrices in associated plaquettes
//Take real trace normalized and return that
//First argument is current link!
double NormalizedRealTrace(std::complex<double>** mat1, std::complex<double>** mat2){
    std::complex<double>** Plaquette = new std::complex<double>*[3];
    for(int a = 0; a<3; a++){
        Plaquette[a] = new std::complex<double>[3];
        for(int b = 0; b<3; b++){
            for(int c = 0; c<3; c++){
                Plaquette[a][b]+=mat1[a][c]*mat2[c][b];//Multiplying matrices 
            }
        }
    }
    std::complex<double> sum = Plaquette[0][0] + Plaquette[1][1] + Plaquette[2][2];
    double trace = real(sum)/3.0;
    deallocate(Plaquette);
    return trace;
}

//This method is intended to initially configure the entire lattice links with random matrices
void SetLinks(){//Starting matrix
    std::complex<double> Start[3][3] = {{1,0,0},{0,1,0},{0,0,1}};//Starting matrix - Cold
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
                        U.set(init,i,j,k,l,m);//Cold start!
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
                        std::complex<double>** mat = singlePoint(i,j,k,l,m);
                        //Hit procedure
                        for(int hit = 0; hit<10; hit++){
                            std::complex<double>** V = GetRandomMatrix();
                            //Multiply random matrix V with current matrix to generate the Monte Carlo configuration
                            std::complex<double>** trialU = new std::complex<double>* [3];
                            for(int p = 0; p<3; p++){
                                trialU[p] = new std::complex<double>[3];
                                for (int q = 0; q < 3; q++){
                                    for(int r = 0; r < 3; r++){
                                        trialU[p][q] += V[p][r]*U.get(i,j,k,l,m)[r][q];
                                    }    
                                }   
                            }
                            //Calculate action (weights) - They are equivalent to Wilson loops
                            //The 6 is because of 6 loops in a plane, and identity matrix for each loop
                            double act_old = 6.0 - 1.0*NormalizedRealTrace(U.get(i,j,k,l,m), mat);
                            double act_new = 6.0 - 1.0*NormalizedRealTrace(trialU,mat) ;
                            double delta_act = act_new - act_old;

                            //Updating links according to the probabilities
                            double prob = (rand()%10000)/10000.0;
                            double compare = exp(-1.0*beta*(delta_act));
                            if(prob < compare){
                                deallocate(U.get(i,j,k,l,m));
                                U.set(trialU,i,j,k,l,m);
                                accept++;
                            }
                            //Memory-efficiency
                            else deallocate(trialU);
                        }                   
                    }
                }
            }
        }
    }
}

//This method takes in a vector of SU(3) matrices, and multiply them in order starting from 0 index
std::complex<double>** Prod(std::vector<std::complex<double>**> a){
    std::complex<double>** temp = new std::complex<double>*[3];
    std::complex<double>** prod = new std::complex<double>*[3];
    for(int i = 0;i<3;i++){
        temp[i] = new std::complex<double> [3];
        prod[i] = new std::complex<double> [3];
    }
    for(int j = 0; j<3; j++){
        for(int k = 0; k<3; k++){
            prod[j][k] = a[0][j][k];
        } 
    }
    for(int x = 1; x<a.size(); x++){
        for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                temp[i][j] = 0;//Reinitializing temp to 0 after each iteration
                for(int k = 0; k<3; k++){
                    temp[i][j] += prod[i][k]*a[x][k][j];//Multiplying SU(3)
                } 
            }
        }
        for(int j = 0; j<3; j++){
            for(int k = 0; k<3; k++){
                prod[j][k] = temp[j][k];//Stroing them in product
            } 
        }
    }
    deallocate(temp);
    return prod;
    }

//Calculates a single RxT 
double singleRxT(int axes1, int axes2, int j, int k, int l, int m, int r , int t){
    if(axes1 == axes2){
        return 0.0;
    }
    std::vector<std::complex<double>**> t1Links;
    std::vector<std::complex<double>**> t2Links;
    std::vector<std::complex<double>**> r1Links;
    std::vector<std::complex<double>**> r2Links;

    //These go along the Wilson loop in the required procedure
    int x[] = {0,0,0,0};

    for (int i = 0; i<t; i++){
        t1Links.push_back(U.get(axes1,(j+x[0])%cPoints,(k+x[1])%cPoints,(l+x[2])%cPoints,(m+x[3])%cPoints));
        x[axes1]+=1;
    }
    
    for (int i = 0; i<r; i++){
        r1Links.push_back(U.get(axes2, (j+x[0])%cPoints, (k+x[1])%cPoints, (l+x[2])%cPoints, (m+x[3])%cPoints));
        x[axes2]+=1;
    }

    //You were earlier doing this is a Gauge - Uninvariant way!!
    for (int i = 0; i<t; i++){
        x[axes1]-=1;
        t2Links.push_back(GetMatrixInverse(U.get(axes1, (j+x[0])%cPoints, (k+x[1])%cPoints, (l+x[2])%cPoints,(m+x[3])%cPoints)));
    }
    for (int i = 0; i<r; i++){
        x[axes2]-=1;
        r2Links.push_back(GetMatrixInverse(U.get(axes2, (j+x[0])%cPoints, (k+x[1])%cPoints, (l+x[2])%cPoints,(m+x[3])%cPoints)));
    }
    std::complex<double>** prod = new std::complex<double>* [3];
    //Product along each arm
    std::complex<double>** prod1= Prod(t1Links);
    std::complex<double>** prod2= Prod(r1Links);
    std::complex<double>** prod3= Prod(t2Links);
    std::complex<double>** prod4= Prod(r2Links);
    for(int i = 0; i<3; i++){
        prod[i] = new std::complex<double>[3];
        for(int j = 0; j<3; j++){
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    for(int m =0; m<3; m++){
                        prod[i][j] += prod1[i][k]*prod2[k][l]*prod3[l][m]*prod4[m][j];//Product of all arms
                    }
                }
            }
        }
    }
    double trace = ((real(prod[0][0])) + (real(prod[1][1])) + (real(prod[2][2])))/3.0;
    deallocate(prod); deallocate(prod1); deallocate(prod2); deallocate(prod3);deallocate(prod4);
    //Memory efficiency
    for (int i = 0; i<t; i++){
        deallocate(t2Links[i]);
    }
    for (int i = 0; i<r; i++){
        deallocate(r2Links[i]);
    }
    return trace;
}

//Calculate an RxT Wilson loop
double RxT(int r, int t){
    double W = 0.0;
    for(int i = 0; i<4; i++){
        for(int j = 0; j<tPoints; j++){
            for (int k = 0; k < xPoints; k++){
                for (int l = 0; l < yPoints; l++){
                    for (int m = 0; m < zPoints; m++){
                        for(int p = 0; p<4; p++){
                            W+= singleRxT(i,p,j,k,l,m,r,t);//Summing over all single Wilson Loops
                        }   
                    }
                }   
            }
        }
    }
    return W/(points*2*6.0);
}


int main(){
    srand( time (NULL) );
    myfile.open("Values2x2.txt", std::ios::out | std::ios::trunc);
    for(int k = 1; k < 9; k++){
        beta = k*0.5;
        //First Fill table, and then hold that table constant for all loops
        FillRandomTable();
        //Set Links - In this Program it is a cold start
        SetLinks();
        for(int i = 0; i< 500; i++){
            accept = 0;
            std::cout<<i<<"\n";
            //Keep doing MonteCarlo
            MonteCarlo();
        }
        std::vector<double> avg;
        for(int i = 500; i< 2500; i++){
            accept = 0;
            std::cout<<i<<"\n";
            //Keep doing MonteCarlo
            MonteCarlo();
            double wLoop = RxT(2,2);
            avg.push_back(wLoop);
        }

        double tot = 0.0;
        for(int i = 0; i <avg.size(); i++){
            tot+= avg[i];
        }
        double mean = tot/(avg.size());
        myfile<<"Results for beta = "<<beta<<"\n";
        myfile<<"Mean is : "<<mean<<"\n";
        double sd = 0.0;
        std::vector<double> bin;
        for(int i = 500;i<avg.size(); i+=20){
            double temp = 0.0;
            for(int j = 0; j<20; j++){
                temp+= avg[i+j];
            }
            bin.push_back(temp/20.0);
        }
        for (int i = 0; i < bin.size(); i++){
            sd += (bin[i] - mean)*(bin[i] - mean);
        }
        sd = sqrt(sd)/sqrt(bin.size());
        myfile<<"Standard Deviation is : "<<sd<<"\n";
        myfile<<"\n";
    }
    myfile.close();
}
