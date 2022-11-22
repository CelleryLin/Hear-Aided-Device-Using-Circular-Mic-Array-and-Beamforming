#include <math.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include  <ctime>
#include <iomanip>
#include<complex.h>
#include<iterator>
#include<random>
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Eigenvalues"

using namespace Eigen;
using namespace std;

# define PI 3.14159265358979323846

int main() {

    int N1 = 200, lambda = 150, d = lambda / 2, snr = 10, M = 6, P = 2;

    double Pmusic = 0, theata = 0, Max = 0, doa[2] = { 20 * PI / 180,60 * PI / 180 }, w[2][1] = {
        {PI / 4},
        {PI / 3}
    };

    complex<double> t = { 0.0,-1.0 }, temp = {}, temp2 = {}, temp3 = {}, T1 = { 1,1 };
    typedef Matrix<complex<double>, Dynamic, Dynamic> MatrixXcd;
    MatrixXcd A(M,M);

    //Create random numbers
    const double mean = 0.0;
    default_random_engine generator;
    normal_distribution<double> norm(mean, 0.1);
    complex<double>  randnumber = { norm(generator),norm(generator) };


    //memory space
    complex<double> DT[M][P] = {}, xx[P][N1] = {}, x[M][N1] = {}, xt[N1][N1] = {}, R[M][M] = {}, SS[1][M] = {}, SST[M][M] = {}, NN[M][4] = {}, NNT[M][M] = {}, PP[1][1] = {};


    //calculate matrix

    for (int m = 0;m < M;m++) {
        for (int p = 0;p < P;p++) {
            DT[m][p] = { cos(2 * PI * d * sin(doa[p]) * m / lambda),  sin((-2) * PI * d * sin(doa[p]) * m / lambda) };
            DT[m][p] = conj(DT[m][p]);
        }       
    }

    //DT matrix(6*2)
    for (int p = 0;p < P;p++) {
        for (int n = 0;n < N1;n++) {
            xx[p][n] = { 2 * cos(w[p][0] * (n + 1)),2 * sin(w[p][0] * (n + 1)) };
            
        }       
    }

    for (int m = 0;m < M;m++) {
        for (int n = 0;n < N1;n++) {
            for (int p = 0;p < P;p++) {
                x[m][n] += DT[m][p] * xx[p][n];
               // x[m][n] = x[m][n] + randnumber;
                xt[m][n] = x[m][n];
            }            
        }     
    }
    // x matrix(6*200)


    for (int m = 0;m < M;m++) {
        for (int n = m;n < N1;n++) {
                temp = conj(xt[m][n]);
                xt[m][n] = conj(xt[n][m]);
                xt[n][m] = temp;  
        }       
    }

    //xt matrix(200*6)
    for (int i = 0;i < 6;i++) {
        for (int m = 0;m < M;m++) {
            for (int n = 0;n < N1;n++) {
                R[i][m] += x[i][n] * xt[n][m];
            }
            //cout << R[i][m] << "\t";
        }
        //cout << endl;
    }

    // R matrix(6*6)

    

    //eigenvalues and eigenvectors
    for (int i = 0;i < M;i++) {
        for (int m = 0; m < M;m++){
            A(i, m) = R[i][m];
        }
    }

    ComplexEigenSolver<MatrixXcd> es(A);
    MatrixXcd D = es.eigenvalues();
    MatrixXcd V = es.eigenvectors();
    for (int i = 0;i < M;i++) {
        for (int m = 0;m < M - P;m++) {
            NN[i][m] = V(i, m);
            NNT[i][m] = NN[i][m];
        }      
    }

    //NN matrix(6*4)
    for (int i = 0;i < M;i++) {
        for (int m = i;m < M - P;m++) {
            temp2 = conj(NNT[i][m]);
            NNT[i][m] = conj(NNT[m][i]);
            NNT[m][i] = temp2;
        }       
    }

    //NNT matrix(4*6)




    //peak search
    for (double ii = -90.0;ii <= 90.0;ii = ii + 0.5) {
        for (int jj = 0;jj < M;jj++) {
            SS[0][jj] = { cos(jj * 2 * PI * d * sin(ii * PI / 180) / lambda),sin(jj * (-2) * PI * d * sin(ii * PI / 180) / lambda) };
            SST[jj][0] = conj(SS[0][jj]);
        }

        //SS & SST matrix
        for (int i = 0;i < 1;i++) {
            for (int j = 0;j < 1;j++) {
                for (int m = 0;m < M;m++) {
                    for (int mp = 0;mp < M - P;mp++) {
                        for (int m2 = 0;m2 < M;m2++) {
                            PP[i][j] += SS[i][m] * NN[m][mp] * NNT[mp][m2] * SST[m2][j];
                            Pmusic = abs(T1/PP[i][j]);
                        }
                    }
                }
            }
        }
       
        if (Pmusic > Max) {
            Max = Pmusic;
            theata = ii;
        }
    }
    cout << theata;
}
