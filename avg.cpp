#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

ofstream ofs("quadratic_avg.csv");

double dt = 0.01;
int N = 40;
double F = 8;

vector<double> lorenz96(const vector<double> &X)
{
    vector<double> dX(N);
    for (int i = 0; i < N; i++)
    {
        int p1 = (i + 1) % N;
        int m1 = (i - 1 + N) % N;
        int m2 = (i - 2 + N) % N;

        dX[i] = (X[p1] - X[m2]) * X[m1] - X[i] + F;
    }
    return dX;
}

double avg(const vector<double> &X)
{
    double sum=0;
    for (int i=0;i<N;i++){
        
        sum+=X[i];
    }
    return sum/N;
}

int main()
{
    ofs << fixed << setprecision(5);
    vector<double> X(N, F);
    vector<double> k1(N), k2(N), k3(N), k4(N), tmp(N);
    vector<double> dX1(N), dX2(N), dX3(N), dX4(N),q(N);

    X[N / 2] += F * 0.01;

    for (int step = 0; step < 40; step++)
    {
// X
        dX1 = lorenz96(X);
        for (int i = 0; i < N; i++)
            k1[i] = dt * dX1[i];

        for (int i = 0; i < N; i++)
            tmp[i] = X[i] + k1[i] / 2;

        dX2 = lorenz96(tmp);
        for (int i = 0; i < N; i++)
            k2[i] = dt * dX2[i];

        for (int i = 0; i < N; i++)
            tmp[i] = X[i] + k2[i] / 2;

        dX3 = lorenz96(tmp);
        for (int i = 0; i < N; i++)
            k3[i] = dt * dX3[i];

        for (int i = 0; i < N; i++)
            tmp[i] = X[i] + k3[i];

        dX4 = lorenz96(tmp);
        for (int i = 0; i < N; i++)
            k4[i] = dt * dX4[i];

        for (int i = 0; i < N; i++)
            X[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;

        for  (int i=0;i<N;i++)
            q[i]=(X[i])*(X[i]);
            
        ofs << log(avg(q))<<endl;
    }
        

    return 0;
}