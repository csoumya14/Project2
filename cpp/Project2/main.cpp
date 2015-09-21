#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#include <iomanip>
using std::setw;
#include "include\armadillo"
using namespace arma;

double find_max(mat A, int n)
{
    // Finding maximum matrix element (and its row and column numbers) of the non-diagonal A
        double max = A(0,1);
        int row_number;
        int column_number;
        row_number = 0;
        column_number = 1;
        for (int i=0; i<n; i++) //Compare entrances in the upper part of the matrix. Choose largest one (abs-value).
        {
            for (int j=i+1; j<n; j++)
            {

            if (fabs(A(i,j)) > fabs(max))
            {
                max = A(i,j);
                row_number = i;
                column_number = j;
            }
            }
        }
        for (int i=1; i<n; i++) //Compare entrances in the lower part of the matrix. Choose largest one (abs-value).
        {
            for (int j=0; j<i; j++)
            {

            if (fabs(A(i,j)) > fabs(max))
            {
                max = A(i,j);
                row_number = i;
                column_number = j;
            }
            }
        }
        double MaxInfo[3] = {max, row_number, column_number}; //0: max value, 1: row number, 2: column number
        //MaxInfo(0) = max;
        //MaxInfo(1) = row_number;
        //MaxInfo(2) = column_number;
        return {MaxInfo[0], MaxInfo[1], MaxInfo[2]};
}

double Jacobi (mat A, int n, double epsilon)
{
    vec MaxInfo(3);
    MaxInfo = find_max(A,n);
    double max = MaxInfo(0);
    int row_number = MaxInfo(1);
    int column_number = MaxInfo(2);

    cout << A << max << setw(10) << row_number << setw(10) << column_number << endl;

    // Use simpler test:

    if (fabs(max) > epsilon)
    {
        // computing tau, tan, cos and sine
        double tau;
        double c;
        double t;
        double s;
        tau = (A(column_number, column_number)-A(row_number, row_number))/(2*A(row_number, column_number));
        if (tau >= 0)
        {
            t = (-tau+sqrt(1+pow(tau,2)));
        }
        else
        {
            t = (-tau-sqrt(1+pow(tau,2)));
        }
        c = 1/sqrt(1+pow(t,2));
        s = t*c;

        cout << A << max << setw(10) << row_number << setw(10) << column_number << setw(10) << epsilon << setw(10) << n << endl;

        // Computing the new matrix A

        mat temp(n,n);

        temp = A;

        for (int i = 0; i<n; i++)
        {
            if (i != row_number && i != column_number)  // determining A(i,k) for new matrix
            {
                temp(i,row_number) = A(i,row_number)*c - A(i,column_number)*s;
                temp(i,column_number) = A(i,column_number)*c + A(i,row_number)*s;
            }
        }
        for (int i = 0; i<n; i++)
        {
            if (i != row_number && i != column_number)  // determining A(i,k) for new matrix
            {
                temp(row_number,i) = temp(i,row_number);
                temp(column_number,i) = temp(i,column_number);
            }
        }
        temp(row_number, row_number) = A(row_number, row_number)*pow(c,2) - 2*A(row_number,column_number)*c*s+A(column_number,column_number)*pow(s,2);
        temp(column_number, column_number) = A(column_number, column_number)*pow(c,2) + 2*A(row_number,column_number)*c*s+A(row_number,row_number)*pow(s,2);
        // temp(row_number, column_number) = (A(row_number, row_number)-A(column_number, column_number))*c*s+A(row_number, column_number)*(pow(c,2)-pow(s,2));
        temp(column_number, row_number) = 0.0;  //By choice of theta
        temp(row_number, column_number) = 0.0;  //By choice of theta

        A = temp;
    }
 //   cout << A << endl;
}

int main(int argc, char *argv[])
{
    double rho_min = 0;
    double rho_max = 10;
    double epsilon;
    epsilon = 1.0e-8;
    int n;
    cout << "Please enter value of n:\n>";
    cin >> n;
    cout << "n = " << n << endl;

    double h = (rho_max - rho_min)/(n+1); //step length


    vec V(n);
    for (int i = 0; i < n; i++)
    {
        V(i) = pow(rho_min + i*h,2);
    }

    mat A(n,n);
    A.zeros();

    for (int i=0; i<n; i++)
    {
        A(i,i) = 2/pow(h,2) + V(i);
    }

    double off_diagonal;
    off_diagonal = -1/pow(h,2);

    for (int i=1; i<n; i++)
    {
        A(i,i-1) = off_diagonal;
    }
    for (int i=0; i<n-1; i++)
    {
        A(i,i+1) = off_diagonal;
    }

    int c=2;
    int d=4;

    // int Comp = Jacobi (c,d);
    Jacobi(A,n,epsilon);

    return 0;
}

