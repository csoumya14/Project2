#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#include <iomanip>
using std::setw;
#include "include\armadillo"
using namespace arma;

void find_max(mat &A, int &n, int &row_number, int &column_number)
// set row_number = 0 and column_number = 1, when running the code. These are the initial guesses for max(A(i,j))
{
    // Finding maximum matrix element (and its row and column numbers) of the non-diagonal A
        double max = A(0,1);

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
//we are only considering the upper part of the matrix when seaching for the largest element
//since the matrix is symmetric
/*
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
*/
        return;
}

void Jacobi (mat &A, int n, double epsilon)
{
    int row_number, column_number;
    row_number = 0;
    column_number = 1;
    find_max(A,n,row_number,column_number);
    double max;
    max = A(row_number, column_number);
    /*
    cout << "First A = " << endl;
    cout << A << fabs(max) << setw(10) << row_number << setw(10) << column_number << endl;
    */
    // Use simpler test:
    int m=0;
    while(pow(fabs(max),2)> epsilon && m<20)
    {
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
        temp(column_number, row_number) = 0.00;  //By choice of theta
        temp(row_number, column_number) = 0.00;  //By choice of theta
/*
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<n; j++)
            {
                if (temp(i,j)<epsilon)
                {
                    temp(i,j) = 0.0;
                }
            }
        }
*/
        A = temp;
        row_number = 0;
        column_number = 1;

        find_max(A,n,row_number,column_number);
        m += 1;
        max = A(row_number,column_number);
/*
        cout << "number of iterations = " << m << endl;
        cout << "new A = " << endl;
        cout << A << max << setw(10) << row_number << setw(10) << column_number << endl;
*/
    }
    }
}

int main()
{
    double rho_min = 0;
    double rho_max = 10;
    double epsilon;
    epsilon = 1.0e-8;
    int n;
    cout << "Please enter value of n:\n>";
    cin >> n;
    cout << "n = " << n << endl;
    cout << "rho_min = " << rho_min << endl;
    cout << "rho_max = " << rho_max << endl;

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
/*
    cout << "A=" << endl;
    cout << A << endl;
*/
    Jacobi(A,n,epsilon);

    vec eigen_values(n);

    for (int i=0; i<n; i++)
    {
        eigen_values(i) = A(i,i);
    }

    cout << "eigen values:" << endl;
    cout << eigen_values << endl;

    return 0;
}

