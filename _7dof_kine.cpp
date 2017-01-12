#include "_7dof_kine.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

//rpy2r
void rpy2r(const double rpy[3],double r[3][3])
{
    r[0][0] = cos(rpy[1])*cos(rpy[2]);
    r[0][1] = -cos(rpy[1])*sin(rpy[2]);
    r[0][2] = sin(rpy[1]);
    r[1][0] = cos(rpy[0])*sin(rpy[2]) + cos(rpy[2])*sin(rpy[0])*sin(rpy[1]);
    r[1][1] = cos(rpy[0])*cos(rpy[2]) - sin(rpy[0])*sin(rpy[1])*sin(rpy[2]);
    r[1][2] = -cos(rpy[1])*sin(rpy[0]);
    r[2][0] = sin(rpy[0])*sin(rpy[2]) - cos(rpy[0])*cos(rpy[2])*sin(rpy[1]);
    r[2][1] = cos(rpy[2])*sin(rpy[0]) + cos(rpy[0])*sin(rpy[1])*sin(rpy[2]);
    r[2][2] = cos(rpy[0])*cos(rpy[1]);
}

//r2rpy
void r2rpy(const double r[3][3], double rpy[3])
{
    double sr,cr;
    if((fabs(r[2][2]) < 1E-8 && fabs(r[1][2]) < 1E-8))
    {
        // singularity
        rpy[0] = 0;  // roll is zero
        rpy[1] = atan2(r[0][2], r[2][2]);  // pitch
        rpy[2] = atan2(r[1][0], r[1][1]);  // yaw is sum of roll + yaw
    }
    else
    {
        rpy[0] = atan2(-r[1][2], r[2][2]);        // roll
        // compute sin / cos of roll angle
        sr = sin(rpy[0]);
        cr = cos(rpy[0]);
        rpy[1] = atan2(r[0][2], cr * r[2][2] - sr * r[1][2]);  // pitch
        rpy[2] = atan2(-r[0][1], r[0][0]);        // yaw
    }
}

//LU decomposition
void LUMethod(const double* A, const int& order, double** LU)
{
    for (int k = 0; k < order; k++)
    {
        for (int j = k; j < order; j++)
        {
            double sumForU = 0;
            for (int t = 0; t <= k - 1; t++)
            {
                sumForU += LU[k][t] * LU[t][j];
            }
            LU[k][j] = A[k*order + j] - sumForU;
        }
        for (int i = k + 1; i < order; i++)
        {
            double sumForL = 0;
            for (int t = 0; t <= k - 1; t++)
            {
                sumForL += LU[i][t] * LU[t][k];
            }
            LU[i][k] = (A[i*order + k] - sumForL) / LU[k][k];
        }
    }
}

//solve Ax=B by LU decomposition
void SolveLinearEquationByLU(double **LUMatrix, const double *b, const int &order, double *x)
{
    double *y = new double[order];
    for (int i = 0; i < order; i++)
    {
        double sumForY = 0;
        for (int t = 0; t <= i - 1; t++)
            sumForY += LUMatrix[i][t] * y[t];
        y[i] = b[i] - sumForY;
    }

    for (int i = order - 1; i >= 0; i--)
    {
        double sumForX = 0;
        for (int t = i + 1; t < order; t++)
            sumForX += LUMatrix[i][t] * x[t];
        x[i] = (y[i] - sumForX) / LUMatrix[i][i];
    }

    delete[] y;
}

//get center eigenvalue by inverse power method
void InversePowerMethod(const double *matrix, const int &order, const double &center, double &outEigenvalue, double *outEigenvector)
{
    //init some variables
    double *u0 = new double[order];
    double *u1 = new double[order];
    double *y0 = new double[order];
    double beta0 = 0, beta1 = 0;
    size_t size = sizeof(u0);
    memset(u0, 0, order*size);

    //init nonzero u0
    u0[0] = 1;

    //translation matrix
    double *matrix_trans = new double[order*order];
    for (int i = 0; i<order*order; i++)
    {
        matrix_trans[i] = matrix[i];
    }
    for (int i = 0; i<order; i++)
        matrix_trans[i*order + i] -= center;

    //matrix_trans LU decomposition
    double **LUMatrix = new double*[order];
    for (int i = 0; i<order; i++)
        LUMatrix[i] = new double[order];
    LUMethod(matrix_trans, order, LUMatrix);

    //****************************************************
    //main process of inverse power method
    do
    {
        beta0 = beta1;		//last beta
        memcpy(u0, u1, order*sizeof(double));


        double lengthOfU = 0.0;
        for (int i = 0; i < order; i++)
            lengthOfU += u0[i] * u0[i];
        lengthOfU = sqrt(lengthOfU);

        for (int i = 0; i<order; i++)
            y0[i] = u0[i] / lengthOfU;

        //A*u1 = y0,get u1
        SolveLinearEquationByLU(LUMatrix, y0, order, u1);

        beta1 = 0;
        for (int i = 0; i < order; i++)
            beta1 += y0[i] * u1[i];
    } while (fabs(beta1 - beta0) / fabs(beta1) > EPS);

    outEigenvalue = (double)1.0 / beta1 + center;
    memcpy(outEigenvector, y0, order*sizeof(double));
    //****************************************************

    delete[](u0, u1, y0);
    for (int i = 0; i<order; i++)
    {
        delete[] LUMatrix[i];
        LUMatrix[i] = 0;
    }
    delete[] LUMatrix;
}

//ψ region of joint_2 or joint_6
void region_cos_type(double a, double b, double c,
                     double joint_u, double joint_l, double* region)
{
    double joint_ul;
    if (joint_u > 0 && joint_l < 0)
    {
        joint_ul = fabs(joint_u) < fabs(joint_l) ? fabs(joint_u) : fabs(joint_l);
    }
    else
        return; //joint region error


    double cos_theta;
    cos_theta = cos(joint_ul);

    double psi_1, psi_2;
    double delta = pow(a, 2) + pow(b, 2) - pow(c - cos_theta, 2);

    if (delta < 0 || fabs(b - c + cos_theta) < ZER)
    {
        region[0] = 1;
        region[1] = -PI;
        region[2] = PI;
    }
    else
    {
        psi_1 = 2 * atan((a - sqrt(delta)) / (b - c + cos_theta));
        psi_2 = 2 * atan((a + sqrt(delta)) / (b - c + cos_theta));
        double psi_u, psi_l;
        psi_u = psi_1 > psi_2 ? psi_1 : psi_2;
        psi_l = psi_1 < psi_2 ? psi_1 : psi_2;

        double psi_mid = (psi_l + psi_u) / 2;
        double theta_at_mid = acos(a*sin(psi_mid) + b*cos(psi_mid) + c);
        if (theta_at_mid < joint_u)
            region[0] = 1;
        else
            region[0] = 2;

        region[1] = psi_l;
        region[2] = psi_u;
    }
}

//get the union of region_1 and region_2
void region_union(const double* region_1, const double* region_2, double* region)
{
    if (fabs(region_1[0] - 1) < ZER && fabs(region_2[0] - 1)<ZER)
    {
        region[0] = 1;
        region[1] = region_1[1] > region_2[1] ? region_1[1] : region_2[1];
        region[2] = region_1[2] < region_2[2] ? region_1[2] : region_2[2];
    }
    else if (fabs(region_1[0] - 2) < ZER && fabs(region_2[0] - 1) < ZER)
    {
        if ((region_2[2]<region_1[1]) || (region_2[1]>region_1[2]))
        {
            region[0] = 1;
            region[1] = region_2[1];
            region[2] = region_2[2];
        }
        else if ((region_2[1]<region_1[1]) && (region_2[2]<region_1[2]))
        {
            region[0] = 1;
            region[1] = region_2[1];
            region[2] = region_1[1];
        }
        else if ((region_2[1] > region_1[1]) && (region_2[2] > region_1[2]))
        {
            region[0] = 1;
            region[1] = region_1[2];
            region[2] = region_2[2];
        }
        else if ((region_2[1] < region_1[1]) && (region_2[2] > region_1[2]))
        {
            region[0] = 2;
            region[1] = region_2[1];
            region[2] = region_1[1];
            region[3] = region_2[2];
            region[4] = region_1[2];
        }
        else region[0] = 3;
    }
    else if (fabs(region_1[0] - 1) < ZER && fabs(region_2[0] - 2) < ZER)
    {
        if ((region_1[2]<region_2[1]) || (region_1[1]>region_2[2]))
        {
            region[0] = 1;
            region[1] = region_1[1];
            region[2] = region_1[2];
        }
        else if ((region_1[1]<region_2[1]) && (region_1[2]<region_2[2]))
        {
            region[0] = 1;
            region[1] = region_1[1];
            region[2] = region_2[1];
        }
        else if ((region_1[1] > region_2[1]) && (region_1[2] > region_2[2]))
        {
            region[0] = 1;
            region[1] = region_2[2];
            region[2] = region_1[2];
        }
        else if ((region_1[1] < region_2[1]) && (region_1[2] > region_2[2]))
        {
            region[0] = 2;
            region[1] = region_1[1];
            region[2] = region_2[1];
            region[3] = region_1[2];
            region[4] = region_2[2];
        }
        else region[0] = 3;
    }
}

//rotate of shoulder
void Rot_shoulder(const double shoulder[3], double R[3][3])
{
    R[0][0] = cos(shoulder[0])*cos(shoulder[1])*cos(shoulder[2]) - sin(shoulder[0])*sin(shoulder[2]);
    R[0][1] = -cos(shoulder[0])*sin(shoulder[1]);
    R[0][2] = -cos(shoulder[2])*sin(shoulder[0]) - cos(shoulder[0])*cos(shoulder[1])*sin(shoulder[2]);

    R[1][0] = cos(shoulder[0])*sin(shoulder[2]) + cos(shoulder[1])*cos(shoulder[2])*sin(shoulder[0]);
    R[1][1] = -sin(shoulder[0])*sin(shoulder[1]);
    R[1][2] = cos(shoulder[0])*cos(shoulder[2]) - cos(shoulder[1])*sin(shoulder[0])*sin(shoulder[2]);

    R[2][0] = -cos(shoulder[2])*sin(shoulder[1]);
    R[2][1] = -cos(shoulder[1]);
    R[2][2] = sin(shoulder[1])*sin(shoulder[2]);
}

//rotate of elbow
void Rot_elbow(double elbow, double R[3][3])
{
    R[0][0] = cos(elbow);	R[0][1] = 0;	R[0][2] = sin(elbow);
    R[1][0] = sin(elbow);	R[1][1] = 0;	R[1][2] = -cos(elbow);
    R[2][0] = 0;			R[2][1] = 1;	R[2][2] = 0;
}

//rotate of wrist
void Rot_wrist(const double wrist[3], double R[3][3])
{
    R[0][0] = cos(wrist[0])*cos(wrist[1])*cos(wrist[2]) - sin(wrist[0])*sin(wrist[2]);
    R[0][1] = -cos(wrist[2])*sin(wrist[0]) - cos(wrist[0])*cos(wrist[1])*sin(wrist[2]);
    R[0][2] = cos(wrist[0])*sin(wrist[1]);

    R[1][0] = cos(wrist[0])*sin(wrist[2]) + cos(wrist[1])*cos(wrist[2])*sin(wrist[0]);
    R[1][1] = cos(wrist[0])*cos(wrist[2]) - cos(wrist[1])*sin(wrist[0])*sin(wrist[2]);
    R[1][2] = sin(wrist[0])*sin(wrist[1]);

    R[2][0] = -cos(wrist[2])*sin(wrist[1]);
    R[2][1] = sin(wrist[1])*sin(wrist[2]);
    R[2][2] = cos(wrist[1]);
}

//forward kinematics
void Fkine(const double joint[7], double T[4][4])
{
    T[0][0] = cos(joint[6])*(sin(joint[5])*(sin(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) - cos(joint[0])*cos(joint[3])*sin(joint[1])) - cos(joint[5])*(cos(joint[4])*(cos(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) + cos(joint[0])*sin(joint[1])*sin(joint[3])) + sin(joint[4])*(cos(joint[2])*sin(joint[0]) + cos(joint[0])*cos(joint[1])*sin(joint[2])))) + sin(joint[6])*(sin(joint[4])*(cos(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) + cos(joint[0])*sin(joint[1])*sin(joint[3])) - cos(joint[4])*(cos(joint[2])*sin(joint[0]) + cos(joint[0])*cos(joint[1])*sin(joint[2])));
    T[0][1] = cos(joint[6])*(sin(joint[4])*(cos(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) + cos(joint[0])*sin(joint[1])*sin(joint[3])) - cos(joint[4])*(cos(joint[2])*sin(joint[0]) + cos(joint[0])*cos(joint[1])*sin(joint[2]))) - sin(joint[6])*(sin(joint[5])*(sin(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) - cos(joint[0])*cos(joint[3])*sin(joint[1])) - cos(joint[5])*(cos(joint[4])*(cos(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) + cos(joint[0])*sin(joint[1])*sin(joint[3])) + sin(joint[4])*(cos(joint[2])*sin(joint[0]) + cos(joint[0])*cos(joint[1])*sin(joint[2]))));
    T[0][2] = -cos(joint[5])*(sin(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) - cos(joint[0])*cos(joint[3])*sin(joint[1])) - sin(joint[5])*(cos(joint[4])*(cos(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) + cos(joint[0])*sin(joint[1])*sin(joint[3])) + sin(joint[4])*(cos(joint[2])*sin(joint[0]) + cos(joint[0])*cos(joint[1])*sin(joint[2])));
    T[0][3] = Dse*cos(joint[0])*sin(joint[1]) - Dwt*(cos(joint[5])*(sin(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) - cos(joint[0])*cos(joint[3])*sin(joint[1])) + sin(joint[5])*(cos(joint[4])*(cos(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) + cos(joint[0])*sin(joint[1])*sin(joint[3])) + sin(joint[4])*(cos(joint[2])*sin(joint[0]) + cos(joint[0])*cos(joint[1])*sin(joint[2])))) - Dew*(sin(joint[3])*(sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[2])) - cos(joint[0])*cos(joint[3])*sin(joint[1]));

    T[1][0] = -cos(joint[6])*(sin(joint[5])*(sin(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) + cos(joint[3])*sin(joint[0])*sin(joint[1])) - cos(joint[5])*(cos(joint[4])*(cos(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) - sin(joint[0])*sin(joint[1])*sin(joint[3])) + sin(joint[4])*(cos(joint[0])*cos(joint[2]) - cos(joint[1])*sin(joint[0])*sin(joint[2])))) - sin(joint[6])*(sin(joint[4])*(cos(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) - sin(joint[0])*sin(joint[1])*sin(joint[3])) - cos(joint[4])*(cos(joint[0])*cos(joint[2]) - cos(joint[1])*sin(joint[0])*sin(joint[2])));
    T[1][1] = sin(joint[6])*(sin(joint[5])*(sin(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) + cos(joint[3])*sin(joint[0])*sin(joint[1])) - cos(joint[5])*(cos(joint[4])*(cos(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) - sin(joint[0])*sin(joint[1])*sin(joint[3])) + sin(joint[4])*(cos(joint[0])*cos(joint[2]) - cos(joint[1])*sin(joint[0])*sin(joint[2])))) - cos(joint[6])*(sin(joint[4])*(cos(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) - sin(joint[0])*sin(joint[1])*sin(joint[3])) - cos(joint[4])*(cos(joint[0])*cos(joint[2]) - cos(joint[1])*sin(joint[0])*sin(joint[2])));
    T[1][2] = cos(joint[5])*(sin(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) + cos(joint[3])*sin(joint[0])*sin(joint[1])) + sin(joint[5])*(cos(joint[4])*(cos(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) - sin(joint[0])*sin(joint[1])*sin(joint[3])) + sin(joint[4])*(cos(joint[0])*cos(joint[2]) - cos(joint[1])*sin(joint[0])*sin(joint[2])));
    T[1][3] = Dew*(sin(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) + cos(joint[3])*sin(joint[0])*sin(joint[1])) + Dwt*(cos(joint[5])*(sin(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) + cos(joint[3])*sin(joint[0])*sin(joint[1])) + sin(joint[5])*(cos(joint[4])*(cos(joint[3])*(cos(joint[0])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])) - sin(joint[0])*sin(joint[1])*sin(joint[3])) + sin(joint[4])*(cos(joint[0])*cos(joint[2]) - cos(joint[1])*sin(joint[0])*sin(joint[2])))) + Dse*sin(joint[0])*sin(joint[1]);

    T[2][0] = sin(joint[6])*(sin(joint[4])*(cos(joint[1])*sin(joint[3]) + cos(joint[2])*cos(joint[3])*sin(joint[1])) + cos(joint[4])*sin(joint[1])*sin(joint[2])) - cos(joint[6])*(cos(joint[5])*(cos(joint[4])*(cos(joint[1])*sin(joint[3]) + cos(joint[2])*cos(joint[3])*sin(joint[1])) - sin(joint[1])*sin(joint[2])*sin(joint[4])) + sin(joint[5])*(cos(joint[1])*cos(joint[3]) - cos(joint[2])*sin(joint[1])*sin(joint[3])));
    T[2][1] = cos(joint[6])*(sin(joint[4])*(cos(joint[1])*sin(joint[3]) + cos(joint[2])*cos(joint[3])*sin(joint[1])) + cos(joint[4])*sin(joint[1])*sin(joint[2])) + sin(joint[6])*(cos(joint[5])*(cos(joint[4])*(cos(joint[1])*sin(joint[3]) + cos(joint[2])*cos(joint[3])*sin(joint[1])) - sin(joint[1])*sin(joint[2])*sin(joint[4])) + sin(joint[5])*(cos(joint[1])*cos(joint[3]) - cos(joint[2])*sin(joint[1])*sin(joint[3])));
    T[2][2] = cos(joint[5])*(cos(joint[1])*cos(joint[3]) - cos(joint[2])*sin(joint[1])*sin(joint[3])) - sin(joint[5])*(cos(joint[4])*(cos(joint[1])*sin(joint[3]) + cos(joint[2])*cos(joint[3])*sin(joint[1])) - sin(joint[1])*sin(joint[2])*sin(joint[4]));
    T[2][3] = Dbs + Dew*(cos(joint[1])*cos(joint[3]) - cos(joint[2])*sin(joint[1])*sin(joint[3])) + Dse*cos(joint[1]) - Dwt*(sin(joint[5])*(cos(joint[4])*(cos(joint[1])*sin(joint[3]) + cos(joint[2])*cos(joint[3])*sin(joint[1])) - sin(joint[1])*sin(joint[2])*sin(joint[4])) - cos(joint[5])*(cos(joint[1])*cos(joint[3]) - cos(joint[2])*sin(joint[1])*sin(joint[3])));

    T[3][0] = 0;
    T[3][1] = 0;
    T[3][2] = 0;
    T[3][3] = 1;
}

void Fkine_rpy(const double joint[7], double end[6])
{

    double T[4][4],r[3][3],rpy[3];
    Fkine(joint,T);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            r[i][j]=T[i][j];
    r2rpy(r,rpy);
    end[0]=T[0][3];
    end[1]=T[1][3];
    end[2]=T[2][3];
    end[3]=rpy[0];
    end[4]=rpy[1];
    end[5]=rpy[2];
}

//inverse kinematics
void Ikine(const double T[4][4], double joint[7],double* psi_last)
{}

void Ikine_rpy(const double end[6], double joint[7], double* psi_last)
{
    double rpy[3],r[3][3],T[4][4];
    rpy[0]=end[3];
    rpy[1]=end[4];
    rpy[2]=end[5];
    rpy2r(rpy,r);

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            T[i][j]=r[i][j];

    T[0][3]=end[0];
    T[1][3]=end[1];
    T[2][3]=end[2];

    T[3][0]=0;
    T[3][1]=0;
    T[3][2]=0;
    T[3][3]=1;

    Ikine(T,joint,psi_last);
}

void Jacob(const double joint[7], double J[6][7])
{
    J[0][0] = -Dwt*sin(joint[5])*sin(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[3])*cos(joint[4]) + Dwt*sin(joint[5])*sin(joint[3])*sin(joint[0])*sin(joint[1])*cos(joint[4]) + Dwt*sin(joint[5])*sin(joint[0])*sin(joint[2])*cos(joint[1])*sin(joint[4]) - Dwt*sin(joint[5])*sin(joint[2])*cos(joint[0])*cos(joint[3])*cos(joint[4]) - Dwt*sin(joint[3])*sin(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[5]) - Dew*sin(joint[3])*sin(joint[0])*cos(joint[1])*cos(joint[2]) - Dwt*sin(joint[5])*cos(joint[0])*cos(joint[2])*sin(joint[4]) - Dwt*sin(joint[3])*sin(joint[2])*cos(joint[0])*cos(joint[5]) - Dwt*sin(joint[0])*cos(joint[3])*sin(joint[1])*cos(joint[5]) - Dew*sin(joint[3])*sin(joint[2])*cos(joint[0]) - Dew*sin(joint[0])*cos(joint[3])*sin(joint[1]) - Dse*sin(joint[0])*sin(joint[1]);
    J[0][1] = -cos(joint[0])*(Dwt*sin(joint[5])*cos(joint[2])*cos(joint[3])*sin(joint[1])*cos(joint[4]) + Dwt*sin(joint[5])*sin(joint[3])*cos(joint[1])*cos(joint[4]) - Dwt*sin(joint[5])*sin(joint[2])*sin(joint[1])*sin(joint[4]) + Dwt*sin(joint[3])*cos(joint[2])*sin(joint[1])*cos(joint[5]) + Dew*sin(joint[3])*cos(joint[2])*sin(joint[1]) - Dwt*cos(joint[1])*cos(joint[3])*cos(joint[5]) - Dew*cos(joint[1])*cos(joint[3]) - Dse*cos(joint[1]));
    J[0][2] = -Dwt*sin(joint[5])*sin(joint[2])*cos(joint[0])*cos(joint[1])*cos(joint[3])*cos(joint[4]) - Dwt*sin(joint[5])*sin(joint[0])*cos(joint[2])*cos(joint[3])*cos(joint[4]) - Dwt*sin(joint[5])*cos(joint[0])*cos(joint[1])*cos(joint[2])*sin(joint[4]) - Dwt*sin(joint[3])*sin(joint[2])*cos(joint[0])*cos(joint[1])*cos(joint[5]) - Dew*sin(joint[3])*sin(joint[2])*cos(joint[0])*cos(joint[1]) + Dwt*sin(joint[5])*sin(joint[0])*sin(joint[2])*sin(joint[4]) - Dwt*sin(joint[3])*sin(joint[0])*cos(joint[2])*cos(joint[5]) - Dew*sin(joint[3])*sin(joint[0])*cos(joint[2]);
    J[0][3] = -Dwt*sin(joint[5])*sin(joint[3])*cos(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[4]) + Dwt*sin(joint[5])*sin(joint[3])*sin(joint[0])*sin(joint[2])*cos(joint[4]) - Dwt*sin(joint[5])*cos(joint[0])*cos(joint[3])*sin(joint[1])*cos(joint[4]) + Dwt*cos(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[3])*cos(joint[5]) + Dew*cos(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[3]) - Dwt*sin(joint[3])*cos(joint[0])*sin(joint[1])*cos(joint[5]) - Dwt*sin(joint[0])*sin(joint[2])*cos(joint[3])*cos(joint[5]) - Dew*sin(joint[3])*cos(joint[0])*sin(joint[1]) - Dew*sin(joint[0])*sin(joint[2])*cos(joint[3]);
    J[0][4] = Dwt*sin(joint[5])*(-cos(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[3])*sin(joint[4]) + cos(joint[0])*sin(joint[1])*sin(joint[3])*sin(joint[4]) + cos(joint[3])*sin(joint[0])*sin(joint[2])*sin(joint[4]) - cos(joint[0])*cos(joint[1])*cos(joint[4])*sin(joint[2]) - cos(joint[2])*cos(joint[4])*sin(joint[0]));
    J[0][5] = Dwt*(cos(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[3])*cos(joint[4])*cos(joint[5]) - cos(joint[0])*cos(joint[1])*cos(joint[2])*sin(joint[3])*sin(joint[5]) - cos(joint[0])*cos(joint[4])*cos(joint[5])*sin(joint[1])*sin(joint[3]) - cos(joint[3])*cos(joint[4])*cos(joint[5])*sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[1])*cos(joint[5])*sin(joint[2])*sin(joint[4]) + sin(joint[0])*sin(joint[2])*sin(joint[3])*sin(joint[5]) - cos(joint[0])*cos(joint[3])*sin(joint[1])*sin(joint[5]) - cos(joint[2])*cos(joint[5])*sin(joint[0])*sin(joint[4]));
    J[0][6] = 0;

    J[1][0] = Dwt*cos(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[3])*sin(joint[5])*cos(joint[4]) + Dwt*cos(joint[5])*sin(joint[3])*cos(joint[0])*cos(joint[1])*cos(joint[2]) - Dwt*sin(joint[3])*cos(joint[0])*sin(joint[1])*sin(joint[5])*cos(joint[4]) - Dwt*cos(joint[0])*sin(joint[2])*cos(joint[1])*sin(joint[5])*sin(joint[4]) - Dwt*sin(joint[2])*sin(joint[0])*cos(joint[3])*sin(joint[5])*cos(joint[4]) + Dew*sin(joint[3])*cos(joint[0])*cos(joint[1])*cos(joint[2]) - Dwt*cos(joint[5])*sin(joint[3])*sin(joint[2])*sin(joint[0]) + Dwt*cos(joint[5])*cos(joint[0])*cos(joint[3])*sin(joint[1]) - Dwt*cos(joint[2])*sin(joint[0])*sin(joint[5])*sin(joint[4]) - Dew*sin(joint[3])*sin(joint[2])*sin(joint[0]) + Dew*cos(joint[0])*cos(joint[3])*sin(joint[1]) + Dse*cos(joint[0])*sin(joint[1]);
    J[1][1] = -sin(joint[0])*(Dwt*cos(joint[2])*cos(joint[3])*sin(joint[1])*sin(joint[5])*cos(joint[4]) + Dwt*cos(joint[5])*sin(joint[3])*cos(joint[2])*sin(joint[1]) + Dwt*sin(joint[3])*cos(joint[1])*sin(joint[5])*cos(joint[4]) - Dwt*sin(joint[2])*sin(joint[1])*sin(joint[5])*sin(joint[4]) + Dew*sin(joint[3])*cos(joint[2])*sin(joint[1]) - Dwt*cos(joint[5])*cos(joint[1])*cos(joint[3]) - Dew*cos(joint[1])*cos(joint[3]) - Dse*cos(joint[1]));
    J[1][2] = -Dwt*sin(joint[2])*cos(joint[1])*sin(joint[0])*cos(joint[3])*sin(joint[5])*cos(joint[4]) - Dwt*cos(joint[5])*sin(joint[3])*sin(joint[2])*cos(joint[1])*sin(joint[0]) + Dwt*cos(joint[0])*cos(joint[2])*cos(joint[3])*sin(joint[5])*cos(joint[4]) - Dwt*cos(joint[1])*cos(joint[2])*sin(joint[0])*sin(joint[5])*sin(joint[4]) - Dew*sin(joint[3])*sin(joint[2])*cos(joint[1])*sin(joint[0]) + Dwt*cos(joint[5])*sin(joint[3])*cos(joint[0])*cos(joint[2]) - Dwt*cos(joint[0])*sin(joint[2])*sin(joint[5])*sin(joint[4]) + Dew*sin(joint[3])*cos(joint[0])*cos(joint[2]);
    J[1][3] = -Dwt*cos(joint[4])*sin(joint[5])*sin(joint[3])*cos(joint[1])*cos(joint[2])*sin(joint[0]) - Dwt*cos(joint[4])*sin(joint[5])*sin(joint[3])*cos(joint[0])*sin(joint[2]) - Dwt*cos(joint[4])*sin(joint[5])*sin(joint[0])*cos(joint[3])*sin(joint[1]) + Dwt*cos(joint[5])*cos(joint[1])*cos(joint[2])*sin(joint[0])*cos(joint[3]) + Dew*cos(joint[1])*cos(joint[2])*sin(joint[0])*cos(joint[3]) - Dwt*cos(joint[5])*sin(joint[3])*sin(joint[0])*sin(joint[1]) + Dwt*cos(joint[5])*cos(joint[0])*sin(joint[2])*cos(joint[3]) - Dew*sin(joint[3])*sin(joint[0])*sin(joint[1]) + Dew*cos(joint[0])*sin(joint[2])*cos(joint[3]);
    J[1][4] = Dwt*sin(joint[5])*(-cos(joint[1])*cos(joint[2])*cos(joint[3])*sin(joint[0])*sin(joint[4]) - cos(joint[1])*cos(joint[4])*sin(joint[0])*sin(joint[2]) + sin(joint[0])*sin(joint[1])*sin(joint[3])*sin(joint[4]) - cos(joint[0])*cos(joint[3])*sin(joint[2])*sin(joint[4]) + cos(joint[0])*cos(joint[2])*cos(joint[4]));
    J[1][5] = -Dwt*(cos(joint[3])*sin(joint[0])*sin(joint[1])*sin(joint[5]) - cos(joint[0])*cos(joint[2])*cos(joint[5])*sin(joint[4]) + cos(joint[0])*sin(joint[2])*sin(joint[3])*sin(joint[5]) - cos(joint[0])*cos(joint[3])*cos(joint[4])*cos(joint[5])*sin(joint[2]) + cos(joint[1])*cos(joint[2])*sin(joint[0])*sin(joint[3])*sin(joint[5]) + cos(joint[1])*cos(joint[5])*sin(joint[0])*sin(joint[2])*sin(joint[4]) + cos(joint[4])*cos(joint[5])*sin(joint[0])*sin(joint[1])*sin(joint[3]) - cos(joint[1])*cos(joint[2])*cos(joint[3])*cos(joint[4])*cos(joint[5])*sin(joint[0]));
    J[1][6] = 0;

    J[2][0] = 0;
    J[2][1] = -Dwt*cos(joint[1])*cos(joint[2])*cos(joint[3])*cos(joint[4])*sin(joint[5]) + Dwt*sin(joint[4])*cos(joint[1])*sin(joint[2])*sin(joint[5]) - Dwt*cos(joint[1])*sin(joint[3])*cos(joint[2])*cos(joint[5]) + Dwt*sin(joint[3])*sin(joint[1])*cos(joint[4])*sin(joint[5]) - Dew*cos(joint[1])*sin(joint[3])*cos(joint[2]) - Dwt*cos(joint[3])*sin(joint[1])*cos(joint[5]) - Dew*cos(joint[3])*sin(joint[1]) - Dse*sin(joint[1]);
    J[2][2] = sin(joint[1])*(Dwt*cos(joint[3])*cos(joint[4])*sin(joint[2])*sin(joint[5]) + Dwt*sin(joint[4])*cos(joint[2])*sin(joint[5]) + Dwt*sin(joint[3])*sin(joint[2])*cos(joint[5]) + Dew*sin(joint[3])*sin(joint[2]));
    J[2][3] = Dwt*sin(joint[3])*cos(joint[2])*sin(joint[1])*cos(joint[4])*sin(joint[5]) - Dwt*cos(joint[1])*cos(joint[3])*cos(joint[4])*sin(joint[5]) - Dwt*cos(joint[2])*cos(joint[3])*sin(joint[1])*cos(joint[5]) - Dew*cos(joint[2])*cos(joint[3])*sin(joint[1]) - Dwt*cos(joint[1])*sin(joint[3])*cos(joint[5]) - Dew*cos(joint[1])*sin(joint[3]);
    J[2][4] = Dwt*sin(joint[5])*(cos(joint[4])*sin(joint[1])*sin(joint[2]) + cos(joint[1])*sin(joint[3])*sin(joint[4]) + cos(joint[2])*cos(joint[3])*sin(joint[1])*sin(joint[4]));
    J[2][5] = Dwt*(-cos(joint[2])*cos(joint[3])*cos(joint[4])*cos(joint[5])*sin(joint[1]) + cos(joint[2])*sin(joint[1])*sin(joint[3])*sin(joint[5]) - cos(joint[1])*cos(joint[4])*cos(joint[5])*sin(joint[3]) + cos(joint[5])*sin(joint[1])*sin(joint[2])*sin(joint[4]) - cos(joint[1])*cos(joint[3])*sin(joint[5]));
    J[2][6] = 0;

    J[3][0] = 0;
    J[3][1] = -sin(joint[0]);
    J[3][2] = cos(joint[0])*sin(joint[1]);
    J[3][3] = -cos(joint[0])*cos(joint[1])*sin(joint[2]) - cos(joint[2])*sin(joint[0]);
    J[3][4] = sin(joint[3])*cos(joint[0])*cos(joint[1])*cos(joint[2]) - sin(joint[3])*sin(joint[0])*sin(joint[2]) + cos(joint[0])*cos(joint[3])*sin(joint[1]);
    J[3][5] = cos(joint[0])*sin(joint[1])*sin(joint[3])*sin(joint[4]) - cos(joint[0])*cos(joint[1])*cos(joint[4])*sin(joint[2]) - cos(joint[2])*cos(joint[4])*sin(joint[0]) + cos(joint[3])*sin(joint[0])*sin(joint[2])*sin(joint[4]) - cos(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[3])*sin(joint[4]);
    J[3][6] = sin(joint[5])*cos(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[3])*cos(joint[4]) - sin(joint[5])*sin(joint[3])*cos(joint[0])*sin(joint[1])*cos(joint[4]) - sin(joint[5])*sin(joint[0])*sin(joint[2])*cos(joint[3])*cos(joint[4]) - sin(joint[5])*sin(joint[2])*cos(joint[0])*cos(joint[1])*sin(joint[4]) + sin(joint[3])*cos(joint[0])*cos(joint[1])*cos(joint[2])*cos(joint[5]) - sin(joint[5])*sin(joint[0])*cos(joint[2])*sin(joint[4]) - sin(joint[3])*sin(joint[0])*sin(joint[2])*cos(joint[5]) + cos(joint[0])*cos(joint[3])*sin(joint[1])*cos(joint[5]);

    J[4][0] = 0;
    J[4][1] = cos(joint[0]);
    J[4][2] = sin(joint[0])*sin(joint[1]);
    J[4][3] = cos(joint[0])*cos(joint[2]) - cos(joint[1])*sin(joint[0])*sin(joint[2]);
    J[4][4] = sin(joint[3])*cos(joint[1])*cos(joint[2])*sin(joint[0]) + sin(joint[3])*cos(joint[0])*sin(joint[2]) + cos(joint[3])*sin(joint[0])*sin(joint[1]);
    J[4][5] = cos(joint[0])*cos(joint[2])*cos(joint[4]) - cos(joint[1])*cos(joint[4])*sin(joint[0])*sin(joint[2]) - cos(joint[0])*cos(joint[3])*sin(joint[2])*sin(joint[4]) + sin(joint[0])*sin(joint[1])*sin(joint[3])*sin(joint[4]) - cos(joint[1])*cos(joint[2])*cos(joint[3])*sin(joint[0])*sin(joint[4]);
    J[4][6] = cos(joint[2])*cos(joint[4])*cos(joint[1])*sin(joint[0])*cos(joint[3])*sin(joint[5]) + cos(joint[0])*cos(joint[4])*sin(joint[2])*cos(joint[3])*sin(joint[5]) + cos(joint[2])*cos(joint[1])*sin(joint[0])*sin(joint[3])*cos(joint[5]) - cos(joint[4])*sin(joint[0])*sin(joint[1])*sin(joint[3])*sin(joint[5]) - cos(joint[1])*sin(joint[0])*sin(joint[2])*sin(joint[4])*sin(joint[5]) + cos(joint[0])*cos(joint[2])*sin(joint[4])*sin(joint[5]) + cos(joint[0])*sin(joint[2])*sin(joint[3])*cos(joint[5]) + sin(joint[0])*cos(joint[3])*sin(joint[1])*cos(joint[5]);

    J[5][0] = 1;
    J[5][1] = 0;
    J[5][2] = cos(joint[1]);
    J[5][3] = sin(joint[1])*sin(joint[2]);
    J[5][4] = cos(joint[1])*cos(joint[3]) - cos(joint[2])*sin(joint[1])*sin(joint[3]);
    J[5][5] = cos(joint[4])*sin(joint[1])*sin(joint[2]) + cos(joint[1])*sin(joint[3])*sin(joint[4]) + cos(joint[2])*cos(joint[3])*sin(joint[1])*sin(joint[4]);
    J[5][6] = cos(joint[4])*sin(joint[1])*cos(joint[2])*cos(joint[3])*sin(joint[5]) - cos(joint[4])*cos(joint[1])*sin(joint[3])*sin(joint[5]) + sin(joint[1])*sin(joint[2])*sin(joint[4])*sin(joint[5]) - sin(joint[1])*sin(joint[3])*cos(joint[2])*cos(joint[5]) + cos(joint[1])*cos(joint[3])*cos(joint[5]);
}

//manipulability
double Manip(const double J[6][7])
{
    double JJT[6][6];
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
        {
            JJT[i][j] = J[i][0] * J[j][0] + J[i][1] * J[j][1] + J[i][2] * J[j][2] + J[i][3] * J[j][3] + J[i][4] * J[j][4] + J[i][5] * J[j][5] + J[i][6] * J[j][6];
        }

    double matrix[6 * 6];
    double center = 0.0;
    double eigenvalue = 0.0;
    double eigenvector[6];

    memcpy(matrix, JJT, 6 * 6 * sizeof(double));
    InversePowerMethod(matrix, 6, center, eigenvalue, eigenvector);
    return eigenvalue;

}

