//	perc3d.cpp													
//	Source code for the simulation of electrical percolation of	CNT/polymer composite
//	developed by Sungmin Jung, Ph.D.
//  University of Cambridge
//  Last update on 18th April 2019
//  Contact: sm0104.jung@gmail.com, sj569@cam.ac.uk, s.jung@eng.cam.ac.uk


#include "cnt.h"
#include "perc3d.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <omp.h>

CPERC3D::CPERC3D()
{
	m_cCNT = 0;
	int m_nNoofCNTs = 0;
	int m_nMNoofCNTsPerc = 0;

	int m_nNoofConductionChannelCNT = 2; // ea.
	double m_dDeltaE = 1.0; // eV
	double m_dDistanceCutoff = 0.0014; // um

	double m_dTolerance = 1.0e-8;
}

CPERC3D::~CPERC3D()
{
	if (m_cCNT != 0) delete[] m_cCNT;
}

void CPERC3D::CopyCNTs(CCNT *cnt, int nNoofCNTs)
{
	m_nNoofCNTs = nNoofCNTs;

	if (m_cCNT != 0) delete[] m_cCNT;

	m_cCNT = new CCNT[m_nNoofCNTs];

	for (int i = 0; i < m_nNoofCNTs; i++)
	{
		m_cCNT[i] = cnt[i];
	}
}

void CPERC3D::SetDomain(eGEOM eGeom, eBCND eBC, double dXmax, double dXmin, double dYmax, double dYmin, double dZmax, double dZmin)
{
	m_eGeom = eGeom;
	m_eBC = eBC;

	m_dXmax = dXmax;
	m_dXmin = dXmin;
	m_dYmax = dYmax;
	m_dYmin = dYmin;
	m_dZmax = dZmax;
	m_dZmin = dZmin;
}

void CPERC3D::SetMaterialParam(double dDensityPolymer, double dAppliedVoltage, int nNoofConductionChannelCNT, double dDeltaE, double dDistanceCutoff)
{
	m_dDensityPolymer = dDensityPolymer;
	m_dAppliedVoltage = dAppliedVoltage;
	m_nNoofConductionChannelCNT = nNoofConductionChannelCNT;
	m_dDeltaE = dDeltaE;
	m_dDistanceCutoff = dDistanceCutoff;
}

void CPERC3D::WeightAndVolumeFraction(double *dWP, double *dVP)
{
	double dWeightTotalCNT = 0;
	double dVolumeTotalCNT = 0;
	for (int i = 0; i < m_nNoofCNTs; i++)
	{

		dWeightTotalCNT += m_cCNT[i].WeightCNT();
		dVolumeTotalCNT += m_cCNT[i].VolumeCNT();
	}

	// Caculation Domain Setting
	double dVolumeTotalPolymer;
	if (m_eGeom == BOX) dVolumeTotalPolymer = (m_dXmax - m_dXmin)*(m_dYmax - m_dYmin)*(m_dZmax - m_dZmin) - dVolumeTotalCNT; // um^3
	if (m_eGeom == CYL) dVolumeTotalPolymer = M_PI * ((m_dXmax - m_dXmin) / 2.)*((m_dYmax - m_dYmin) / 2.)*(m_dZmax - m_dZmin) - dVolumeTotalCNT; // um^3

	double dWeightTotalPolymer = m_dDensityPolymer * dVolumeTotalPolymer; // gram

	*dWP = dWeightTotalCNT / (dWeightTotalPolymer + dWeightTotalCNT);
	*dVP = dVolumeTotalCNT / (dVolumeTotalPolymer + dVolumeTotalCNT);

	printf("Calculation of Fraction\n");
	printf("- Weight Percent: %f [%%]\n", *dWP * 100.);
	printf("- Volume Percent: %f [%%]\n", *dVP * 100.);
}

void CPERC3D::Stretching(double dPoissonPolymer, double dPoissonCNT, double dStrainPolymerLong, double *dXmax, double *dXmin, double *dYmax, double *dYmin, double *dZmax, double *dZmin)
{
	double dConstLong = 1. + dStrainPolymerLong;
	double dConstTran = pow(dConstLong, -dPoissonPolymer);

	m_dXmax = dConstTran*m_dXmax;
	m_dYmax = dConstTran*m_dYmax;
	m_dZmax = dConstLong*m_dZmax;

	m_dXmin = dConstTran*m_dXmin;
	m_dYmin = dConstTran*m_dYmin;
	m_dZmin = dConstLong*m_dZmin;

	*dXmax = dConstTran * (*dXmax);
	*dYmax = dConstTran * (*dYmax);
	*dZmax = dConstLong * (*dZmax);

	*dXmin = dConstTran * (*dXmin);
	*dYmin = dConstTran * (*dYmin);
	*dZmin = dConstLong * (*dZmin);

	int i;
#pragma omp parallel for private(i) schedule(dynamic)
	for (i = 0; i < m_nNoofCNTs; i++)
	{
		m_cCNT[i].Stretching(dPoissonPolymer, dPoissonCNT, dStrainPolymerLong);
		m_cCNT[i].SetCNTRange();
	}

}


void CPERC3D::Stretching2(double dPoissonPolymer, double dStrainPolymerLong, double *dXmax, double *dXmin, double *dYmax, double *dYmin, double *dZmax, double *dZmin)
{
	double dConstLong = 1. + dStrainPolymerLong;
	double dConstTran = pow(dConstLong, -dPoissonPolymer);

	m_dXmax = dConstTran * m_dXmax;
	m_dYmax = dConstTran * m_dYmax;
	m_dZmax = dConstLong * m_dZmax;

	m_dXmin = dConstTran * m_dXmin;
	m_dYmin = dConstTran * m_dYmin;
	m_dZmin = dConstLong * m_dZmin;

	*dXmax = dConstTran * (*dXmax);
	*dYmax = dConstTran * (*dYmax);
	*dZmax = dConstLong * (*dZmax);

	*dXmin = dConstTran * (*dXmin);
	*dYmin = dConstTran * (*dYmin);
	*dZmin = dConstLong * (*dZmin);
}

void CPERC3D::Stretching3(double dPoissonPolymer, double dStrainPolymerLong)
{
//	double dConstTran = 1. + dStrainPolymerTran;
	double dConstLong = 1. + dStrainPolymerLong;
	double dConstTran = pow(dConstLong, -dPoissonPolymer);

	m_dXmax = dConstTran * m_dXmax;
	m_dYmax = dConstTran * m_dYmax;
	m_dZmax = dConstLong * m_dZmax;

	m_dXmin = dConstTran * m_dXmin;
	m_dYmin = dConstTran * m_dYmin;
	m_dZmin = dConstLong * m_dZmin;

}



void CPERC3D::GlobalNodeNumbering(int nNoofCNTsPerc, int *plist, int *nNoofTotalNodePerc)
{
	// Global Numbering of Nodal points on CNTs in the percolation network
	*nNoofTotalNodePerc = 0;
	for (int i = 0; i < nNoofCNTsPerc; i++)
	{
		int I = plist[i];
		for (int p = 0; p < m_cCNT[I].m_nNoofNodes; p++)
		{
			if (m_cCNT[I].m_eNodalProperty[p] == MED)
			{
				m_cCNT[I].m_nGIndexList[p] = *nNoofTotalNodePerc;
				(*nNoofTotalNodePerc)++;
			}
			if (m_cCNT[I].m_eNodalProperty[p] == TPE)
			{
				m_cCNT[I].m_nGIndexList[p] = -1;
			}
			if (m_cCNT[I].m_eNodalProperty[p] == BTE)
			{
				m_cCNT[I].m_nGIndexList[p] = -2;
			}
		}
	}
}

void CPERC3D::InsertNodesOnPercNetwork(SLinkedCond *linkperc, int nNoofLinkPerc, int *nNoofTotalNodePerc)
{

	for (int l = 0; l < nNoofLinkPerc; l++)
	{
		int i = linkperc[l].m_nCNTNo1;
		int j = linkperc[l].m_nCNTNo2;

		double U = linkperc[l].m_dU;
		double V = linkperc[l].m_dV;

		eMEDIUM eMed1, eMed2;

		linkperc[l].m_nGIndex1 = m_cCNT[i].InsertNode(U, nNoofTotalNodePerc, &eMed1);
		linkperc[l].m_nGIndex2 = m_cCNT[j].InsertNode(V, nNoofTotalNodePerc, &eMed2);

		linkperc[l].m_eMedium1 = eMed1;
		linkperc[l].m_eMedium2 = eMed2;
	}
}

void CPERC3D::GlobalNodeNumberingReOrder(int nNoofCNTsPerc, int *plist, int nNoofTotalNodePerc, int nNoofLinkPerc, SLinkedCond *linkperc)
{
	int *index_change_table = new int[nNoofTotalNodePerc];
	// Global Numbering of Nodal points on CNTs in the percolation network
	int nGIndex = 0;
	for (int i = 0; i < nNoofCNTsPerc; i++)
	{
		int I = plist[i];
		for (int p = 0; p < m_cCNT[I].m_nNoofNodes; p++)
		{
			if (m_cCNT[I].m_eNodalProperty[p] == MED)
			{
				index_change_table[m_cCNT[I].m_nGIndexList[p]] = nGIndex;
				m_cCNT[I].m_nGIndexList[p] = nGIndex;
				(nGIndex)++;
			}
			if (m_cCNT[I].m_eNodalProperty[p] == TPE)
			{
				m_cCNT[I].m_nGIndexList[p] = -1;
			}
			if (m_cCNT[I].m_eNodalProperty[p] == BTE)
			{
				m_cCNT[I].m_nGIndexList[p] = -2;
			}
		}
	}

	for (int l = 0; l < nNoofLinkPerc; l++)
	{
		if (linkperc[l].m_eMedium1 == MED) linkperc[l].m_nGIndex1 = index_change_table[linkperc[l].m_nGIndex1];
		if (linkperc[l].m_eMedium2 == MED) linkperc[l].m_nGIndex2 = index_change_table[linkperc[l].m_nGIndex2];
	}

	delete[] index_change_table;
}


bool CPERC3D::CompareDistanceCutoff(int i, int j, double distance)
{
	double dR1 = m_cCNT[i].m_dRadius + (m_cCNT[i].m_nNoofWall - 1)*dDistance_VanderWaals;
	double dR2 = m_cCNT[j].m_dRadius + (m_cCNT[j].m_nNoofWall - 1)*dDistance_VanderWaals;
	if (distance <= dR1 + dR2 + m_dDistanceCutoff)
		return true;
	else
		return false;
}

double CPERC3D::DistanceSegToSeg(double x1i, double y1i, double z1i,
	double x1f, double y1f, double z1f,
	double x2i, double y2i, double z2i,
	double x2f, double y2f, double z2f,
	double *u, double *v)
{
	double P1x = x1f - x1i;
	double P1y = y1f - y1i;
	double P1z = z1f - z1i;

	double P2x = x2f - x2i;
	double P2y = y2f - y2i;
	double P2z = z2f - z2i;

	double a11 = P1x * P1x + P1y * P1y + P1z * P1z;
	double a12 = -(P1x * P2x + P1y * P2y + P1z * P2z);
	double a21 = (P1x * P2x + P1y * P2y + P1z * P2z);
	double a22 = -(P2x * P2x + P2y * P2y + P2z * P2z);

	double b1 = P1x * (x2i - x1i) + P1y * (y2i - y1i) + P1z * (z2i - z1i);
	double b2 = P2x * (x2i - x1i) + P2y * (y2i - y1i) + P2z * (z2i - z1i);
	double det = a11 * a22 - a12 * a21;

	if (det != 0.0)
	{
		*u = (a22 * b1 - a12 * b2) / det;
		*v = (-a21 * b1 + a11 * b2) / det;
	}

	if (*u <= 0) *u = 0.0;
	if (*u >= 1) *u = 1.0;

	if (*v <= 0) *v = 0.0;
	if (*v >= 1) *v = 1.0;

	double dx = *v * (x2f - x2i) - *u * (x1f - x1i) + (x2i - x1i);
	double dy = *v * (y2f - y2i) - *u * (y1f - y1i) + (y2i - y1i);
	double dz = *v * (z2f - z2i) - *u * (z1f - z1i) + (z2i - z1i);

	double distance = sqrt(dx*dx + dy * dy + dz * dz);

	return distance;
}



double CPERC3D::Gcontact(double d_um, double D_um)
{
	double G = 0;
	double dTransProbability = 0.0;
	double dtunnel = 1.0e6*(dPlankConst / (2.*M_PI)) / sqrt(8.0 * dElectronMess*m_dDeltaE*dElectronCharge);
	

	if (d_um <= (D_um + dDistance_VanderWaals))
	{
		dTransProbability = exp(-dDistance_VanderWaals / dtunnel);
	}
	else if ((D_um + dDistance_VanderWaals) < d_um && d_um <= (D_um + m_dDistanceCutoff))
	{
		dTransProbability = exp(-(d_um - D_um) / dtunnel);
	}
	else
	{
		dTransProbability = 0.0;
	}
	return dG0 * m_nNoofConductionChannelCNT*dTransProbability;
}

void CPERC3D::FilePrintCNTs(char *fn)
{
	FILE *fpcnt = fopen(fn, "w");

	for (int i = 0; i < m_nNoofCNTs; i++)
	{
		m_cCNT[i].PrintNodalInfo(fpcnt);
	}

	fclose(fpcnt);

}

void CPERC3D::FilePrintCNTs(char *fn, SInput sInParam, int *nPercolationList, int nNoofPercolationNetwork )
{

	FILE *fpcnt = fopen(fn, "w");

	FilePrintInput(fpcnt, sInParam);
	fprintf(fpcnt, "Total No. of Perc. Network is %d\n", nNoofPercolationNetwork);

	for (int n = 0; n < nNoofPercolationNetwork; n++)
	{
		fprintf(fpcnt, "%d\n", nPercolationList[n]);
	}
	fprintf(fpcnt, "\n");

	fprintf(fpcnt, "#CNT Data \n\n");

	for (int i = 0; i < m_nNoofCNTs; i++)
	{
		m_cCNT[i].PrintNodalInfo(fpcnt);
	}

	fclose(fpcnt);
}

void CPERC3D::FilePrintInput(FILE *fp, SInput sInParam)
{
	if (sInParam.m_eGeom == BOX) fprintf(fp, "DOMAIN: BOX\n");
	if (sInParam.m_eGeom == CYL) fprintf(fp, "DOMAIN: CYL\n");
	if (sInParam.m_eBoundaryCond == FBC) fprintf(fp, "BNDC: FBC\n");
	if (sInParam.m_eBoundaryCond == PBC) fprintf(fp, "BNDC: PBC\n");

	fprintf(fp, "MIN:(%3.3f, %3.3f, %3.3f), MAX:(%3.3f, %3.3f, %3.3f)\n", sInParam.m_dXmin, sInParam.m_dYmin, sInParam.m_dZmin, sInParam.m_dXmax, sInParam.m_dYmax, sInParam.m_dZmax);
	fprintf(fp, "No. of CNTs: %d\n", sInParam.m_nNoofCNTs);
	fprintf(fp, "Strain: %1.3f\n", sInParam.m_dStrain);
	fprintf(fp, "Poisson ratio: %1.3f\n", sInParam.m_dPoissonPolymer);
}

void CPERC3D::FilePrintCNTs(char *fn, int nNoofCNTsPerc, int *plist)
{
	FILE *fpcnt = fopen(fn, "w");
	for (int i = 0; i < nNoofCNTsPerc; i++)
	{
		m_cCNT[plist[i]].PrintNodalInfo(fpcnt);
	}
	fclose(fpcnt);
}

void CPERC3D::FilePrintLink(char *fn, SLinkedCond *link, int nNoofLink)
{
	FILE *fptxt_w = fopen(fn, "wt");

	for (int l = 0; l < nNoofLink; l++)
	{
		int i = link[l].m_nCNTNo1;
		int j = link[l].m_nCNTNo2;

		double g = link[l].m_dConductance;
		double dist = link[l].m_dDistance;
		int nG1 = link[l].m_nGIndex1;
		int nG2 = link[l].m_nGIndex2;
		double U = link[l].m_dU;
		double V = link[l].m_dV;
		eMEDIUM eMed1 = link[l].m_eMedium1;
		eMEDIUM eMed2 = link[l].m_eMedium2;

		char s1[255], s2[255];
		if (eMed1 == MED) sprintf(s1, "MED");
		if (eMed1 == TPE) sprintf(s1, "TPE");
		if (eMed1 == BTE) sprintf(s1, "BTE");

		if (eMed2 == MED) sprintf(s2, "MED");
		if (eMed2 == TPE) sprintf(s2, "TPE");
		if (eMed2 == BTE) sprintf(s2, "BTE");

		fprintf(fptxt_w, "[%7d] (%5d)\t%f\t(%5d)\t%f\t%f\t%f\t%5d[%s]\t%5d[%s]\n", l, i, U, j, V, dist, g, nG1, s1, nG2, s2);
	}

	fclose(fptxt_w);


}

void CPERC3D::FilePrintMatrix(char *fn, MatrixXd A)
{
	int row = A.rows();
	int col = A.cols();

	FILE *fp = fopen(fn, "w");

	for (int r = 0; r < row; r++)
	{
		for (int c = 0; c < col; c++)
		{
			fprintf(fp, "%1.12g\t", A(r, c));
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void CPERC3D::FilePrintVector(char *fn, VectorXd A)
{
	int row = A.rows();

	FILE *fp = fopen(fn, "w");

	for (int r = 0; r < row; r++)
	{
		fprintf(fp, "%1.12g\n", A(r));
	}

	fclose(fp);
}

