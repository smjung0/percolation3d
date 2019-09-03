//	perc3d.h													
//	Source code for the simulation of electrical percolation of	CNT/polymer composite
//	developed by Sungmin Jung, Ph.D.
//  University of Cambridge
//  Last update on 18th April 2019
//  Contact: sm0104.jung@gmail.com, sj569@cam.ac.uk, s.jung@eng.cam.ac.uk

#pragma once
#include "cnt.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include<iostream>
#include <algorithm>
#include <stdlib.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

#ifndef PERC3D_H
#define PERC3D_H 1

struct SInput
{
	eGEOM m_eGeom;// = BOX/CYL;
	eBCND m_eBoundaryCond;// = FBC/PBC;

	double m_dXmin;// = -10.0; // in um
	double m_dXmax;// = +10.0; // in um 
	double m_dYmin;// = -10.0; // in um
	double m_dYmax;// = +10.0; // in um
	double m_dZmin;// = -10.0; // in um
	double m_dZmax;// = +10.0; // in um

	// Physical Parameters on CNT/Polymer Network
	int m_nNoofCNTs; // number of CNTs in polymer

	bool m_bWeibull;// = true;
	double m_dLengthBar;// = 5.64; // um // Weibull Distribution scale
	double m_dLengthPower;// = 2.4;          // Weibull Distribution shape parameter
	int m_nRandSeedWeibull;
	bool m_bLognormal; // = true;
	double m_dMean;  // 3.2 for mean of 25nm length, 
	double m_dStdev; // deviation of lognormal, 0.0 ~0 .5 
	int m_nRandSeedLogNormal;

	int m_nNoofWall;// l = 1;
	double m_dConductivity_cnt;// = 1.0e-2; //[S/um]

	int m_nNoofNodesInitial;// = 2;
	double m_dTheta_dev_dgr;// = 0.;

	double m_dDensityPolymer;// = 1.25e-12; // 1.25e-12g/um^3
	double m_dPoissonPolymer;// = 0.5;
	double m_dAppliedVoltage;// = 1.0;

	int m_nNoofConductionChannel;// = 400; // ea
	double m_dDeltaE;// = 5.0; // eV
	double m_dDistanceCutoff;// = 0.0008; // um

	// numerical parameter for solving matrix
	double m_dTolerance; // = 1.0e-8;

	// strain parameter
	double m_dStrain; // = 0.00;
};


class CPERC3D
{
	public:

		eGEOM m_eGeom;
		eBCND m_eBC;

		double m_dXmax;
		double m_dXmin;
		double m_dYmax;
		double m_dYmin;
		double m_dZmax;
		double m_dZmin;

		double m_dDensityPolymer;

		double m_dTolerance;

		int m_nNoofConductionChannelCNT;
		double m_dDeltaE;
		double m_dDistanceCutoff;
		double m_dAppliedVoltage;

		CCNT *m_cCNT;
		int m_nNoofCNTs;


		// constructor and destructor in perc3d.cpp 
		CPERC3D();
		~CPERC3D();

		void CopyCNTs(CCNT *cnt, int nNoofCNTs);
		void SetDomain(eGEOM eGeom, eBCND eBC, double dXmax, double dXmin, double dYmax, double dYmin, double dZmax, double dZmin);
		void SetMaterialParam(double dDensityPolymer, double dAppliedVoltage, int m_nNoofConductionChannelCNT, double m_dDeltaE, double m_dDistanceCutoff);
		void WeightAndVolumeFraction(double *dWP, double *dVP);

		void Stretching(double dPoissonPolymer, double dPoissonCNT, double dStrainPolymerLong, double *dXmax, double *dXmin, double *dYmax, double *dYmin, double *dZmax, double *dZmin);
		void Stretching2(double dPoissonPolymer, double dStrainPolymerLong, double *dXmax, double *dXmin, double *dYmax, double *dYmin, double *dZmax, double *dZmin);
		void Stretching3(double dPoissonPolymer, double dStrainPolymerLong);


		void SetTolerance(double dTolerance) { m_dTolerance = dTolerance;  }

		// Finding Percolation Network in perc3d_structure.cpp
		int UnionCNTs(void);
			int UnionTwoCNTs(FILE *fpbin, int i, int j);
			double UnionCheck(int i, int j, int *pmin, int *qmin, double *umin, double *vmin);
				void Union(int i, int j);
				int Find(int i);
		void FindingGroup(int *glist);
		void ListingCNTsInGroup(int *glist, bool *b_perc);

		// Node Numbering for global index of matrix in perc3d.cpp
		void GlobalNodeNumbering(int nNoofCNTsPerc, int *plist, int *nNoofTotalNodePerc);
		void InsertNodesOnPercNetwork(SLinkedCond *linkperc, int nNoofLinkPerc, int *nNoofTotalNodePerc);
		void GlobalNodeNumberingReOrder(int nNoofCNTsPerc, int *plist, int nNoofTotalNodePerc, int nNoofLinkPerc, SLinkedCond *linkperc);

		// Solving global conductivity in perc3d_solve.cpp
		double SolveConductivityGlobal(int *nPercolationList, int *nNoofPercolationNetwork, double *dConductance, bool *error);
			double SolveConductanceIndPercNetwork(int Index, int *glist, int nNoofLink, bool *error);
				void SetMatrices(int nNoofPercCNTs, int *plist, SLinkedCond *slinkcond, int nNoofLinkedCond,
								SpMat *A, VectorXd *B, int nNoofTotalNodePerc);
				double ConductanceIndPercNetwork(int nNoofPercCNTs, int *plist, VectorXd mV, bool *error);

	
		void SetVoltagesToCNT(int nNoofPercCNTs, int *plist, VectorXd mV);
		void SetVoltagesToLink(int nNoofLinkedCond, SLinkedCond *slinkedcond, VectorXd mV);
		double CalPowerDissipationCNT(int nNoofPercCNT, int *plist);
		double CalPowerDissipationLink(int nNoofLink, SLinkedCond *slinkedcond);

		// Physical Calculation for distance in perc3d.cpp
		bool CompareDistanceCutoff(int i, int j, double distance);
		double DistanceSegToSeg(double x1i, double y1i, double z1i, double x1f, double y1f, double z1f,
								double x2i, double y2i, double z2i, double x2f, double y2f, double z2f, double *u, double *v);
		double Gcontact(double d_um, double D_um);

		// File out for debugging in perc3d.cpp
		void FilePrintMatrix(char *fn, MatrixXd A);
		void FilePrintVector(char *fn, VectorXd A);
		void FilePrintCNTs(char *fn);
		void FilePrintCNTs(char *fn, int nNoofCNTsPerc, int *plist);
		void FilePrintCNTs(char *fn, SInput sInParam, int *nPercolationList, int nNoofPercolationNetwork);
		void FilePrintLink(char *fn, SLinkedCond *link, int nNoofLink); // think later


		void FilePrintInput(FILE *fp, SInput sInParam);


};


#endif