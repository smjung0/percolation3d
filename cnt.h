//	cnt.h													
//	Source code for the simulation of electrical percolation of	CNT/polymer composite
//	developed by Sungmin Jung, Ph.D. 
//  University of Cambridge
//  Contact: sm0104.jung@gmail.com, sj569@cam.ac.uk, s.jung@eng.cam.ac.uk
//  Copyright (C) 2019


#pragma once

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
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

#ifndef CNT_H
#define CNT_H 1

//const double M_PI = acos(-1.);
const double dDistance_VanderWaals = 0.00034;
const double dPlankConst = 6.62607e-34; // m2 kg / s
const double dElectronMess = 9.10938356e-31; // kg
const double dElectronCharge = 1.60217646e-19;  // C
const double dSurfaceDensityCNT = (1.0 / 1315.0)*1.0e-12;
const double dG0 = 2.*dElectronCharge*dElectronCharge / dPlankConst;//0.774626e-4;

enum eMEDIUM { MED = 0, TPE, BTE }; // medium of node or CNTs
enum eGEOM { BOX = 0, CYL };
enum eBCND { FBC = 0, PBC };
enum eCALC { ONCE = 0, VOLF, STRN };

struct nLoopParam
{
	int nInit;
	int nNoofSteps;
	int nDelta;
};

struct dLoopParam
{
	double dInit;
	int nNoofSteps;
	double dDelta;
};

struct SLinkedCond
{
	int m_nCNTNo1 = 0;
	int m_nCNTNo2 = 0;

	double m_dU = 0.0;
	double m_dV = 0.0;

	int m_nGIndex1 = 0;
	int m_nGIndex2 = 0;
	eMEDIUM m_eMedium1 = MED;
	eMEDIUM m_eMedium2 = MED;

	double m_dDistance = 0.0;
	double m_dConductance = 0.0;

	double m_dV1 = 0.0;
	double m_dV2 = 0.0;
	double m_dP = 0.0;
};

class CCNT
{

public:
	int m_nParent;  // parent adress for union find algorithm
	int m_nAddress; // its own address
	int m_nRank;    // level in the network

	int m_nNoofNodes; // number of nodes on the single CNT
	Vector3d *m_vP; // position vector

	double *m_dL; // position along the line
	double *m_dV; // voltage at node

	eMEDIUM *m_eNodalProperty;
	eMEDIUM m_eCNTProperty;
	int *m_nGIndexList;

	double m_dXmax, m_dXmin;
	double m_dYmax, m_dYmin;
	double m_dZmax, m_dZmin;

	double m_dXcntmax, m_dXcntmin;
	double m_dYcntmax, m_dYcntmin;
	double m_dZcntmax, m_dZcntmin;

	Vector3d m_vPc;

	double m_dLength; // length of CNT
	double m_dRadius; // radius of CNT
	int m_nNoofWall;  // number of wall of CNT

	double m_dConductivity; 
	double *m_dConductanceList;

	double *m_dP; // powe at segment
//	double *m_dR; // radius of segment

public:
	CCNT();
	~CCNT();

	void SetNoofWall(int nNoofWall);
	double WeightCNT(void);
	double VolumeCNT(void);
	double dtor(double dgr) { return dgr * M_PI / 180.; }
	double rtod(double rad) { return rad * 180. / M_PI; }

	void MaksSet(int index);
	void SetDomain(double dXmin, double dXmax, double dYmin, double dYmax, double dZmin, double dZmax);
	void SetInitLineSegments(eGEOM eGeom, int nNoofNodes, double dLength, double dRadius, double dtheta_dev_dgr);
	void SetCNTBoundary(void);
	void SetNodalProperties(void);

	void PrintNodalInfo(FILE *fp);
	int InsertNode(int index, double u, int *nGI, eMEDIUM *eMed);
	int InsertNode(double L, int *nGI, eMEDIUM *eMed);
	void Stretching(double dPoissionPolymer, double dPoissionCNT, double dStrainLong);
	void Stretching2(double dPoissionPolymer, double dStrainLong);
	void SetCNTRange(void);
	eMEDIUM ReturnNodalProperty(double z);

	void SetConductance(double dConductivity);
//	void SetConductance(void);
	void SetVoltage(double dAppliedVoltage, VectorXd *mV);
	double SetPowerDissipation(void);

	CCNT& operator = (const CCNT& x);

	void GetAngles(Vector3d n, double *theta, double *phi);
	Vector3d GetUnitVector(double theta, double phi);
	void AddAngles(double *theta, double *phi, double dtheta, double dphi);
};

#endif
