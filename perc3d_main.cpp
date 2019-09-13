//	perc3d_main.cpp													
//	Source code for the simulation of electrical percolation of	CNT/polymer composite
//	developed by Sungmin Jung, Ph.D.
//  University of Cambridge
//  Last update on 29th May 2019
//  Contact: sm0104.jung@gmail.com, sj569@cam.ac.uk, s.jung@eng.cam.ac.uk

#include "perc3d.h"
#include "cnt.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <iostream>
#include <random>
#include <chrono>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

void ConfigureCNTsAndDomain(eGEOM eGeom, int nRandomSeed,
	CCNT *CNTs, int nNoofCNTs, 
	double dXMin, double dXMax, double dYMin, double dYMax, double dZMin, double dZMax,
	double dConductivity_cnt, int nNoofNodesInitial,
	double *dLength, double *dRadius, int nNoofWall, double dtheta_dev_dgr, 
	double dStrainPLM, double dPoissonRatioPLM);

void WeibullDistribution(bool bWeibul, int nRandSeedWeibull, double *dLengthCNT, int nNoofCNTs, double dPower, double dLengthBar);
void LogNormDistribution(bool bLogNorm, int nRandSeedLogNormal, double *dRadiusCNT, int nNoofCNTs, double dMean, double dStdev);

void ReadInput(char *fn, SInput *sInParam, eCALC *eCalculate, int *nRandomSeed, int *nRandTrial, nLoopParam *nNCNT, dLoopParam *dSTRN, bool *bSaveCNTData );
void SingleCalculation(SInput sInParam, int nRandomSeed, double *dVolumeFraction, double *dWeightFraction, double *dConductance, double *dCondcutivity, bool *bError, char *fn_cnt_out, bool bSaveCNTData);

int main(int argc, char *argv[]) 
{
	int year = 2019;
	int date = 12;
	char month[255] = "SEP";
	double version = 0.44;

	
	if ( argv[1] == NULL )
	{
		printf("Version %2.2f: %2d %s, %d\n", version, date, month, year);
		exit(-1);
	}

	//////////////////////////////////////////////////////////////////////////
	//																		//
	// Parameter Setting for Simulation										//
	//																		//
	//////////////////////////////////////////////////////////////////////////

	char fninput[255] = "input.txt";
	SInput sInParam;
	eCALC eCalculate;
	int nRandTrial;
	int nRandomSeed;
	nLoopParam nNCNT; 
	dLoopParam dSTRN;
	bool bSaveCNTData = false;

	ReadInput(/*fninput*/argv[1], &sInParam, &eCalculate, &nRandomSeed, &nRandTrial, &nNCNT, &dSTRN, &bSaveCNTData);

//	omp_set_num_threads(28);


	int *nRandSeedArray = new int[nRandTrial];
	if (nRandomSeed < 0)
		for (int i = 0; i < nRandTrial; i++) 	nRandSeedArray[i] = (int) (RAND_MAX*(double) (i+1)/(double) nRandTrial);
	else 
		for (int i = 0; i < nRandTrial; i++) 	nRandSeedArray[i] = nRandomSeed;

	char directory[256] = "./cntdata";


	if (eCalculate == VOLF)
	{
		double dVolumeFraction;
		double dWeightFraction;
		double dConductance;
		double dConductivity;
		bool bError;

		FILE *fp_conductance;
		fp_conductance = fopen("conductance_VOLF.txt", "w");
		FILE *fp_conductivity;
		fp_conductivity = fopen("conductivity_VOLF.txt", "w");
		
		if( bSaveCNTData )
		{
			mkdir(directory, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}

		fprintf(fp_conductance, "VF[vol%%]\tWF[wt%%]");
		fprintf(fp_conductivity, "VF[vol%%]\tWF[wt%%]");
		for (int r = 0; r < nRandTrial; r++)
		{
			fprintf(fp_conductance, "\tTRY(%d)", r + 1);
			fprintf(fp_conductivity, "\tTRY(%d)", r + 1);
		}
		fprintf(fp_conductance, "\n");
		fprintf(fp_conductivity, "\n");

		for (int n = 0; n < nNCNT.nNoofSteps; n++)
		{
			sInParam.m_nNoofCNTs = nNCNT.nDelta*n + nNCNT.nInit;

			for (int r = 0; r < nRandTrial; r++)
			{

				if (nRandSeedArray[r] == 0)
				{
					sleep(1000); //
					nRandSeedArray[r] = time(NULL);
				}

				printf("Number of trial: %d of %d\n", r+1, nRandTrial);
				printf("Random Seed: %d\n", nRandSeedArray[r]);

				char fn_cnt_out[255];

				if(bSaveCNTData)
				{
					sprintf(fn_cnt_out, "cnt_n%d_s%2.2f_r%d.txt", sInParam.m_nNoofCNTs, sInParam.m_dStrain, nRandSeedArray[r]);
					chdir(directory);
				}
				SingleCalculation(sInParam, nRandSeedArray[r], &dVolumeFraction, &dWeightFraction, &dConductance, &dConductivity, &bError, fn_cnt_out, bSaveCNTData);
				if(bSaveCNTData)
				{
					chdir("..");
				}

				printf("\n");

				if (r == 0) fprintf(fp_conductance, "%2.5e\t%2.5e\t%2.5e", dVolumeFraction*100., dWeightFraction * 100, dConductance);
				else fprintf(fp_conductance, "\t%5.2e", dConductance);

				if (r == 0) fprintf(fp_conductivity, "%2.5e\t%2.5e\t%2.5e", dVolumeFraction*100., dWeightFraction * 100, dConductivity);
				else fprintf(fp_conductivity, "\t%5.2e", dConductivity);

			}
			fprintf(fp_conductance, "\n");
			fprintf(fp_conductivity, "\n");

		}

		fclose(fp_conductivity);
		fclose(fp_conductance);

	}

	if (eCalculate == STRN)
	{
		double dVolumeFraction;
		double dWeightFraction;
		double dConductance;
		double dConductivity;
		bool bError;

		FILE *fp_conductance = fopen("conductance_STRN.txt", "w");
		FILE *fp_conductivity = fopen("conductivity_STRN.txt", "w");
		
		if( bSaveCNTData )
		{
			mkdir(directory,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}

		fprintf(fp_conductance, "STRAIN[%%]\tVF[vol%%]\tWF[wt%%]");
		fprintf(fp_conductivity, "STRAIN[%%]\tVF[vol%%]\tWF[wt%%]");
		for (int r = 0; r < nRandTrial; r++)
		{
			fprintf(fp_conductance, "\tTRY(%d)", r + 1);
			fprintf(fp_conductivity, "\tTRY(%d)", r + 1);
		}
		fprintf(fp_conductance, "\n");
		fprintf(fp_conductivity, "\n");



		for (int s = 0; s < dSTRN.nNoofSteps; s++)
		{
			sInParam.m_dStrain = dSTRN.dInit + dSTRN.dDelta*s;


			for (int r = 0; r < nRandTrial; r++)
			{
				printf("Strain : %f\n", sInParam.m_dStrain);

				if (nRandSeedArray[r] == 0)
				{
					sleep(1000); //
					nRandSeedArray[r] = time(NULL);
				}

				printf("Number of trial: %d of %d\n", r+1, nRandTrial);
				printf("Random Seed: %d\n", nRandSeedArray[r]);

				char fn_cnt_out[256];

				if( bSaveCNTData ) 
				{
					sprintf(fn_cnt_out, "cnt_n%d_s%2.2f_r%d.txt", sInParam.m_nNoofCNTs, sInParam.m_dStrain, nRandSeedArray[r]);
					chdir(directory);
				}	
				SingleCalculation(sInParam, nRandSeedArray[r], &dVolumeFraction, &dWeightFraction, &dConductance, &dConductivity, &bError, fn_cnt_out, bSaveCNTData );
				
				if( bSaveCNTData )
				{
					chdir("..");
				}
				printf("\n");

				if (r == 0) fprintf(fp_conductance, "%2.5e\t%2.5e\t%2.5e\t%2.5e", sInParam.m_dStrain, dVolumeFraction*100., dWeightFraction*100., dConductance);
				else fprintf(fp_conductance, "\t%5.2e", dConductance);

				if (r == 0) fprintf(fp_conductivity, "%2.5e\t%2.5e\t%2.5e\t%2.5e", sInParam.m_dStrain, dVolumeFraction*100., dWeightFraction*100., dConductivity);
				else fprintf(fp_conductivity, "\t%5.2e", dConductivity);

			}
			fprintf(fp_conductance, "\n");
			fprintf(fp_conductivity, "\n");

		}

		fclose(fp_conductivity);
		fclose(fp_conductance);

	}


	delete[] nRandSeedArray;

	return 1;
}

void WeibullDistribution(bool bWeibul, int nRandSeedWeibull, double *dLengthCNT, int nNoofCNTs, double dPower, double dLengthBar)
{
	unsigned nRandSeedforLength;

	if (nRandSeedWeibull <= 0)
		nRandSeedforLength = std::chrono::system_clock::now().time_since_epoch().count();

	else 
		nRandSeedforLength = nRandSeedWeibull;
	
	std::default_random_engine generator(nRandSeedforLength);
	std::weibull_distribution<double> distribution(dPower, dLengthBar);

	int p[20] = {};

	for (int i = 0; i < nNoofCNTs; i++)
	{
		if( bWeibul == true )	dLengthCNT[i] = distribution(generator);
		else dLengthCNT[i] = dLengthBar;

		if (dLengthCNT[i] < 20) ++p[(int)dLengthCNT[i]];
	}

	FILE *fp = fopen("weibull.txt", "w");
	for (int n = 0; n < 20; n++)
	{
		fprintf(fp, "%f\t%d\n", n + 0.5, p[n]);
	}
	fclose(fp);

}

void LogNormDistribution(bool bLogNorm, int nRandSeedLogNormal, double *dRadiusCNT, int nNoofCNTs, double dMean, double dStdev)
{
	unsigned nRandSeedforRadius;

	if (nRandSeedLogNormal <= 0)
		nRandSeedforRadius = std::chrono::system_clock::now().time_since_epoch().count();
	else
		nRandSeedforRadius = nRandSeedLogNormal;

	std::default_random_engine generator(nRandSeedforRadius);
	std::lognormal_distribution<double> distribution(dMean, dStdev);

	int p[200] = {};

	for (int i = 0; i < nNoofCNTs; i++)
	{
		if (bLogNorm == true)
		{
			dRadiusCNT[i] = distribution(generator);
			if (dRadiusCNT[i] < 200) ++p[(int)(dRadiusCNT[i])];
			dRadiusCNT[i] = dRadiusCNT[i] / 1000.;
		}

		else dRadiusCNT[i] = dMean;
	}

	FILE *fp = fopen("lognormal.txt", "w");
	for (int n = 0; n < 200; n++)
	{
		fprintf(fp, "%f\t%d\n", (n + 0.5), p[n]);
	}
	fclose(fp);

}


void ConfigureCNTsAndDomain(eGEOM eGeom, int nRandomSeed,
	CCNT *CNTs, int nNoofCNTs,
	double dXMin, double dXMax, double dYMin, double dYMax, double dZMin, double dZMax,
	double dConductivity_cnt, int nNoofNodesInitial,
	double *dLength, double *dRadius, int nNoofWall, double dtheta_dev_dgr, 
	double dStrainPLM, double dPoissonRatioPLM)
{
	srand(nRandomSeed);

	std::cout << "Initializing CNTs... ";

	int i;
//#pragma omp parallel for private(i) //schedule(dynamic)
	for (i = 0; i < nNoofCNTs; i++)
	{
		CNTs[i].MaksSet(i);
		CNTs[i].SetDomain(dXMin, dXMax, dYMin, dYMax, dZMin, dZMax);
		CNTs[i].SetInitLineSegments(eGeom, nNoofNodesInitial, dLength[i], dRadius[i], dtheta_dev_dgr);
		CNTs[i].Stretching2(dPoissonRatioPLM, dStrainPLM);
		CNTs[i].SetNoofWall(nNoofWall);
		CNTs[i].SetCNTBoundary();
		CNTs[i].SetNodalProperties();
		CNTs[i].SetConductance(dConductivity_cnt);
		CNTs[i].SetCNTRange();
	}

	std::cout << "Finished!\n";
}


void ReadInput(char *fn, SInput *sInParam, eCALC *eCalculate, int *nRandomSeed, int *nRandTrial, nLoopParam *nNCNT, dLoopParam *dSTRN, bool *bSaveCNTData)
{
	FILE *fp = fopen(fn, "r");

	char dummy[255], dummy2[255], dummy3[255];
	char str[200] = { 0, };
	string sstr;

	string keywd[] = {	"STRUCTURE", "DIMENSION", "VOLTAGE", "LENGTH", "RADIUS", "TOLERANCE", 
						"CNT", "POLYMER" , "RANDTRIAL", "CALCULATE", "END" };

	char *KEYCK;
	char *CMTCK;
	char *line;
	int seq;

	do {
		line = fgets(str, sizeof(str), fp);

		seq = 100;
		CMTCK = strstr(str, "#");
		if (CMTCK == NULL)
		{
			for (seq = 0; seq < 11; seq++)
			{
				sstr = str;
				if (sstr.find(keywd[seq]) != string::npos)
				{
					break;
				}
			}

		}

		switch (seq)
		{

		case 0:
			// STRUCTURE
			sscanf(str, "%s %s %s", dummy, dummy2, dummy3);
			if (strstr(dummy2, "BOX")) sInParam->m_eGeom = BOX;
			else sInParam->m_eGeom = CYL;
			if (strstr(dummy3, "FBC")) sInParam->m_eBoundaryCond = FBC;
			else sInParam->m_eBoundaryCond = PBC;
			break;

		case 1:
			// DIMENSION
			sscanf(str, "%s x(%lf [um], %lf [um]) y(%lf [um], %lf [um]) z(%lf [um], %lf [um])", dummy,	&(sInParam->m_dXmin), &(sInParam->m_dXmax),
																										&(sInParam->m_dYmin), &(sInParam->m_dYmax),
																										&(sInParam->m_dZmin), &(sInParam->m_dZmax) );
			break;

		case 2:
			// VOLTAGE
			sscanf(str, "%s %lf [V]", dummy, &(sInParam->m_dAppliedVoltage));
			break;

		case 3:
			// LENGTH 
			// LENGTH WEIBULL 5.64 [um] 2.4
			if (strstr(str, "WEIBULL")) {

				sscanf(str, "%s %s %d %lf %lf", dummy, dummy2, &(sInParam->m_nRandSeedWeibull), &(sInParam->m_dLengthBar), &(sInParam->m_dLengthPower) );
				sInParam->m_bWeibull = true;
			}
			// LENGTH CONST 5.64 [um]
			else
			{
				sscanf(str, "%s %s %lf", dummy, dummy2, &(sInParam->m_dLengthBar) );
				sInParam->m_dLengthPower = 0.0;
				sInParam->m_bWeibull = false;
				sInParam->m_nRandSeedWeibull = -1;
			}
			break;

		case 4:
			// RADIUS 
			//RADIUS LOGNORMAL 0.025[um] 0.3363[um]
			if (strstr(str, "LOGNORMAL"))
			{
				sscanf(str, "%s %s %d %lf %lf", dummy, dummy2, &(sInParam->m_nRandSeedLogNormal), &(sInParam->m_dMean), &(sInParam->m_dStdev) );
				sInParam->m_bLognormal = true;
			}
			//RADIUS CONST
			else
			{
				sscanf(str, "%s %s %lf", dummy, dummy2, &(sInParam->m_dMean) );
				sInParam->m_dStdev = 1.0;
				sInParam->m_bLognormal = false;
				sInParam->m_nRandSeedLogNormal = -1;
			}
			break;

		case 5:
			// TOLERANCE
			sscanf(str, "%s %lf", dummy, &(sInParam->m_dTolerance) );
			break;

		case 6:
			// CNT
			sscanf(str, "%s save=%s noofcnt=%d wall=%d sigma=%lf [S/um] node=%d tht_max=%lf [dgr]", dummy, dummy2, &(sInParam->m_nNoofCNTs), &(sInParam->m_nNoofWall), &(sInParam->m_dConductivity_cnt), &(sInParam->m_nNoofNodesInitial), &(sInParam->m_dTheta_dev_dgr) );
			if (strstr(dummy2, "yes")|| strstr(dummy2, "YES") || strstr(dummy2, "Yes")) *bSaveCNTData = true;
			else if (strstr(dummy2, "no") || strstr(dummy2, "No") || strstr(dummy2, "NO")) *bSaveCNTData = false;
			else *bSaveCNTData = false;
			break;

		case 7:
			// POLYMER
			sscanf(str, "%s strain=%lf mass=%lf [g/um^3] poisson=%lf deltaE=%lf [eV] M=%d Dcutoff=%lf [um]", dummy, &(sInParam->m_dStrain), &(sInParam->m_dDensityPolymer), &(sInParam->m_dPoissonPolymer), &(sInParam->m_dDeltaE), &(sInParam->m_nNoofConductionChannel), &(sInParam->m_dDistanceCutoff) );
			break;

		case 8:
			// RANDTRIAL
			sscanf(str, "%s n=%d seed=%d", dummy, nRandTrial, nRandomSeed);
			break;

		case 9:
			// CALCULATE
			if (strstr(str, "VOLF"))
			{
				*eCalculate = VOLF;
				sscanf(str, "%s %s init=%d noofsteps=%d delta=%d", dummy, dummy2, &(nNCNT->nInit), &(nNCNT->nNoofSteps), &(nNCNT->nDelta) );
			}
			if (strstr(str, "STRN"))
			{
				*eCalculate = STRN;
				// double *dSInit, int *nSStep, double *dSDelta,
				sscanf(str, "%s %s init=%lf noofsteps=%d delta=%lf", dummy, dummy2, &(dSTRN->dInit), &(dSTRN->nNoofSteps), &(dSTRN->dDelta)  );
			}
			if (strstr(str, "ONCE"))
			{
				*eCalculate = ONCE;
			}
			break;

		case 10:
			// END
			break;

		default: break;
		}

	} while (line != NULL);

	fclose(fp);



}


void SingleCalculation(SInput sInParam, int nRandomSeed,
	double *dFractionVolume, double *dFractionWeight, 
	double *dConductance, double *dConductivity, 
	bool *bError, char *fn_cnt_out, bool bSaveCNTData)
{
	CCNT *CNTs = new CCNT[sInParam.m_nNoofCNTs];
	double *dLengthCNT = new double[sInParam.m_nNoofCNTs];
	double *dRadiusCNT = new double[sInParam.m_nNoofCNTs];

	std::cout << "No. of CNTs is " << sInParam.m_nNoofCNTs << "." << endl;

	WeibullDistribution(sInParam.m_bWeibull, sInParam.m_nRandSeedWeibull, dLengthCNT, sInParam.m_nNoofCNTs, sInParam.m_dLengthPower, sInParam.m_dLengthBar);
	LogNormDistribution(sInParam.m_bLognormal, sInParam.m_nRandSeedLogNormal, dRadiusCNT, sInParam.m_nNoofCNTs, sInParam.m_dMean, sInParam.m_dStdev);

	ConfigureCNTsAndDomain(sInParam.m_eGeom, nRandomSeed, CNTs, sInParam.m_nNoofCNTs,
		sInParam.m_dXmin, sInParam.m_dXmax, sInParam.m_dYmin, sInParam.m_dYmax, sInParam.m_dZmin, sInParam.m_dZmax,
		sInParam.m_dConductivity_cnt, sInParam.m_nNoofNodesInitial, dLengthCNT, dRadiusCNT, sInParam.m_nNoofWall, sInParam.m_dTheta_dev_dgr,
		sInParam.m_dStrain, sInParam.m_dPoissonPolymer);

	int *nPercolationList = new int[sInParam.m_nNoofCNTs];
	int nNoofPercolationNetwork = 0;

	CPERC3D cPerc3D;
	cPerc3D.SetDomain(sInParam.m_eGeom, sInParam.m_eBoundaryCond, sInParam.m_dXmax, sInParam.m_dXmin, sInParam.m_dYmax, sInParam.m_dYmin, sInParam.m_dZmax, sInParam.m_dZmin);
	cPerc3D.SetMaterialParam(sInParam.m_dDensityPolymer, sInParam.m_dAppliedVoltage, sInParam.m_nNoofConductionChannel, sInParam.m_dDeltaE, sInParam.m_dDistanceCutoff);
	cPerc3D.CopyCNTs(CNTs, sInParam.m_nNoofCNTs);
	cPerc3D.SetTolerance(sInParam.m_dTolerance);

	cPerc3D.Stretching3(sInParam.m_dPoissonPolymer, sInParam.m_dStrain);

	cPerc3D.WeightAndVolumeFraction(dFractionWeight, dFractionVolume);

	*dConductivity = cPerc3D.SolveConductivityGlobal(nPercolationList, &nNoofPercolationNetwork, dConductance, bError);

	if(bSaveCNTData)
	{
		cPerc3D.FilePrintCNTs(fn_cnt_out, sInParam, nPercolationList, nNoofPercolationNetwork);
	}

	delete[] CNTs;
	delete[] nPercolationList;
	delete[] dLengthCNT;
	delete[] dRadiusCNT;
}
