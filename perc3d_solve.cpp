//	perc3d_solve.cpp													
//	Source code for the simulation of electrical percolation of	CNT/polymer composite
//	developed by Sungmin Jung, Ph.D. 
//  University of Cambridge
//  Last update on 18th April 2019
//  Contact: sm0104.jung@gmail.com, sj569@cam.ac.uk, s.jung@eng.cam.ac.uk


#include "perc3d.h"

double CPERC3D::SolveConductivityGlobal(int *nPercolationList, int *nNoofPercolationNetwork, double *dConductance, bool *error)
{
	std::cout << "Union-Find Started... ";
	clock_t t;
	t = clock();
	int nNoofLink = UnionCNTs();


	int *glist = new int[m_nNoofCNTs];
	bool *b_perc = new bool[m_nNoofCNTs];

	std::cout << "Finished!\n";
	std::cout << "Finding Percolation Network Started... ";

	*nNoofPercolationNetwork = 0;
	FindingGroup(glist);
	ListingCNTsInGroup(glist, b_perc);

	t = clock() - t;

	std::cout << "Finished!\n";
	printf("- Calculation Time for Union-Find Algorithm is %f [sec].\n", (double)t / (double)CLOCKS_PER_SEC);

	for (int i = 0; i < m_nNoofCNTs; i++)
	{
		if (b_perc[i] == true)
		{
			nPercolationList[*nNoofPercolationNetwork] = i; // Individual Percolation network
			(*nNoofPercolationNetwork)++;
			std::cout << "- Representative CNT in Percolation Network: " << i << endl;
		}
	}

	double dGlobalConductance = 0;

	bool error_ind = true;
	
	for (int i = 0; i < *nNoofPercolationNetwork; i++)
	{
		dGlobalConductance += SolveConductanceIndPercNetwork(nPercolationList[i], glist, nNoofLink, &error_ind);
		if (error_ind == false) *error = false;
	}

	double dGlobalConductivity = dGlobalConductance * (m_dZmax - m_dZmin) / ((m_dXmax - m_dXmin)*(m_dYmax - m_dYmin));

	*dConductance = dGlobalConductance;
	printf("Conductance: %1.6e [S]\n", dGlobalConductance);
	printf("Conductivity: %1.6g [S/m]\n", dGlobalConductivity*1.0e6);
//	std::cout << endl;

	delete[] glist;
	delete[] b_perc;

	return dGlobalConductivity * 1.0e6;

}


double CPERC3D::SolveConductanceIndPercNetwork(int Index, int *glist, int nNoofLink, bool *error)
{
	// making percolation list
	int *plist;
	int nNoofCNTsPerc = glist[Index];
	plist = new int[nNoofCNTsPerc];

	int pindex = 0;
	for (int i = 0; i < m_nNoofCNTs; i++)
	{
		if (m_cCNT[i].m_nParent == Index)
		{
			plist[pindex] = i;
			pindex++;
		}
	}

	int nNoofTotalNodePerc = 0;
	int nNoofLinkPerc = 0;

	SLinkedCond *link, *linkperc;
	link = new SLinkedCond[nNoofLink];
	linkperc = new SLinkedCond[nNoofLink];

	FILE *fpbin_r = fopen("link_total.bin", "rb");
	int lp = 0;
	for (int l = 0; l < nNoofLink; l++)
	{
		fread(&(link[l]), sizeof(link[l]), 1, fpbin_r);
		if (m_cCNT[link[l].m_nCNTNo1].m_nParent == Index && m_cCNT[link[l].m_nCNTNo2].m_nParent == Index)
		{
			linkperc[lp] = link[l];
			lp++;
		}
	}
	fclose(fpbin_r);

	nNoofLinkPerc = lp;

	printf("- No. of link in percolation network is %d.\n", nNoofLinkPerc);

	GlobalNodeNumbering(nNoofCNTsPerc, plist, &nNoofTotalNodePerc);
	printf("- No. of Total Nodes in Perc. Network before linking is: %d.\n", nNoofTotalNodePerc);

	InsertNodesOnPercNetwork(linkperc, nNoofLinkPerc, &nNoofTotalNodePerc);
	printf("- No. of Total Nodes in Perc. Network after  linking is: %d.\n", nNoofTotalNodePerc);

	printf("Global Node Number Setting... ");
	GlobalNodeNumberingReOrder(nNoofCNTsPerc, plist, nNoofTotalNodePerc, nNoofLinkPerc, linkperc);
	printf("Finished!\n");

	// debugging line
	char fn_link[255] = "link_perc.txt";
	FilePrintLink(fn_link, linkperc, nNoofLinkPerc);

	printf("Matrix Setting Started... ");
	clock_t t;
	t = clock();
	SpMat mA(nNoofTotalNodePerc, nNoofTotalNodePerc);
	VectorXd mB(nNoofTotalNodePerc);
	VectorXd mX(nNoofTotalNodePerc);

	mA.setZero();
	mB.setZero();
	mX.setZero();

	SetMatrices(nNoofCNTsPerc, plist, linkperc, nNoofLinkPerc, &mA, &mB, nNoofTotalNodePerc);
	printf("Finished!\n");
/*
	mA = mA;
	mB = mB;
*/
//	char fnA[255] = "A.txt";
//	char fnB[255] = "B.txt";
//	FilePrintMatrix(fnA, mA);
//	FilePrintVector(fnB, mB);

	printf("Solver Tolerance: %1.9g\n", m_dTolerance);

	printf("Matrix Solving Started... ");
	// LeastSquaresConjugateGradient<SparseMatrix<double> > solver;  
	// BiCGSTAB<SparseMatrix<double> > solver;
	// 	ConjugateGradient<SparseMatrix<double>, Lower|Upper, IdentityPreconditioner > solver; // works but Lower/Upper only parallel processing which has no meaning
	ConjugateGradient<SparseMatrix<double>, Lower | Upper > solver;
//	 	ConjugateGradient<SparseMatrix<double> > solver;
		// Eigen does not support full multiple processing for Conjugate Gradient
		// For the Conjugate Gradient in Eigen, Lower|Upper decomposition shows 
		// the best speed of 3.x sec for 20000 number of CNTs under the condition of 
		// tolerance to be 1.0e-9
		// 
	//	printf("Tol: %1.9g\n", solver.tolerance());

	int n_th = Eigen::nbThreads();
	//	printf("Parallel Processing Thread: %d\n", n_th);
	Eigen::setNbThreads(n_th);
	Eigen::initParallel();
	
	solver.setTolerance(m_dTolerance);

	solver.compute(mA);

	if (solver.info() != Success)
	{
		printf("\n- Decomposition Failed\n");
		exit(-1);
	}
	mX = solver.solve(mB);
	if (solver.info() != Success)
	{
		printf("\n- Solving Failed\n");
		exit(-1);
	}

	t = clock() - t;
	printf("Finished!\n");
	printf("- Calculation Time for Matrix Solving is %f [sec].\n", (double)t / CLOCKS_PER_SEC);

	printf("Conductance Calculation Startd\n");

	SetVoltagesToCNT(nNoofCNTsPerc, plist, mX);
	SetVoltagesToLink(nNoofLinkPerc, linkperc, mX);

	ConductanceIndPercNetwork(nNoofCNTsPerc, plist, mX, error);

	double dPowerCNTTotal = CalPowerDissipationCNT(nNoofCNTsPerc, plist);
	double dPowerLinkTotal = CalPowerDissipationLink(nNoofLink, linkperc);

	double dPowerTotal = (dPowerCNTTotal + dPowerLinkTotal);
	printf("- Total Power Dissipation in CNTs is %2.12e [J]\n", dPowerCNTTotal);
	printf("- Total Power Dissipation in LNKs is %2.12e [J]\n", dPowerLinkTotal);
	printf("- Total Power Dissipation is %2.12e [J]\n", dPowerTotal);

	//	double conductance = Conductance(nNoofCNTsPerc, plist, mX);
	double conductance = dPowerTotal / (m_dAppliedVoltage*m_dAppliedVoltage);
	

	delete[] plist;
	delete[] link;
	delete[] linkperc;

	return conductance;
}



void CPERC3D::SetMatrices(
	int nNoofPercCNTs, int *plist,
	SLinkedCond *slinkedcond, int nNoofLinkedCond,
	SpMat *mA, VectorXd *mB, int nNoofTotalNodePerc)
{
	std::vector<T> tripleList;
	tripleList.reserve(nNoofTotalNodePerc * 100);

//	std::cout << "Matrix Setting of Direct Conductance Terms\n";
	int i;

	for (i = 0; i < nNoofPercCNTs; i++)
	{
		int pID = plist[i];
		int nNodes = m_cCNT[pID].m_nNoofNodes;
		int P, Q, R;
/*
		if (pID == 3374)
		{
			printf("..");
		}
*/
		for (int m = 0; m < nNodes; m++)
		{
			Q = m_cCNT[pID].m_nGIndexList[m];

			if (m_cCNT[pID].m_eNodalProperty[m] == MED)
			{
				if (m == 0)
				{
					R = m_cCNT[pID].m_nGIndexList[m + 1];
					eMEDIUM eM = m_cCNT[pID].m_eNodalProperty[m + 1];
					double Cond = m_cCNT[pID].m_dConductanceList[m];
					if (eM == MED)
					{
						tripleList.push_back(T(Q, Q, -Cond));
						tripleList.push_back(T(Q, R, +Cond));
						(*mB)(Q) += 0.0;
					}
					if (eM == TPE)
					{
						tripleList.push_back(T(Q, Q, -Cond));
						(*mB)(Q) += -Cond * m_dAppliedVoltage;
					}
					if (eM == BTE)
					{
						tripleList.push_back(T(Q, Q, -Cond));
						(*mB)(Q) += 0.0;
					}
				}
				else if (1 <= m && m < nNodes - 1)
				{
					P = m_cCNT[pID].m_nGIndexList[m - 1];
					R = m_cCNT[pID].m_nGIndexList[m + 1];
					eMEDIUM eM1 = m_cCNT[pID].m_eNodalProperty[m - 1];
					eMEDIUM eM2 = m_cCNT[pID].m_eNodalProperty[m + 1];

					double Cond1 = m_cCNT[pID].m_dConductanceList[m - 1];
					double Cond2 = m_cCNT[pID].m_dConductanceList[m];

					if (eM1 == BTE)
					{
						tripleList.push_back(T(Q, Q, -Cond1));
						(*mB)(Q) += 0.0;
					}
					else if (eM1 == TPE)
					{
						tripleList.push_back(T(Q, Q, -Cond1));
						(*mB)(Q) += -Cond1 * m_dAppliedVoltage;
					}
					else
					{
						tripleList.push_back(T(Q, P, +Cond1));
						tripleList.push_back(T(Q, Q, -Cond1));
						(*mB)(Q) += 0.0;
					}

					if (eM2 == BTE)
					{
						tripleList.push_back(T(Q, Q, -Cond2));
						(*mB)(Q) += 0.0;
					}
					else if (eM2 == TPE)
					{
						tripleList.push_back(T(Q, Q, -Cond2));
						(*mB)(Q) += -Cond2 * m_dAppliedVoltage;
					}
					else
					{
						tripleList.push_back(T(Q, R, +Cond2));
						tripleList.push_back(T(Q, Q, -Cond2));
						(*mB)(Q) += 0.0;
					}

				}
				else
				{
					P = m_cCNT[pID].m_nGIndexList[m - 1];
					eMEDIUM eM = m_cCNT[pID].m_eNodalProperty[m - 1];
					double Cond = m_cCNT[pID].m_dConductanceList[m - 1];
					if (eM == BTE)
					{
						tripleList.push_back(T(Q, Q, -Cond));
						(*mB)(Q) += 0.0;
					}
					else if (eM == TPE)
					{
						tripleList.push_back(T(Q, Q, -Cond));
						(*mB)(Q) += -Cond * m_dAppliedVoltage;
					}
					else
					{
						tripleList.push_back(T(Q, P, +Cond));
						tripleList.push_back(T(Q, Q, -Cond));
						(*mB)(Q) += 0.0;
					}
				}
			}



		}
	}
	
	for (int i = 0; i < nNoofLinkedCond; i++)
	{
		int Q = slinkedcond[i].m_nGIndex1;
		int R = slinkedcond[i].m_nGIndex2;

		int eMed1 = slinkedcond[i].m_eMedium1;
		int eMed2 = slinkedcond[i].m_eMedium2;

		double Cond = slinkedcond[i].m_dConductance;

		if (eMed1 == MED)
		{
			if (eMed2 == BTE)
			{
				tripleList.push_back(T(Q, Q, -Cond));
				(*mB)(Q) += 0.0;
			}
			if (eMed2 == TPE)
			{
				tripleList.push_back(T(Q, Q, -Cond));
				(*mB)(Q) += -Cond * m_dAppliedVoltage;
			}
			if (eMed2 == MED)
			{
				tripleList.push_back(T(Q, Q, -Cond));
				tripleList.push_back(T(Q, R, +Cond));
				(*mB)(Q) += 0.0;
			}
		}

		if (eMed2 == MED)
		{
			if (eMed1 == BTE)
			{
				tripleList.push_back(T(R, R, -Cond));
				(*mB)(R) += 0.0;
			}
			if (eMed1 == TPE)
			{
				tripleList.push_back(T(R, R, -Cond));
				(*mB)(R) += -Cond * m_dAppliedVoltage;
			}
			if (eMed1 == MED)
			{
				tripleList.push_back(T(R, R, -Cond));
				tripleList.push_back(T(R, Q, +Cond));
				(*mB)(R) += 0.0;
			}

		}

	}

	// Matrix element setting
	mA->setFromTriplets(tripleList.begin(), tripleList.end());
	mA->makeCompressed();
}

double CPERC3D::ConductanceIndPercNetwork(int nNoofPercCNTs, int *plist, VectorXd mV, bool *error)
{
	double cond;

	double Iout = 0;
	double Iin = 0;

	for (int i = 0; i < nNoofPercCNTs; i++)
	{
		int pID = plist[i];
		int nNoofNode = m_cCNT[pID].m_nNoofNodes;

		for (int m = 0; m < nNoofNode - 1; m++)
		{
			if (m_cCNT[pID].m_eNodalProperty[m] == MED && m_cCNT[pID].m_eNodalProperty[m + 1] == TPE)
			{
				int giv = m_cCNT[pID].m_nGIndexList[m];
				double voltage = m_dAppliedVoltage - mV(giv);
				double conduct = m_cCNT[pID].m_dConductanceList[m];
				Iout += voltage * conduct;
			}

			if (m_cCNT[pID].m_eNodalProperty[m] == MED && m_cCNT[pID].m_eNodalProperty[m + 1] == BTE)
			{
				int giv = m_cCNT[pID].m_nGIndexList[m];
				double voltage = mV(giv);
				double conduct = m_cCNT[pID].m_dConductanceList[m];
				Iin += voltage * conduct;
			}
		}

		for (int m = 0; m < nNoofNode - 1; m++)
		{
			if (m_cCNT[pID].m_eNodalProperty[m] == TPE && m_cCNT[pID].m_eNodalProperty[m + 1] == MED)
			{
				int giv = m_cCNT[pID].m_nGIndexList[m + 1];
				double voltage = m_dAppliedVoltage - mV(giv);
				double conduct = m_cCNT[pID].m_dConductanceList[m];
				Iout += voltage * conduct;
			}
			if (m_cCNT[pID].m_eNodalProperty[m] == BTE && m_cCNT[pID].m_eNodalProperty[m + 1] == MED)
			{
				int giv = m_cCNT[pID].m_nGIndexList[m + 1];
				double voltage = mV(giv);
				double conduct = m_cCNT[pID].m_dConductanceList[m];
				Iin += voltage * conduct;
			}


		}

	}

	double accuracy = 100. - (100.*2.*fabs(Iout - Iin) / (Iout + Iin));

	printf("- Iout = %2.12e\n", Iout);
	printf("- Iinp = %2.12e\n", Iin);
	printf("- |Iout - Iin| = %2.12e\n", fabs(Iout - Iin));
	printf("- Accuracy: %2.3f[%%]\n", accuracy);

	double cond_out = Iout / m_dAppliedVoltage;
	double cond_in = Iin / m_dAppliedVoltage;

	cond = (cond_out + cond_in) / 2.;

	if (accuracy <= 90) *error = false;

	return cond;
}


void CPERC3D::SetVoltagesToCNT(int nNoofPercCNTs, int *plist, VectorXd mV)
{
	for (int i = 0; i < nNoofPercCNTs; i++)
	{
		int PID = plist[i];
		m_cCNT[PID].SetVoltage(m_dAppliedVoltage, &mV);
	}
}

void CPERC3D::SetVoltagesToLink(int nNoofLinkedCond, SLinkedCond *slinkedcond, VectorXd mV)
{
	for (int i = 0; i < nNoofLinkedCond; i++)
	{
		int gi1 = slinkedcond[i].m_nGIndex1;
		int gi2 = slinkedcond[i].m_nGIndex2;

		slinkedcond[i].m_dV1 = mV(gi1);
		slinkedcond[i].m_dV2 = mV(gi2);
	}
}

double CPERC3D::CalPowerDissipationCNT(int nNoofPercCNTs, int *plist)
{
	double ret = 0.0;
	for (int i = 0; i < nNoofPercCNTs; i++)
	{
		int PID = plist[i];
		ret += m_cCNT[PID].SetPowerDissipation();
	}

	return ret;
}

double CPERC3D::CalPowerDissipationLink(int nNoofLink, SLinkedCond *slinkedcond)
{
	double ret = 0.0;
	for (int i = 0; i < nNoofLink; i++)
	{
		double v = slinkedcond[i].m_dV1 - slinkedcond[i].m_dV2;
		slinkedcond[i].m_dP = v * v*slinkedcond[i].m_dConductance;
		ret += slinkedcond[i].m_dP;
	}

	return ret;
}

