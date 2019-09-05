//	perc3d_structure.cpp													
//	Source code for the simulation of electrical percolation of	CNT/polymer composite
//	developed by Sungmin Jung, Ph.D.
//  University of Cambridge 
//  Contact: sm0104.jung@gmail.com, sj569@cam.ac.uk, s.jung@eng.cam.ac.uk
//  Copyright (C) 2019

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

// union and finding functions
int CPERC3D::Find(int i) // compressed with rank
{
	if (m_cCNT[i].m_nParent != i)
	{
		m_cCNT[i].m_nParent = Find(m_cCNT[i].m_nParent);
	}
	return m_cCNT[i].m_nParent;
}

void CPERC3D::Union(int i, int j) // compressed with rank
{
	int iroot = Find(i);
	int jroot = Find(j);

	if (iroot == jroot)
		return;

	if (m_cCNT[iroot].m_nRank < m_cCNT[jroot].m_nRank)
		m_cCNT[iroot].m_nParent = jroot;
	else if (m_cCNT[iroot].m_nRank > m_cCNT[jroot].m_nRank)
		m_cCNT[jroot].m_nParent = iroot;
	else
	{
		m_cCNT[jroot].m_nParent = iroot;
		m_cCNT[iroot].m_nRank = m_cCNT[iroot].m_nRank + 1;
	}
}

int CPERC3D::UnionTwoCNTs(FILE *fpbin, int i, int j)
{
	int pmin = 0, qmin = 0;
	double umin = 0.0, vmin = 0.0;
	double dmin = UnionCheck(i, j, &pmin, &qmin, &umin, &vmin);
	double ri = m_cCNT[i].m_dRadius + (m_cCNT[i].m_nNoofWall - 1)*dDistance_VanderWaals;
	double rj = m_cCNT[j].m_dRadius + (m_cCNT[j].m_nNoofWall - 1)*dDistance_VanderWaals;
	double D = ri + rj;

	int nlink = 0;

	if (CompareDistanceCutoff(i, j, dmin))
	{
		double zi = (1 - umin)*m_cCNT[i].m_vP[pmin](2) + umin * m_cCNT[i].m_vP[pmin+1](2);
		double zj = (1 - vmin)*m_cCNT[j].m_vP[qmin](2) + vmin * m_cCNT[j].m_vP[qmin+1](2);

		if (m_dZmin <= zi && zi <= m_dZmax && m_dZmin <= zj && zj <= m_dZmax )
		{
			SLinkedCond link;
			link.m_nCNTNo1 = i;
			link.m_nCNTNo2 = j;

			link.m_dU = m_cCNT[i].m_dL[pmin] + (m_cCNT[i].m_dL[pmin + 1] - m_cCNT[i].m_dL[pmin])*umin;
			link.m_dV = m_cCNT[j].m_dL[qmin] + (m_cCNT[j].m_dL[qmin + 1] - m_cCNT[j].m_dL[qmin])*vmin;

			link.m_dConductance = Gcontact(dmin, D);
			link.m_dDistance = dmin;
			link.m_nGIndex1 = 0;
			link.m_nGIndex2 = 0;

			link.m_eMedium1 = m_cCNT[i].ReturnNodalProperty(zi);
			link.m_eMedium2 = m_cCNT[j].ReturnNodalProperty(zj);

			if (link.m_eMedium1 == MED && link.m_eMedium2 == MED)
			{
				fwrite(&link, sizeof(link), 1, fpbin);
				Union(i, j);
				nlink = 1;
			}
					   			 		  
		}

//		printf("%d, %d\n", i, j);
	}

	return nlink;
}

double CPERC3D::UnionCheck(int i, int j, int *pmin, int *qmin, double *umin, double *vmin)
{
	int nNodeI = m_cCNT[i].m_nNoofNodes;
	int nNodeJ = m_cCNT[j].m_nNoofNodes;
	double ri = m_cCNT[i].m_dRadius + (m_cCNT[i].m_nNoofWall - 1)*dDistance_VanderWaals;
	double rj = m_cCNT[j].m_dRadius + (m_cCNT[j].m_nNoofWall - 1)*dDistance_VanderWaals;
	double D = ri + rj;

	double dmin = 100000000;

	for (int p = 0; p < nNodeI - 1; p++)
	{
		for (int q = 0; q < nNodeJ - 1; q++)
		{
			if (m_cCNT[i].m_eNodalProperty[p] == MED &&
				m_cCNT[j].m_eNodalProperty[q] == MED)
			{
				double u, v;

				double x1i = m_cCNT[i].m_vP[p](0);
				double y1i = m_cCNT[i].m_vP[p](1);
				double z1i = m_cCNT[i].m_vP[p](2);

				double x1f = m_cCNT[i].m_vP[p + 1](0);
				double y1f = m_cCNT[i].m_vP[p + 1](1);
				double z1f = m_cCNT[i].m_vP[p + 1](2);

				double x2i = m_cCNT[j].m_vP[q](0);
				double y2i = m_cCNT[j].m_vP[q](1);
				double z2i = m_cCNT[j].m_vP[q](2);

				double x2f = m_cCNT[j].m_vP[q+1](0);
				double y2f = m_cCNT[j].m_vP[q+1](1);
				double z2f = m_cCNT[j].m_vP[q+1](2);

				if (m_eBC == PBC)
				{
					if (x1i < m_dXmin || x1f < m_dXmin)
					{
						x1i = x1i + (m_dXmax - m_dXmin);
						x1f = x1f + (m_dXmax - m_dXmin);
					}
					if (x1i > m_dXmax || x1f > m_dXmax)
					{
						x1i = x1i - (m_dXmax - m_dXmin);
						x1f = x1f - (m_dXmax - m_dXmin);
					}

					if (y1i < m_dYmin || y1f < m_dYmin)
					{
						y1i = y1i + (m_dYmax - m_dYmin);
						y1f = y1f + (m_dYmax - m_dYmin);
					}
					if (y1i > m_dYmax || y1f > m_dYmax)
					{
						y1i = y1i - (m_dYmax - m_dYmin);
						y1f = y1f - (m_dYmax - m_dYmin);
					}
				}

				if (m_eBC == FBC)
				{

				}

				double dist = DistanceSegToSeg(x1i, y1i, z1i, x1f, y1f, z1f, x2i, y2i, z2i, x2f, y2f, z2f, &u, &v);

				if (dist < dmin)
				{
					dmin = dist;
					*pmin = p;
					*qmin = q;
					*umin = u;
					*vmin = v;
				}

			}

		}
	}

	return dmin;
}


int CPERC3D::UnionCNTs(void)
{
	int nCount = 0;
	FILE *fpbin = fopen("link_total.bin", "wb");

	int i;
	int a = 0;

	#pragma omp parallel for private(i) schedule(dynamic) reduction(+:nCount)
	for (i = 0; i < m_nNoofCNTs; i++)
	{
		for (int j = i + 1; j < m_nNoofCNTs; j++)
		{
			if ( (m_cCNT[i].m_dXcntmin <= m_cCNT[j].m_dXcntmax && m_cCNT[i].m_dXcntmax >= m_cCNT[j].m_dXcntmin) ||
				 (m_cCNT[i].m_dYcntmin <= m_cCNT[j].m_dYcntmax && m_cCNT[i].m_dYcntmax >= m_cCNT[j].m_dYcntmin) ||
				 (m_cCNT[i].m_dZcntmin <= m_cCNT[j].m_dZcntmax && m_cCNT[i].m_dZcntmax >= m_cCNT[j].m_dZcntmin) ) 
			{
				int nlink_unit = UnionTwoCNTs(fpbin, i, j);
				nCount += nlink_unit;
			}
		}
	}

	fclose(fpbin);

	FILE *fpbin_r = fopen("link_total.bin", "rb");
	SLinkedCond *link = new SLinkedCond[nCount];
	for (int l = 0; l < nCount; l++)
	{
		fread(&link[l], sizeof(link[l]), 1, fpbin_r);
	}
	fclose(fpbin_r);

	return nCount;
}

void CPERC3D::FindingGroup(int *glist)
{
	int i;

#pragma omp parallel for private(i) schedule(dynamic) //reduction(+:nCount)
	for (i = 0; i < m_nNoofCNTs; i++)
	{
		m_cCNT[i].m_nParent = Find(i);
	}

#pragma omp parallel for private(i) schedule(dynamic) //reduction(+:nCount)
	for (i = 0; i < m_nNoofCNTs; i++)
	{
		glist[i] = 0;

		for (int j = 0; j < m_nNoofCNTs; j++)
		{
			if (m_cCNT[j].m_nParent == i)
			{
				glist[i]++;
			}

		}
	}
}

void CPERC3D::ListingCNTsInGroup(int *glist, bool *b_perc)
{

	bool check_top = false;
	bool check_bot = false;

	int i;

#pragma omp parallel for private(i) schedule(dynamic) //reduction(+:nCount)
	for (i = 0; i < m_nNoofCNTs; i++)
	{
		b_perc[i] = false;

		check_top = false;
		check_bot = false;

		for (int j = 0; j < m_nNoofCNTs; j++)
		{
			if (m_cCNT[j].m_nParent == i)
			{
				if (m_cCNT[j].m_eCNTProperty == TPE) check_top = true;
				if (m_cCNT[j].m_eCNTProperty == BTE) check_bot = true;
			}
		}
		if (check_top && check_bot) b_perc[i] = true;
	}

}
