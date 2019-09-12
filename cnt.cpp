//	cnt.cpp													
//	Source code for the simulation of electrical percolation of	CNT/polymer composite
//	developed by Sungmin Jung, Ph.D. 
//  University of Cambridge
//  Last update on 18th April 2019
//  Contact: sm0104.jung@gmail.com, sj569@cam.ac.uk, s.jung@eng.cam.ac.uk



#include "cnt.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <stdlib.h>


using namespace std;

CCNT::CCNT()
{
	m_nParent = -1;
	m_nAddress = -1;
	m_nRank = 0;

	m_eCNTProperty = MED;

	m_vP = 0;

	m_dL = 0;
	m_dV = 0;
	m_dP = 0;
//	m_dR = 0;

	m_eNodalProperty = 0;
	m_nGIndexList = 0;
	m_dConductanceList = 0;
}

CCNT::~CCNT()
{
	if (m_vP != 0) delete[] m_vP;

	if (m_dL != 0) delete[] m_dL;
	if (m_dV != 0) delete[] m_dV;

	if (m_eNodalProperty != 0) delete[] m_eNodalProperty;
	if (m_nGIndexList != 0) delete[] m_nGIndexList;
	if (m_dConductanceList != 0) delete[] m_dConductanceList;
	if (m_dP != 0) delete[] m_dP;
//	if (m_dR != 0) delete[] m_dR;
}

void CCNT::SetNoofWall(int nNoofWall)
{
	m_nNoofWall = nNoofWall;
}

double CCNT::WeightCNT(void)
{
	double dWeightCNT = 0.0;

	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		for (int n = 0; n < m_nNoofWall; n++)
		{
			double radius = m_dRadius - dDistance_VanderWaals * (double) n;
			dWeightCNT += dSurfaceDensityCNT*2.*M_PI*radius*(m_dL[p + 1] - m_dL[p]);
		}
	}

	return dWeightCNT;
}

double CCNT::VolumeCNT(void)
{
	double dVolumeCNT = 0.0;
	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		double r_W = m_dRadius;
		dVolumeCNT += M_PI * r_W*r_W*(m_dL[p + 1] - m_dL[p]);
	}

	return dVolumeCNT;
}

CCNT& CCNT::operator = (const CCNT& x)
{
	m_nParent = x.m_nParent;
	m_nAddress = x.m_nAddress;
	m_nRank = x.m_nRank;

	m_eCNTProperty = x.m_eCNTProperty;

	m_dXmax = x.m_dXmax;
	m_dXmin = x.m_dXmin;

	m_dYmax = x.m_dYmax;
	m_dYmin = x.m_dYmin;

	m_dZmax = x.m_dZmax;
	m_dZmin = x.m_dZmin;

	m_dXcntmax = x.m_dXcntmax;
	m_dXcntmin = x.m_dXcntmin;

	m_dYcntmax = x.m_dYcntmax;
	m_dYcntmin = x.m_dYcntmin;

	m_dZcntmax = x.m_dZcntmax;
	m_dZcntmin = x.m_dZcntmin;


	m_dLength = x.m_dLength;
	m_dRadius = x.m_dRadius;
	m_nNoofWall = x.m_nNoofWall;

	m_dConductivity = x.m_dConductivity;

	m_nNoofNodes = x.m_nNoofNodes;

	m_vP = new Vector3d[m_nNoofNodes];

	m_dL = new double[m_nNoofNodes];

	m_dV = new double[m_nNoofNodes];
	
	m_eNodalProperty = new eMEDIUM[m_nNoofNodes];
	m_nGIndexList = new int[m_nNoofNodes];
	m_dConductanceList = new double[m_nNoofNodes - 1];
	m_dP = new double[m_nNoofNodes - 1];
//	m_dR = new double[m_nNoofNodes - 1];

	memcpy(m_vP, x.m_vP, m_nNoofNodes * sizeof(Vector3d));
	memcpy(m_dL, x.m_dL, m_nNoofNodes * sizeof(double));
	memcpy(m_dV, x.m_dV, m_nNoofNodes * sizeof(double));

	memcpy(m_eNodalProperty, x.m_eNodalProperty, m_nNoofNodes * sizeof(eMEDIUM));
	memcpy(m_nGIndexList, x.m_nGIndexList, m_nNoofNodes * sizeof(int));
	memcpy(m_dConductanceList, x.m_dConductanceList, (m_nNoofNodes - 1) * sizeof(double));
	memcpy(m_dP, x.m_dP, (m_nNoofNodes - 1) * sizeof(double));
//	memcpy(m_dR, x.m_dR, (m_nNoofNodes - 1) * sizeof(double));

	return *this;
}

void CCNT::MaksSet(int index)
{
	m_nParent = index;
	m_nAddress = index;
}

void CCNT::SetDomain(double dXmin, double dXmax, double dYmin, double dYmax, double dZmin, double dZmax)
{
	m_dXmin = dXmin;
	m_dYmin = dYmin;
	m_dZmin = dZmin;

	m_dXmax = dXmax;
	m_dYmax = dYmax;
	m_dZmax = dZmax;
}
/*
void CCNT::SetConductance(void)
{
	m_dConductivity = dG0 * m_dLength / (M_PI*m_dRadius*m_dRadius);

	if (m_dConductanceList != 0) delete[] m_dConductanceList;
	m_dConductanceList = new double[m_nNoofNodes - 1];

	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		Vector3d vP = m_vP[p + 1] - m_vP[p];

		double dl = vP.norm();

		double dPartialConductance = m_dConductivity * (M_PI*m_dRadius* m_dRadius) / dl ;
		m_dConductanceList[p] = dPartialConductance;
		if (m_eNodalProperty[p] == TPE && m_eNodalProperty[p + 1] == TPE) m_dConductanceList[p] = 10000000;
		if (m_eNodalProperty[p] == BTE && m_eNodalProperty[p + 1] == BTE) m_dConductanceList[p] = 10000000;
	}


}
*/
void CCNT::SetConductance(double dConductivity) 
{
	m_dConductivity = dConductivity;

	if (m_dConductanceList != 0) delete[] m_dConductanceList;
	m_dConductanceList = new double[m_nNoofNodes - 1];

	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		Vector3d vP = m_vP[p + 1] - m_vP[p];

		double dl = vP.norm();

		double dPartialConductance = m_dConductivity * M_PI*m_dRadius*m_dRadius / dl;

		m_dConductanceList[p] = dPartialConductance;
		if (m_eNodalProperty[p] == TPE && m_eNodalProperty[p + 1] == TPE) m_dConductanceList[p] = 10000000;
		if (m_eNodalProperty[p] == BTE && m_eNodalProperty[p + 1] == BTE) m_dConductanceList[p] = 10000000;
	}
}

void CCNT::SetVoltage(double dAppliedVoltage, VectorXd *mV)
{
	if (m_dV != 0) delete[] m_dV;

	m_dV = new double[m_nNoofNodes];

	for (int p = 0; p < m_nNoofNodes; p++)
	{
		int gi = m_nGIndexList[p];
		if (gi == -1) m_dV[p] = dAppliedVoltage;
		else if (gi == -2) m_dV[p] = 0.0;
		else
		{
			m_dV[p] = (*mV)(gi);
		}
	}
}

double CCNT::SetPowerDissipation(void)
{
	double ret = 0.;

	if (m_dP != 0) delete[] m_dP;

	m_dP = new double[m_nNoofNodes - 1];

	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		double v = (m_dV[p + 1] - m_dV[p]);
		m_dP[p] = v * v*m_dConductanceList[p];
		ret += m_dP[p];
	}
	return ret;
}

void CCNT::PrintNodalInfo(FILE *fp)
{
	fprintf(fp, "Address: %d\n", m_nAddress);
	fprintf(fp, "Parent:  %d\n", m_nParent);
	fprintf(fp, "Length: %f, Radius: %f\n", m_dLength, m_dRadius);
	fprintf(fp, "No. of Wall: %d\n", m_nNoofWall);
	fprintf(fp, "Total No. of Nodes: %d\n", m_nNoofNodes);
	fprintf(fp, "Node: X\t\tY\tZ(MEDIUM) \tL\tIndex\tVTG\tCND\t\tPower\n");

	for (int i = 0; i < m_nNoofNodes; i++)
	{
		char cmed[255];
		if (m_eNodalProperty[i] == MED) sprintf(cmed, "MED");
		if (m_eNodalProperty[i] == TPE) sprintf(cmed, "TPE");
		if (m_eNodalProperty[i] == BTE) sprintf(cmed, "BTE");

		if (i < m_nNoofNodes - 1)
		{
			fprintf(fp, "%4d: %2.3f\t%2.3f\t%2.3f(%s)\t%2.3f\t%d\t%2.3f\t%2.3e\t%2.3e\n", i, m_vP[i](0), m_vP[i](1), m_vP[i](2), cmed, m_dL[i], m_nGIndexList[i], m_dV[i], m_dConductanceList[i], m_dP[i]);
		}
		else
		{
			fprintf(fp, "%4d: %2.3f\t%2.3f\t%2.3f(%s)\t%2.3f\t%d\t%2.3f\n", i, m_vP[i](0), m_vP[i](1), m_vP[i](2), cmed,  m_dL[i], m_nGIndexList[i], m_dV[i]);
		}
	}

	fprintf(fp, "\n");

}

eMEDIUM CCNT::ReturnNodalProperty(double dZ)
{
	eMEDIUM ret;
	if (dZ <= m_dZmin) ret = BTE;
	else if (dZ >= m_dZmax) ret = TPE;
	else ret = MED;

	return ret;
}

int CCNT::InsertNode(double L, int *nGI, eMEDIUM *eMed)
{
	int p;
	double u;

	if (L < m_dL[0])
	{
		p = 0;
		u = 0.0;
	}
	else if (L >= m_dL[m_nNoofNodes - 1])
	{
		p = m_nNoofNodes - 1;
		u = 0.0;
	}
	else if (m_dL[0] <= L && L < m_dL[m_nNoofNodes - 1])
	{
		for (p = 0; p < m_nNoofNodes - 1; p++)
		{
			if (m_dL[p] <= L && L < m_dL[p + 1])
			{
				break;
			}
		}
		u = (L - m_dL[p]) / (m_dL[p + 1] - m_dL[p]);
		if (fabs(u) < 1.0e-9)
		{
			u = 0.0;
		}
		if (fabs(1.0 - u) < 1.0e-9)
		{
			p = p + 1;
			u = 0.0;
		}
	}
	else
	{
		p = m_nNoofNodes - 1;
		u = 0.0;

	}

	int ret = InsertNode(p, u, nGI, eMed);

	return ret;
}

int CCNT::InsertNode(int index, double u, int *nGI, eMEDIUM *eMed)
{
	int nRetGlobalNode = 0;

	if (0 <= index && index < m_nNoofNodes - 1)
	{
		if (u <= 0)
		{
			nRetGlobalNode = m_nGIndexList[index];
			*eMed = m_eNodalProperty[index];
		}
		else if (u >= 1)
		{
			nRetGlobalNode = m_nGIndexList[index + 1];
			*eMed = m_eNodalProperty[index + 1];
		}
		else
		{
			Vector3d *vP = new Vector3d[m_nNoofNodes + 1];

			eMEDIUM *eNodalProperty = new eMEDIUM[m_nNoofNodes + 1];
			int *nGIndexList = new int[m_nNoofNodes + 1];

			for (int i = 0; i <= index; i++)
			{
				vP[i] = m_vP[i];
				eNodalProperty[i] = m_eNodalProperty[i];
				nGIndexList[i] = m_nGIndexList[i];
			}

			vP[index + 1] = u * (m_vP[index + 1] - m_vP[index]) + m_vP[index];

			eNodalProperty[index + 1] = ReturnNodalProperty(vP[index + 1](2));
			*eMed = eNodalProperty[index + 1];
			nGIndexList[index + 1] = *nGI;
			nRetGlobalNode = *nGI;
			(*nGI)++;

			for (int i = index + 2; i <= m_nNoofNodes; i++)
			{
				vP[i] = m_vP[i - 1];
				eNodalProperty[i] = m_eNodalProperty[i - 1];
				nGIndexList[i] = m_nGIndexList[i - 1];
			}

			delete[] m_vP;
			delete[] m_eNodalProperty;
			delete[] m_nGIndexList;

			m_nNoofNodes++;

			m_vP = new Vector3d[m_nNoofNodes];
			m_eNodalProperty = new eMEDIUM[m_nNoofNodes];
			m_nGIndexList = new int[m_nNoofNodes];
			for (int i = 0; i < m_nNoofNodes; i++)
			{
				m_vP[i] = vP[i];
				m_eNodalProperty[i] = eNodalProperty[i];
				m_nGIndexList[i] = nGIndexList[i];
			}

			// voltage resetting
			if (m_dV != 0) delete[] m_dV;
			m_dV = new double[m_nNoofNodes];
			for (int p = 0; p < m_nNoofNodes; p++) m_dV[p] = 0.0;

			// conductance recalculation
			SetConductance(m_dConductivity);

			// power dissipation resetting
			if (m_dP != 0) delete[] m_dP;
			m_dP = new double[m_nNoofNodes - 1];
			for (int p = 0; p < m_nNoofNodes-1; p++) m_dP[p] = 0.0;

			// length paramerter recalculation
			if (m_dL != 0) delete[] m_dL;
			m_dL = new double[m_nNoofNodes];

			m_dL[0] = 0;
			for (int p = 1; p < m_nNoofNodes; p++)
			{
				Vector3d dP = m_vP[p] - m_vP[p - 1];
				double l = dP.norm();
				m_dL[p] = m_dL[p - 1] + l;
			}

			delete[] vP;
			delete[] eNodalProperty;
			delete[] nGIndexList;
		}

	}
	else
	{
		nRetGlobalNode = m_nGIndexList[m_nNoofNodes - 1];
		*eMed = m_eNodalProperty[m_nNoofNodes - 1];
	}

	return nRetGlobalNode;

}

void CCNT::SetCNTBoundary(void)
{
	// Searching Electrode Boundary and inserting nodes
	int nz = 0;
	double *dLz = new double[m_nNoofNodes];

	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		if (    ( m_vP[p](2) < m_dZmax && m_dZmax < m_vP[p+1](2) ) || ( m_vP[p+1](2) < m_dZmax && m_dZmax < m_vP[p](2) )     )
		{
			double uz = (m_dZmax - m_vP[p](2)) / (m_vP[p+1](2) - m_vP[p](2));
			dLz[nz] = m_dL[p] + uz*(m_dL[p + 1] - m_dL[p]);
			nz++;
		}

		if (    ( m_vP[p](2) < m_dZmin && m_dZmin < m_vP[p+1](2) ) || ( m_vP[p+1](2) < m_dZmin && m_dZmin < m_vP[p](2) )      )  
		{
			double uz = (m_dZmin - m_vP[p](2)) / (m_vP[p+1](2) - m_vP[p](2));
			dLz[nz] = m_dL[p] + uz * (m_dL[p + 1] - m_dL[p]);
			nz++;
		}

	}

	// Insertion of Electrode Boundary for Z-direction
	for (int l = 0; l < nz; l++)
	{
		int gn_z = 0;
		eMEDIUM eMed_z;
		InsertNode(dLz[l], &gn_z, &eMed_z);
	}

	delete[] dLz;

	// x-direction...
	int nx = 0;
	double *dLx = new double[m_nNoofNodes];

	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		if (    ( m_vP[p](0) < m_dXmax && m_dXmax < m_vP[p+1](0) ) || ( m_vP[p+1](0) < m_dXmax && m_dXmax < m_vP[p](0) )     )
		{
			double ux = (m_dXmax - m_vP[p](0)) / (m_vP[p+1](0) - m_vP[p](0));
			dLx[nx] = m_dL[p] + ux * (m_dL[p + 1] - m_dL[p]);
			nx++;
		}
		if (    ( m_vP[p](0) < m_dXmin && m_dXmin < m_vP[p+1](0) ) || ( m_vP[p+1](0) < m_dXmin && m_dXmin < m_vP[p](0) )     )
		{
			double ux = (m_dXmin - m_vP[p](0)) / (m_vP[p+1](0) - m_vP[p](0));
			dLx[nx] = m_dL[p] + ux * (m_dL[p + 1] - m_dL[p]);
			nx++;
		}

	}

	for (int l = 0; l < nx; l++)
	{
		int gn_x = 0;
		eMEDIUM eMed_x;
		InsertNode(dLx[l], &gn_x, &eMed_x);
	}

	delete[] dLx;

	// y-direction...
	int ny = 0;
	double *dLy = new double[m_nNoofNodes];

	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		if (    ( m_vP[p](1) < m_dYmax && m_dYmax < m_vP[p+1](1) ) || ( m_vP[p+1](1) < m_dYmax && m_dYmax < m_vP[p](1) )     )
		{
			double uy = (m_dYmax - m_vP[p](1)) / (m_vP[p+1](1) - m_vP[p](1));
			dLy[ny] = m_dL[p] + uy * (m_dL[p + 1] - m_dL[p]);
			ny++;
		}

		if ((m_vP[p](1) < m_dYmin && m_dYmin < m_vP[p + 1](1)) ||
			(m_vP[p + 1](1) < m_dYmin && m_dYmin < m_vP[p](1)))
		{
			double uy = (m_dYmin - m_vP[p](1)) / (m_vP[p + 1](1) - m_vP[p](1));
			dLy[ny] = m_dL[p] + uy * (m_dL[p + 1] - m_dL[p]);
			ny++;
		}
	}

	for (int l = 0; l < ny; l++)
	{
		int gn_y = 0;
		eMEDIUM eMed_y;
		InsertNode(dLy[l], &gn_y, &eMed_y);
	}
	delete[] dLy;

	// CNT property 
	for (int l = 0; l < m_nNoofNodes; l++)
	{
		if (m_eNodalProperty[l] == TPE)
		{
			m_eCNTProperty = TPE;
		}

		if (m_eNodalProperty[l] == BTE)
		{
			m_eCNTProperty = BTE;
		}
	}

}

void CCNT::SetNodalProperties(void)
{
	for (int i = 0; i < m_nNoofNodes; i++)
	{
		m_eNodalProperty[i] = ReturnNodalProperty(m_vP[i](2));
	}
}

void CCNT::SetInitLineSegments(eGEOM eGeom, int nNoofNodes, double dLength, double dRadius, double dtheta_dev_dgr)
{
	m_nNoofNodes = nNoofNodes;

	m_dLength = dLength;

	m_dRadius = dRadius;

	m_vP = new Vector3d[m_nNoofNodes];
	m_dL = new double[m_nNoofNodes];

	m_eNodalProperty = new eMEDIUM[m_nNoofNodes];
	m_nGIndexList = new int[m_nNoofNodes];

	for (int i = 0; i < m_nNoofNodes; i++)
	{
		m_eNodalProperty[i] = MED;
		m_nGIndexList[i] = 0;
	}

	double dl = m_dLength / (double)(m_nNoofNodes - 1);

	bool check = true;

	do {

		check = true;

		m_vP[0](0) = ((double)rand() / (double)RAND_MAX)*(m_dXmax - m_dXmin) + m_dXmin;
		m_vP[0](1) = ((double)rand() / (double)RAND_MAX)*(m_dYmax - m_dYmin) + m_dYmin;
		m_vP[0](2) = ((double)rand() / (double)RAND_MAX)*(m_dZmax - m_dZmin) + m_dZmin;

		m_eNodalProperty[0] = ReturnNodalProperty(m_vP[0](2));

		double theta0 = 90.*(1. - 2.* (double)rand() / (double)RAND_MAX);
		double phi0 = 360.*(double)rand() / (double)RAND_MAX;


		Vector3d n0 = GetUnitVector(theta0, phi0);
		
		m_vP[1] = m_vP[0] + dl * n0;

		m_eNodalProperty[1] = ReturnNodalProperty(m_vP[1](2));

		for (int i = 2; i < m_nNoofNodes; i++)
		{
			Vector3d vp = m_vP[i - 1] - m_vP[i - 2];
			vp.normalize();
			double tht1, phi1;
			GetAngles(vp, &tht1, &phi1);

			double dtht = dtheta_dev_dgr * (1. - 2.*(double)rand() / (double)RAND_MAX);
			double phi2 = dtor(360.*((double)rand() / (double)RAND_MAX));

			Vector3d vq = GetUnitVector(tht1+dtht, phi1);

			Matrix3d mR;
			double u = vp(0), v = vp(1), w = vp(2);

			mR(0,0) = u * u + (v*v + w * w)*cos(phi2);
			mR(0,1) = u * v*(1. - cos(phi2)) - w * sin(phi2);
			mR(0,2) = u * w*(1. - cos(phi2)) + v * sin(phi2);
			mR(1,0) = u * v*(1. - cos(phi2)) + w * sin(phi2);
			mR(1,1) = v * v + (u*u + w * w)*cos(phi2);
			mR(1,2) = v * w*(1. - cos(phi2)) - u * sin(phi2);
			mR(2,0) = u * w*(1. - cos(phi2)) - v * sin(phi2);
			mR(2,1) = v * w*(1. - cos(phi2)) + u * sin(phi2);
			mR(2,2) = w * w + (u*u + v * v)*cos(phi2);
			
			Vector3d n = mR * vq;

			n.normalize();

			if (n.dot(vp) < 0.0)
			{
				n = -n;
			}

			m_vP[i] = m_vP[i - 1] + dl * n;

			m_eNodalProperty[i] = ReturnNodalProperty(m_vP[i](2));
		}

		if (eGeom == CYL)
		{
			double xc = 0.5*(m_dXmax + m_dXmin);
			double yc = 0.5*(m_dYmax + m_dYmin);
			double rx = (m_dXmax - m_dXmin) / 2.;
			double ry = (m_dYmax - m_dYmin) / 2.;

			for (int i = 0; i < m_nNoofNodes; i++)
			{
				double f = (m_vP[i](0) - xc) * (m_vP[i](0) - xc) / (rx*rx) + (m_vP[i](1) - yc) * (m_vP[i](1) - yc) / (ry*ry);

				if (f > 1.) check = false;
			}

		}
	} while (!check);

	// length paramerter recalculation
	m_dL = new double[m_nNoofNodes];

	m_dL[0] = 0;

	for (int p = 1; p < m_nNoofNodes; p++)
	{
		Vector3d dp = m_vP[p] - m_vP[p - 1];
		
		double l = dp.norm();

		m_dL[p] = m_dL[p - 1] + l;
	}

	m_dV = new double[m_nNoofNodes];

	for (int p = 0; p < m_nNoofNodes; p++)
	{
		m_dV[p] = 0.0;
	}

	m_dConductanceList = new double[m_nNoofNodes - 1];
	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		m_dConductanceList[p] = 0.0;
	}

	m_dP = new double[m_nNoofNodes - 1];
	for (int p = 0; p < m_nNoofNodes - 1; p++)
	{
		m_dP[p] = 0.0;
	}

	if (m_nNoofNodes % 2 == 0)
	{
		int a = m_nNoofNodes / 2 - 1;
		int b = a + 1;

		m_vPc = 0.5*(m_vP[a] + m_vP[b]);
	}
	else
	{
		int a = (m_nNoofNodes-1) / 2;

		m_vPc = m_vP[a];
	}

}

void CCNT::Stretching(double dPoissonPolymer, double dPoissonCNT, double dStrainPolymerLong)
{
	double dConstLong = (1. + dStrainPolymerLong);
	double dStrainPolymerTran = - dPoissonPolymer * dStrainPolymerLong;
	double dConstTran = (1. + dStrainPolymerTran);

	m_dXmax = dConstTran * m_dXmax;
	m_dYmax = dConstTran * m_dYmax;
	m_dZmax = dConstLong * m_dZmax;

	m_dXmin = dConstTran * m_dXmin;
	m_dYmin = dConstTran * m_dYmin;
	m_dZmin = dConstLong * m_dZmin;

	for (int p = 0; p < m_nNoofNodes; p++)
	{
		m_vP[p](0) = dConstTran * m_vP[p](0);
		m_vP[p](1) = dConstTran * m_vP[p](1);
		m_vP[p](2) = dConstLong * m_vP[p](2);
	}

	double *dLnew = new double[m_nNoofNodes];

	dLnew[0] = 0;

	for (int p = 1; p < m_nNoofNodes; p++)
	{
		Vector3d dp = m_vP[p] - m_vP[p-1];
		double length = dp.norm();
	
		dLnew[p] = dLnew[p-1] + length;		
	}

	double lnew = dLnew[m_nNoofNodes - 1];
	double lold = m_dL[m_nNoofNodes - 1];
	double epsi_long = (lnew - lold) / lold;
	double epsi_tran = -dPoissonCNT * epsi_long;

	m_dRadius = (1. + epsi_tran)*m_dRadius;


	for (int p = 0; p < m_nNoofNodes; p++)
	{
		m_dL[p] = dLnew[p];
	}

	m_dLength = m_dL[m_nNoofNodes - 1];

	delete[] dLnew;

	SetConductance(m_dConductivity);
}


void CCNT::Stretching2(double dPoissonPolymer, double dStrainPolymerLong)
{

	double dConstLong = (1. + dStrainPolymerLong);
//	double dConstTran = (1. - dPoissonPolymer * dStrainPolymerLong);
	double dConstTran = pow(dConstLong, -dPoissonPolymer);

	m_dXmax = dConstTran * m_dXmax;
	m_dYmax = dConstTran * m_dYmax;		
	m_dZmax = dConstLong * m_dZmax;

	m_dXmin = dConstTran * m_dXmin;
	m_dYmin = dConstTran * m_dYmin;
	m_dZmin = dConstLong * m_dZmin;	

	Vector3d n = Vector3d(0., 0., 0.);
	Vector3d vPc0 = m_vPc;
	Vector3d vPc1;
	double dTheta0, dPhi0;
	double dTheta1, dPhi1;

	for (int p = 1; p < m_nNoofNodes; p++)
	{
		n += (m_vP[p] - m_vP[p - 1]);
	}
	n = n / m_dLength;

	GetAngles(n, &dTheta0, &dPhi0);

	dPhi1 = dPhi0;
//	dTheta1 = rtod(    atan(tan(dtor(dTheta0))*(1 + dStrainPolymerLong) / (1. + dStrainPolymerTran))    );
	dTheta1 = rtod(    atan(tan(dtor(dTheta0))*dConstLong/ dConstTran));

	vPc1(0) = dConstTran*vPc0(0);
	vPc1(1) = dConstTran*vPc0(1);
	vPc1(2) = dConstLong*vPc0(2);

	double dDeltaTheta = dTheta1 - dTheta0;

	double u = sin(dtor(dPhi0));
	double v = -cos(dtor(dPhi0));
	double cdt = cos(dtor(dDeltaTheta));
	double sdt = sin(dtor(dDeltaTheta));

	Matrix3d mR;

	mR(0, 0) = cdt + u * u*(1. - cdt);
	mR(0, 1) = u * v*(1. - cdt);
	mR(0, 2) = v * sdt;
	mR(1, 0) = u * v*(1. - cdt);
	mR(1, 1) = cdt + v * v*(1 - cdt);
	mR(1, 2) = -u * sdt;
	mR(2, 0) = -v * sdt;
	mR(2, 1) = u * sdt;
	mR(2, 2) = cdt;

	Vector3d *vPnew = new Vector3d[m_nNoofNodes];
	
	for (int p = 0; p < m_nNoofNodes; p++)
	{
		Vector3d vp = m_vP[p] - vPc0;

		vPnew[p] = mR * vp + vPc1;
	}


	for (int p = 0; p < m_nNoofNodes; p++)
	{
		m_vP[p] = vPnew[p];
	}

	m_vPc = vPc1;

	m_dL = new double[m_nNoofNodes];

	m_dL[0] = 0;

	for (int p = 1; p < m_nNoofNodes; p++)
	{
		Vector3d dp = m_vP[p] - m_vP[p - 1];

		double l = dp.norm();

		m_dL[p] = m_dL[p - 1] + l;
	}


	SetConductance(m_dConductivity);

	delete[] vPnew;
}

void CCNT::GetAngles(Vector3d n, double *theta_dgr, double *phi_dgr)
{
	double nx = n(0);
	double ny = n(1);
	double nz = n(2);

	if ( -1.0 <= nz && nz <= 1.0)
	{
		*theta_dgr = rtod(asin(nz));
	}
	else if (nz > 1.0)
	{
		*theta_dgr = 90.0;
	}
	else 
	{
		*theta_dgr = -90.0;
	}

	if (nx >= 0.0 && ny >= 0.0)
	{
		if (nx > 0)
		{
			*phi_dgr = rtod(atan(ny / nx));
		}
		else
		{
			*phi_dgr = 90.0;
		}
	}

	if (nx >= 0.0 && ny < 0.0)
	{
		if (nx > 0.0)
		{
			*phi_dgr = 360. + rtod(atan(ny / nx));
		}
		else
		{
			*phi_dgr = 270.;
		}
	}

	if (nx < 0.0 && ny >= 0.0)
	{
		*phi_dgr = 180. + rtod(atan(ny / nx));
	}

	if (nx < 0.0 && ny < 0.0)
	{
		*phi_dgr = 180. + rtod(atan(ny / nx));
	}

}

Vector3d CCNT::GetUnitVector(double theta_dgr, double phi_dgr)
{
	double theta_rad = dtor(theta_dgr);
	double phi_rad = dtor(phi_dgr);

	Vector3d n;

	n(0) = cos(theta_rad)*cos(phi_rad);
	n(1) = cos(theta_rad)*sin(phi_rad);
	n(2) = sin(theta_rad);

	return n;
}

void CCNT::AddAngles(double *theta_dgr, double *phi_dgr, double dtheta_dgr, double dphi_dgr)
{
	double tht = *theta_dgr + dtheta_dgr;
	double phi = *phi_dgr + dphi_dgr;

	if (tht < -90.)
	{
		tht = tht + 180.;
		phi = phi + 180.;
	}
	else if (tht > 90.)
	{
		tht = tht - 180.;
		phi = phi + 180.;
	}

	if (phi < 0)
	{
		phi = phi + 360.;
	}

	else if (phi > 360)
	{
		phi = phi - 360.;
	}

	*theta_dgr = tht;
	*phi_dgr = phi;
}

void CCNT::SetCNTRange(void)
{
	m_dXcntmax = -100000; 	m_dXcntmin = +100000;
	m_dYcntmax = -100000; 	m_dYcntmin = +100000;
	m_dZcntmax = -100000;	m_dZcntmin = +100000;

	for (int p = 0; p < m_nNoofNodes; p++)
	{
		if (m_vP[p](0) >= m_dXcntmax) m_dXcntmax = m_vP[p](0);
		if (m_vP[p](0) <= m_dXcntmin) m_dXcntmin = m_vP[p](0);
		if (m_vP[p](1) >= m_dYcntmax) m_dYcntmax = m_vP[p](1);
		if (m_vP[p](1) <= m_dYcntmin) m_dYcntmin = m_vP[p](1);
		if (m_vP[p](2) >= m_dZcntmax) m_dZcntmax = m_vP[p](2);
		if (m_vP[p](2) <= m_dZcntmin) m_dZcntmin = m_vP[p](2);
	}

	m_dXcntmax += m_dRadius;
	m_dXcntmin -= m_dRadius;

	m_dYcntmax += m_dRadius;
	m_dYcntmin -= m_dRadius;

	m_dZcntmax += m_dRadius;
	m_dZcntmin -= m_dRadius;
}
