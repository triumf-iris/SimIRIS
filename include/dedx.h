#ifndef __DEDX_H
#define __DEDX_H

#include <fstream>
#include <string>
#include <array>
#include <utility>
#include <map>
#include "IrisMaterial.h"

class dedx{
 	public: 	
		
		double eFoil[100];	
		double eAl[100];	
		double eB[100];	
		double eC4H10[100];	
		double eCsI[100];
		double eTgt[100];	
		double eMy[100];	
		double eP[100];	
		double eSi[100];	
		double eSi3N4[100];	
		double eSiO2[100];	

		double dedxFoil[100];	
		double dedxAl[100];	
		double dedxB[100];	
		double dedxC4H10[100];	
		double dedxCsI[100];
		double dedxTgt[100];	
		double dedxMy[100];	
		double dedxP[100];	
		double dedxSi[100];	
		double dedxSi3N4[100];	
		double dedxSiO2[100];	

 	public:
  		dedx();//! Create
  		virtual ~dedx() {} //!
		
		void loadIncomingELoss(std::string, std::string, std::string, std::string, double);
		void loadOutgoingELoss(std::string, std::string, std::string, std::string, double);
		void loadELoss(std::string, IrisMaterial, double);
		std::array<double, 100> GetE(IrisMaterial material) { return tables[material].first; }
		std::array<double, 100> GetDeDx(IrisMaterial material) { return tables[material].second; }
  		void Clear();  //!
	protected:

 	private:
		std::map<IrisMaterial, std::pair<std::array<double, 100>, std::array<double, 100> > > tables;
  
};

#endif
