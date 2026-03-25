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
