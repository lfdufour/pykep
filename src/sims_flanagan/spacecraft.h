/*****************************************************************************
 *   Copyright (C) 2004-2015 The PyKEP development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef KEP_TOOLBOX_SPACECRAFT_H
#define KEP_TOOLBOX_SPACECRAFT_H

#include <iostream>
#include <algorithm>    // std::max

// Serialization code
#include "../serialization.h"
// Serialization code (END)
#include "../config.h"

namespace kep_toolbox {
namespace sims_flanagan{

/// Spacecraft
/**
 * A container for system design parameters of a spacecraft.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */


class __KEP_TOOL_VISIBLE spacecraft
{
	friend std::ostream &operator<<(std::ostream &s, const spacecraft &in );
public:
	spacecraft():m_mass(0),m_thrust(0),m_isp(0), m_AT(0), m_BT(0), m_maxP(0), m_minP(0), m_Pmargin(0) {}
	spacecraft(const double &mass_, const double &thrust_, const double &isp_) 
		: m_mass(mass_), m_thrust(thrust_), m_isp(isp_), m_AT(0), m_BT(0), m_maxP(0), m_minP(0), m_Pmargin(0) {}
	spacecraft(const double &mass_, const double &thrust_, const double &isp_,
		  const double &AT_, const double &BT_, const double &maxP_,const double &minP_,
		  const double &P1AU_, const double &Pmargin_) 
		: m_mass(mass_),m_thrust(thrust_),m_isp(isp_), m_AT(AT_), m_BT(BT_), m_maxP(maxP_), m_minP(minP_), m_Pmargin(Pmargin_) {}
	double get_mass() const {return m_mass;}
	double get_thrust() const {return m_thrust;}
	double get_thrust_electricSolar(double distanceSun) const {
		double P = std::max(m_P1AU/(distanceSun-1.0)/(distanceSun - 1.0)-m_Pmargin,0.0);
		double thrustElectric = 0.0;
		if(P<=m_minP){
			thrustElectric = 0.0;
			return thrustElectric;
		}else if(P>m_maxP){
			P = m_Pmax;
		}else{
			P = P;
		}
		thrustElectric  = m_AT + m_BT*P;
		return thrustElectric ; 
	}
	double get_isp() const {return m_isp;}
	void set_mass(const double _mass) {m_mass=_mass;}
	void set_thrust(const double _thrust) {m_thrust=_thrust;}
	void set_isp(const double _isp) {m_isp=_isp;}
	std::string human_readable() const;
private:
// Serialization code
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & m_mass;
		ar & m_thrust;
		ar & m_isp;
		ar &  m_AT;
		ar &  m_BT;
		ar &  m_maxP;
		ar &  m_minP;
		ar &  m_P1AU;
		ar &  m_Pmargin
	}
// Serialization code (END)
	double m_mass;
	double m_thrust;
	double m_isp;
	double m_AT;
	double m_BT;
	double m_maxP;
	double m_minP;
	double m_P1AU;
	double m_Pmargin;
};

__KEP_TOOL_VISIBLE std::ostream &operator<<(std::ostream &s, const spacecraft &in );


}} //Namespaces

#endif // KEP_TOOLBOX_SPACECRAFT_H
