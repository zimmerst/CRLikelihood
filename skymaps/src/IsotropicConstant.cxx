/** @file IsotropicConstant.cxx
    @brief implement IsotropicConstant
$Header: /nfs/slac/g/glast/ground/cvs/skymaps/src/IsotropicConstant.cxx,v 1.2 2011/06/22 02:47:46 lande Exp $
*/

#include "skymaps/IsotropicConstant.h"
using namespace skymaps;

IsotropicConstant::IsotropicConstant(double constant)
:  m_constant(constant)
{
    setName("constant");
}

IsotropicConstant::~IsotropicConstant()
{
}

double IsotropicConstant::value(const astro::SkyDir& /*dir*/, double energy)const
{
    return m_constant;
}

double IsotropicConstant::integral(const astro::SkyDir& /*dir*/, double a, double b)const
{
    return m_constant*(b-a);
}
