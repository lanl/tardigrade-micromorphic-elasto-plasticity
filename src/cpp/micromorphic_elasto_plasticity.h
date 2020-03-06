/*!
 * micromorphic_elasto_plasticity.h
 *
 * An implementation of a elasto-plastic micromorphic constitutive model 
 * following the derivations of Farhad Shahabi in his dissertation.
 */

#ifndef MICROMORPHIC_ELASTO_PLASTICITY_H
#define MICROMORPHIC_ELASTO_PLASTICITY_H

#include<error_tools.h>
#define USE_EIGEN
#include<vector_tools.h>
#include<constitutive_tools.h>
#include<micromorphic_tools.h>
#include<micromorphic_linear_elasticity.h>

namespace micromorphicElastoPlasticity{

    typedef micromorphicTools::variableType variableType;
    typedef micromorphicTools::variableVector variableVector;
    typedef micromorphicTools::variableMatrix variableMatrix;

    typedef micromorphicTools::parameterType parameterType;
    typedef micromorphicTools::parameterVector parameterVector;
    typedef micromorphicTools::parameterMatrix parameterMatrix;

    typedef micromorphicTools::constantType constantType;
    typedef micromorphicTools::constantVector constantVector;
    typedef micromorphicTools::constantMatrix constantMatrix;

    typedef errorTools::Node errorNode;
    typedef errorNode* errorOut;

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue );

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &elasticRightCauchyGreenInverse );

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdElasticRCG );

}

#endif
