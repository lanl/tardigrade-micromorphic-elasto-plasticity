/*!
 * tardigrade_micromorphic_elasto_plasticity.h
 *
 * An implementation of a elasto-plastic micromorphic constitutive model 
 * following the derivations of Farhad Shahabi in his dissertation.
 */

#ifndef TARDIGRADE_MICROMORPHIC_ELASTO_PLASTICITY_H
#define TARDIGRADE_MICROMORPHIC_ELASTO_PLASTICITY_H

#include<tardigrade_error_tools.h>
#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>
#include<tardigrade_micromorphic_linear_elasticity.h>
#include<tardigrade_solver_tools.h>
//#include<micromorphic_material_library.h>

namespace tardigradeMicromorphicElastoPlasticity{

    typedef tardigradeMicromorphicTools::variableType variableType;
    typedef tardigradeMicromorphicTools::variableVector variableVector;
    typedef tardigradeMicromorphicTools::variableMatrix variableMatrix;

    typedef tardigradeMicromorphicTools::parameterType parameterType;
    typedef tardigradeMicromorphicTools::parameterVector parameterVector;
    typedef tardigradeMicromorphicTools::parameterMatrix parameterMatrix;

    typedef tardigradeMicromorphicTools::constantType constantType;
    typedef tardigradeMicromorphicTools::constantVector constantVector;
    typedef tardigradeMicromorphicTools::constantMatrix constantMatrix;

    typedef tardigradeErrorTools::Node errorNode;
    typedef errorNode* errorOut;

    struct cout_redirect{
        cout_redirect( std::streambuf * new_buffer)
            : old( std::cout.rdbuf( new_buffer ) )
        { }

        ~cout_redirect( ) {
            std::cout.rdbuf( old );
        }

        private:
            std::streambuf * old;
    };

    struct cerr_redirect{
        cerr_redirect( std::streambuf * new_buffer)
            : old( std::cerr.rdbuf( new_buffer ) )
        { }

        ~cerr_redirect( ) {
            std::cerr.rdbuf( old );
        }

        private:
            std::streambuf * old;
    };

    errorOut computeDruckerPragerInternalParameters( const parameterType &frictionAngle, const parameterType &beta,
                                                     parameterType &A, parameterType &B );

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue );

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdElasticRCG, double tol = 1e-9 );

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdElasticRCG, variableMatrix &d2FdStress2,
                                                           variableMatrix &d2FdStressdElasticRCG, double tol = 1e-9 );

    errorOut computeHigherOrderDruckerPragerYieldEquation( const variableVector &referenceHigherOrderStress,
                                                           const variableVector &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue );

    errorOut computeHigherOrderDruckerPragerYieldEquation( const variableVector &referenceHigherOrderStress,
                                                           const variableVector &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                           variableMatrix &dFdElasticRCG );

    errorOut computeHigherOrderDruckerPragerYieldEquation( const variableVector &referenceHigherOrderStress,
                                                           const variableVector &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                           variableMatrix &dFdElasticRCG, variableMatrix &d2FdStress2,
                                                           variableMatrix &d2FdStressdElasticRCG );

    errorOut computeElasticPartOfDeformation( const variableVector &deformationGradient, const variableVector &microDeformation,
                                              const variableVector &gradientMicroDeformation,
                                              const variableVector &plasticDeformationGradient,
                                              const variableVector &plasticMicroDeformation,
                                              const variableVector &plasticGradientMicroDeformation,
                                              variableVector &inversePlasticDeformationGradient,
                                              variableVector &inversePlasticMicroDeformation,
                                              variableVector &elasticDeformationGradient, variableVector &elasticMicroDeformation,
                                              variableVector &elasticGradientMicroDeformation );

    errorOut computeElasticPartOfDeformation( const variableVector &deformationGradient, const variableVector &microDeformation,
                                              const variableVector &gradientMicroDeformation,
                                              const variableVector &plasticDeformationGradient,
                                              const variableVector &plasticMicroDeformation,
                                              const variableVector &plasticGradientMicroDeformation,
                                              variableVector &elasticDeformationGradient, variableVector &elasticMicroDeformation,
                                              variableVector &elasticGradientMicroDeformation );

    errorOut computeElasticPartOfDeformation( const variableVector &deformationGradient, const variableVector &microDeformation,
                                              const variableVector &gradientMicroDeformation,
                                              const variableVector &plasticDeformationGradient,
                                              const variableVector &plasticMicroDeformation,
                                              const variableVector &plasticGradientMicroDeformation,
                                              variableVector &elasticDeformationGradient, variableVector &elasticMicroDeformation,
                                              variableVector &elasticGradientMicroDeformation,
                                              variableMatrix &dElasticFdF, variableMatrix &dElasticFdPlasticF,
                                              variableMatrix &dElasticChidChi, variableMatrix &dElasticChidPlasticChi,
                                              variableMatrix &dElasticGradChidGradChi, variableMatrix &dElasticGradChidPlasticGradChi,
                                              variableMatrix &dElasticGradChidPlasticF, variableMatrix &dElasticGradChidChi,
                                              variableMatrix &dElasticGradChidPlasticChi );

    errorOut computeElasticDeformationMeasures( const variableVector &elasticDeformationGradient,
                                                const variableVector &elasticMicroDeformation,
                                                const variableVector &elasticGradientMicroDeformation,
                                                variableVector &elasticRightCauchyGreen,
                                                variableVector &elasticMicroRightCauchyGreen,
                                                variableVector &elasticPsi, variableVector &elasticGamma );

    errorOut computeElasticDeformationMeasures( const variableVector &elasticDeformationGradient,
                                                const variableVector &elasticMicroDeformation,
                                                const variableVector &elasticGradientMicroDeformation,
                                                variableVector &elasticRightCauchyGreen,
                                                variableVector &elasticMicroRightCauchyGreen,
                                                variableVector &elasticPsi, variableVector &elasticGamma,
                                                variableMatrix &dElasticRCGdElasticF, variableMatrix &dElasticMicroRCGdElasticChi,
                                                variableMatrix &dElasticPsidElasticF, variableMatrix &dElasticPsidElasticChi,
                                                variableMatrix &dElasticGammadElasticF, variableMatrix &dElasticGammadElasticGradChi );

    errorOut computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &macroPlasticVelocityGradient );

    errorOut computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &macroPlasticVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma );

    errorOut computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma,
                                                  variableMatrix &dPlasticMacroLdElasticRCG,
                                                  variableMatrix &dPlasticMacroLdMacroFlowDirection,
                                                  variableMatrix &dPlasticMacroLdMicroFlowDirection );

    errorOut computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient );

    errorOut computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma );

    errorOut computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma,
                                                  variableMatrix &dPlasticMicroLdElasticMicroRCG,
                                                  variableMatrix &dPlasticMicroLdElasticPsi,
                                                  variableMatrix &dPlasticMicroLdMicroFlowDirection );

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient );

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm );

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL );

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL );

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL,
                                                          variableMatrix &dPlasticMicroGradientLdElasticPsi,
                                                          variableMatrix &dPlasticMicroGradientLdElasticGamma,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection );

    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma, 
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableVector &microGradientFlowDirection,
                                              variableVector &plasticMacroVelocityGradient, variableVector &plasticMicroVelocityGradient,
                                              variableVector &plasticMicroGradientVelocityGradient );

    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma, 
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableVector &microGradientFlowDirection,
                                              variableVector &macroPlasticVelocityGradient, variableVector &microPlasticVelocityGradient,
                                              variableVector &microGradientPlasticVelocityGradient, variableVector &dMacroLpdMacroGamma,
                                              variableVector &dMacroLpdMicroGamma, variableVector &dMicroLpdMicroGamma,
                                              variableVector &dMicroGradientLpdMicroGamma,
                                              variableMatrix &dMicroGradientLpdMicroGradientGamma );

    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma, 
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableVector &microGradientFlowDirection,
                                              variableVector &macroPlasticVelocityGradient, variableVector &microPlasticVelocityGradient,
                                              variableVector &microGradientPlasticVelocityGradient, variableVector &dMacroLpdMacroGamma,
                                              variableVector &dMacroLpdMicroGamma, variableVector &dMicroLpdMicroGamma,
                                              variableVector &dPlasticMicroGradientLdMicroGamma,
                                              variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                              variableMatrix &dPlasticMacroLdElasticRCG,
                                              variableMatrix &dPlasticMacroLdMacroFlowDirection,
                                              variableMatrix &dPlasticMacroLdMicroFlowDirection,
                                              variableMatrix &dPlasticMicroLdElasticMicroRCG,
                                              variableMatrix &dPlasticMicroLdElasticPsi,
                                              variableMatrix &dPlasticMicroLdMicroFlowDirection,
                                              variableMatrix &dPlasticMicroGradientLdElasticMicroRCG,
                                              variableMatrix &dPlasticMicroGradientLdElasticPsi,
                                              variableMatrix &dPlasticMicroGradientLdElasticGamma,
                                              variableMatrix &dPlasticMicroGradientLdMicroFlowDirection,
                                              variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection );

    errorOut evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        const parameterType alpha = 0.5 );

    errorOut evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableMatrix &LHS,
                                        const parameterType alpha = 0.5 );

    errorOut evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                        const parameterType alpha = 0.5 );

    errorOut evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       const parameterType alphaMacro = 0.5,
                                       const parameterType alphaMicro = 0.5,
                                       const parameterType alphaMicroGradient = 0.5 );

    errorOut evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       variableMatrix &dPlasticFdPlasticMacroL,
                                       variableMatrix &dPlasticMicroDeformationdPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMacroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMicroGradientL,
                                       const parameterType alphaMacro = 0.5,
                                       const parameterType alphaMicro = 0.5,
                                       const parameterType alphaMicroGradient = 0.5 );

    errorOut evolveStrainStateVariables( const constantType &Dt, const variableType &currentMacroGamma,
                                         const variableType &currentMicroGamma, const variableVector &currentMicroGradientGamma,
                                         const variableType &currentdMacroGdMacroC, const variableType &currentdMicroGdMicroC,
                                         const variableMatrix &currentdMicroGradientGdMicroGradientC,
                                         const variableType &previousMacroStrainISV, const variableType &previousMicroStrainISV,
                                         const variableVector &previousMicroGradientStrainISV,
                                         const variableType &previousMacroGamma, const variableType &previousMicroGamma,
                                         const variableVector &previousMicroGradientGamma, const variableType &previousdMacroGdMacroC,
                                         const variableType &previousdMicroGdMicroC,
                                         const variableMatrix &previousdMicroGradientGdMicroGradientC,
                                         variableType &currentMacroStrainISV, variableType &currentMicroStrainISV,
                                         variableVector &currentMicroGradientStrainISV,
                                         const parameterType alphaMacro = 0.5,
                                         const parameterType alphaMicro = 0.5,
                                         const parameterType alphaMicroGradient = 0.5 );

    errorOut evolveStrainStateVariables( const constantType &Dt, const variableType &currentMacroGamma,
                                         const variableType &currentMicroGamma, const variableVector &currentMicroGradientGamma,
                                         const variableType &currentdMacroGdMacroC, const variableType &currentdMicroGdMicroC,
                                         const variableMatrix &currentdMicroGradientGdMicroGradientC,
                                         const variableType &previousMacroStrainISV, const variableType &previousMicroStrainISV,
                                         const variableVector &previousMicroGradientStrainISV,
                                         const variableType &previousMacroGamma, const variableType &previousMicroGamma,
                                         const variableVector &previousMicroGradientGamma, const variableType &previousdMacroGdMacroC,
                                         const variableType &previousdMicroGdMicroC,
                                         const variableMatrix &previousdMicroGradientGdMicroGradientC,
                                         variableType &currentMacroStrainISV, variableType &currentMicroStrainISV,
                                         variableVector &currentMicroGradientStrainISV,
                                         variableType &dCurrentMacroISVdCurrentMacroGamma, variableType &dCurrentMacroISVddMacroGdMacroC,
                                         variableType &dCurrentMicroISVdCurrentMicroGamma, variableType &dCurrentMicroISVddMicroGdMicroC,
                                         variableMatrix &dCurrentMicroGradISVdCurrentMicroGradGamma,
                                         variableMatrix &dCurrentMicroGradISVddMicroGradGdMicroGradC,
                                         const parameterType alphaMacro = 0.5,
                                         const parameterType alphaMicro = 0.5,
                                         const parameterType alphaMicroGradient = 0.5 );

    errorOut computeFlowDirections( const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                    const variableVector &referenceHigherOrderStress, const variableType &macroCohesion,
                                    const variableType &microCohesion, const variableVector &microGradientCohesion,
                                    const variableVector &elasticRightCauchyGreen, const parameterVector &macroFlowParameters,
                                    const parameterVector &microFlowParameters, const parameterVector &microGradientFlowParameters,
                                    variableVector &macroFlowDirection, variableVector &microFlowDirection,
                                    variableVector &microGradientFlowDirection, variableType &dGdMacroCohesion,
                                    variableType &dGdMicroCohesion, variableMatrix &dGdMicroGradientCohesion );

    errorOut computeFlowDirections( const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                    const variableVector &referenceHigherOrderStress, const variableType &macroCohesion,
                                    const variableType &microCohesion, const variableVector &microGradientCohesion,
                                    const variableVector &elasticRightCauchyGreen, const parameterVector &macroFlowParameters,
                                    const parameterVector &microFlowParameters, const parameterVector &microGradientFlowParameters,
                                    variableVector &macroFlowDirection, variableVector &microFlowDirection,
                                    variableVector &microGradientFlowDirection, variableType &dGdMacroCohesion,
                                    variableType &dGdMicroCohesion, variableMatrix &dGdMicroGradientCohesion,
                                    variableMatrix &dMacroFlowDirectiondPK2Stress, variableMatrix &dMacroFlowDirectiondElasticRCG,
                                    variableMatrix &dMicroFlowDirectiondReferenceMicroStress,
                                    variableMatrix &dMicroFlowDirectiondElasticRCG,
                                    variableMatrix &dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                    variableMatrix &dMicroGradientFlowDirectiondElasticRCG );

    errorOut evaluateYieldFunctions( const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                     const variableVector &referenceHigherOrderStress, const variableType &macroCohesion,
                                     const variableType &microCohesion, const variableVector &microGradientCohesion,
                                     const variableVector &elasticRightCauchyGreen,
                                     const parameterVector &macroYieldParameters, const parameterVector &microYieldParameters,
                                     const parameterVector &microGradientYieldParameters, variableVector &yieldFunctionValues
#ifdef DEBUG_MODE
                                     , tardigradeSolverTools::debugMap &DEBUG
#endif
                                   );

    errorOut evaluateYieldFunctions( const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                     const variableVector &referenceHigherOrderStress, const variableType &macroCohesion,
                                     const variableType &microCohesion, const variableVector &microGradientCohesion,
                                     const variableVector &elasticRightCauchyGreen,
                                     const parameterVector &macroYieldParameters, const parameterVector &microYieldParameters,
                                     const parameterVector &microGradientYieldParameters, variableVector &yieldFunctionValues,
                                     variableVector &dMacroFdPK2, variableType &dMacroFdMacroC, variableVector &dMacroFdElasticRCG,
                                     variableVector &dMicroFdSigma, variableType &dMicroFdMicroC, variableVector &dMicroFdElasticRCG,
                                     variableMatrix &dMicroGradientFdM, variableMatrix &dMicroGradientFdMicroGradientC,
                                     variableMatrix &dMicroGradientFdElasticRCG
#ifdef DEBUG_MODE
                                     , tardigradeSolverTools::debugMap &DEBUG
#endif
                                    );

    errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation );

    errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation,
                                                     variableMatrix &dFdGradU, variableMatrix &dChidPhi,
                                                     variableMatrix &dGradChidGradPhi );

    errorOut extractMaterialParameters( const std::vector< double > &fparams,
                                        parameterVector &macroHardeningParameters, parameterVector &microHardeningParameters,
                                        parameterVector &microGradientHardeningParameters,
                                        parameterVector &macroFlowParameters, parameterVector &microFlowParameters,
                                        parameterVector &microGradientFlowParameters,
                                        parameterVector &macroYieldParameters, parameterVector &microYieldParameters,
                                        parameterVector &microGradientYieldParameters,
                                        parameterVector &Amatrix, parameterVector &Bmatrix,
                                        parameterVector &Cmatrix, parameterVector &Dmatrix,
                                        constantType &alphaMacro, constantType &alphaMicro, constantType &alphaMicroGradient,
                                        constantType &relativeTolerance, constantType &absoluteTolerance );

    errorOut extractStateVariables( std::vector< double > &SDVS,
                                    variableType &previousMacroStrainISV, variableType &previousMicroStrainISV,
                                    variableVector &previousMicroGradientStrainISV,
                                    variableType &previousMacroGamma, variableType &previousMicroGamma,
                                    variableVector &previousMicroGradientGamma,
                                    variableVector &previousPlasticDeformationGradient,
                                    variableVector &previousPlasticMicroDeformation,
                                    variableVector &previousPlasticGradientMicroDeformation );

    errorOut computeCohesion( const variableType &macroStrainISV, const variableType &microStrainISV,
                              const variableVector &microGradientStrainISV,
                              const parameterVector &macroHardeningParameters, const parameterVector &microHardeningParameters,
                              const parameterVector &microGradientHardeningParameters,
                              variableType &macroCohesion, variableType &microCohesion,
                              variableVector &microGradientCohesion );

    errorOut computeCohesion( const variableType &macroStrainISV, const variableType &microStrainISV,
                              const variableVector &microGradientStrainISV,
                              const parameterVector &macroHardeningParameters, const parameterVector &microHardeningParameters,
                              const parameterVector &microGradientHardeningParameters,
                              variableType &macroCohesion, variableType &microCohesion,
                              variableVector &microGradientCohesion,
                              variableType &dMacroCdMacroStrainISV, variableType &dMicroCdMicroStrainISV,
                              variableMatrix &dMicroGradientCdMicroGradientStrainISV );

    errorOut computePlasticDeformationResidual( const tardigradeSolverTools::floatVector &x, const tardigradeSolverTools::floatMatrix &floatArgs,
                                                const tardigradeSolverTools::intMatrix &intArgs, tardigradeSolverTools::floatVector &residual,
                                                tardigradeSolverTools::floatMatrix &jacobian, tardigradeSolverTools::floatMatrix &floatOuts,
                                                tardigradeSolverTools::intMatrix &intOuts
#ifdef DEBUG_MODE
                                                , tardigradeSolverTools::debugMap &DEBUG
#endif
                                              );

    errorOut computePlasticMultiplierResidual( const tardigradeSolverTools::floatVector &x, const tardigradeSolverTools::floatMatrix &floatArgs,
                                               const tardigradeSolverTools::intMatrix &intArgs, tardigradeSolverTools::floatVector &residual,
                                               tardigradeSolverTools::floatMatrix &jacobian, tardigradeSolverTools::floatMatrix &floatOuts,
                                               tardigradeSolverTools::intMatrix &intOuts
#ifdef DEBUG_MODE
                                               , tardigradeSolverTools::debugMap &DEBUG
#endif
                                             );

    errorOut computePlasticMultiplierLagrangian( const tardigradeSolverTools::floatVector &x, const tardigradeSolverTools::floatMatrix &floatArgs,
                                                 const tardigradeSolverTools::intMatrix &intArgs, tardigradeSolverTools::floatType &lagrangian,
                                                 tardigradeSolverTools::floatVector &jacobian, tardigradeSolverTools::floatMatrix &floatOuts,
                                                 tardigradeSolverTools::intMatrix &intOuts
#ifdef DEBUG_MODE
                                                 , tardigradeSolverTools::debugMap &DEBUG
#endif
                                               ); 
 
    int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                        const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                        const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                        const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                        std::vector< double > &SDVS,
                        const std::vector< double > &current_ADD_DOF,
                        const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                        const std::vector< double > &previous_ADD_DOF,
                        const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                        std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                        std::vector< std::vector< double > > &ADD_TERMS,
                        std::string &output_message
#ifdef DEBUG_MODE
                        , tardigradeSolverTools::homotopyMap &DEBUG
#endif
                      );

    int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                        const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                        const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                        const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                        std::vector< double > &SDVS,
                        const std::vector< double > &current_ADD_DOF,
                        const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                        const std::vector< double > &previous_ADD_DOF,
                        const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                        std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                        std::vector< std::vector< double > > &DPK2Dgrad_u,   std::vector< std::vector< double > > &DPK2Dphi,
                        std::vector< std::vector< double > > &DPK2Dgrad_phi,
                        std::vector< std::vector< double > > &DSIGMADgrad_u, std::vector< std::vector< double > > &DSIGMADphi,
                        std::vector< std::vector< double > > &DSIGMADgrad_phi,
                        std::vector< std::vector< double > > &DMDgrad_u,     std::vector< std::vector< double > > &DMDphi,
                        std::vector< std::vector< double > > &DMDgrad_phi,
                        std::vector< std::vector< double > > &ADD_TERMS,
                        std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS,
                        std::string &output_message
#ifdef DEBUG_MODE
                        , tardigradeSolverTools::homotopyMap &DEBUG
#endif
                      );

}

#endif
