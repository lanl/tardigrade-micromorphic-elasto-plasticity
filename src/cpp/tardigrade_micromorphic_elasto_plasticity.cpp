/*!
 * tardigrade_micromorphic_elasto_plasticity.cpp
 *
 * An implementation of a elasto-plastic micromorphic constitutive model 
 * following the derivations of Farhad Shahabi in his dissertation.
 */

#include<tardigrade_micromorphic_elasto_plasticity.h>

namespace tardigradeMicromorphicElastoPlasticity{

    errorOut computeDruckerPragerInternalParameters( const parameterType &frictionAngle, const parameterType &beta,
                                                     parameterType &A, parameterType &B ){
        /*!
         * Compute the Drucker-Prager internal parameters
         *
         * :param const parameterType &frictionAngle: The material friction angle ( 0 < frictionAngle < pi / 2 );
         * :param const parameterType &beta: The beta parameter.
         * :param parameterType &A: The A parameter.
         * :param parameterType &B: The B parameter.
         */

        //Make sure the parameters are within bounds
        if ( ( 0 > frictionAngle ) || ( frictionAngle > 1.570796 ) ){
            return new errorNode( "computeSecondOrderDruckerPragerYieldEquation",
                                  "The friction angle must between 0 and pi / 2" );
        }

        if ( abs( beta ) > 1 ){
            return new errorNode( "computeSecondOrderDruckerPragerYieldEquation",
                                  "Beta must be between -1 and 1" );
        }

        //Compute the parameters
        parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );

        A = betaAngle * std::cos( frictionAngle );

        B = betaAngle * std::sin( frictionAngle );

        return NULL;
    }

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue ){
        /*!
         * Compute the second-order Drucker Prager Yield equation
         *
         * F = ||dev ( stressMeasure ) || - \left( A^{\phi} \bar{c} - B^{\phi} \bar{p} \right) \leq 0
         * 
         * || dev ( stressMeasure ) || = \sqrt{ dev( referenceStressMeasure ) : dev( referenceStressMeasure ) }
         *  dev( referenceStressMeasure ) : dev( referenceStressMeasure ) = dev( referenceStressMeasure )_{IJ} dev( referenceStressMeasure )_{IJ}
         *  dev( referenceStressMeasure )_{IJ} = referenceStressMeasure_{IJ} - \bar{p} elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration
         * :param const variableType &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableType &yieldValue: The yield value.
         */

        parameterType AAngle, BAngle;
        errorOut error = computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation",
                                             "Error in computation of the Drucker-Prager internal parameters" );
            result->addNext( error );
            return result;
        }

        //Compute the decomposition of the stress
        variableType pressure;
        variableVector deviatoricReferenceStress;
        
        error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation",
                                             "Error in computation of second-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        return NULL;
    }

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdElasticRCG, double tol ){
        /*!
         * Compute the second-order Drucker Prager Yield equation
         *
         * F = ||dev ( stressMeasure ) || - \left( A^{\phi} \bar{c} - B^{\phi} \bar{p} \right) \leq 0
         * 
         * || dev ( stressMeasure ) || = \sqrt{ dev( referenceStressMeasure ) : dev( referenceStressMeasure ) }
         *  dev( referenceStressMeasure ) : dev( referenceStressMeasure ) = dev( referenceStressMeasure )_{IJ} dev( referenceStressMeasure )_{IJ}
         *  dev( referenceStressMeasure )_{IJ} = referenceStressMeasure_{IJ} - \bar{p} elasticRightCauchyGreen_{IJ}^{-1}
         *
         *  Also compute the Jacobians
         * \frac{ \partial F }{ \partial stressMeasure_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure ) \frac{ \partial dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial stressMeasure_{IJ} }
         * \frac{ \partial F }{ \partial \bar{c} } = -A^{\phi}
         * \frac{ \partial F }{ \partial C_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure ) \frac{ \partial dev( stressMeasure )_{AB} }{ \partial C_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial C_{IJ} }
         *
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration
         * :param const variableType &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableType &yieldValue: The yield value.
         * :param variableVector &dFdStress: The Jacobian of the yield surface w.r.t. the stress measure.
         * :param variableType &dFdc: The Jacobian of the yield surface w.r.t. the cohesion.
         * :param variableVector &dFdElasticRCG: The Jacobian of the yield surface w.r.t. the elastic 
         * :param double tol: The tolerance used to prevent nans in the Jacobians
         */

        parameterType AAngle, BAngle;
        errorOut error = computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of the Drucker-Prager internal parameters" );
            result->addNext( error );
            return result;
        }

        //Compute the decomposition of the stress
        variableType pressure;
        variableVector deviatoricReferenceStress;

        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableVector dPressuredStress, dPressuredRCG;
        
        error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of second-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Evaluate the jacobians
        variableVector devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );

        dFdStress = tardigradeVectorTools::Tdot( dDevStressdStress, devStressDirection )
                  + BAngle * dPressuredStress;

        dFdc = - AAngle;

        dFdElasticRCG = tardigradeVectorTools::Tdot( dDevStressdRCG, devStressDirection )
                      + BAngle * dPressuredRCG;

        return NULL;
    }

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdElasticRCG, variableMatrix &d2FdStress2,
                                                           variableMatrix &d2FdStressdElasticRCG, double tol ){
        /*!
         * Compute the second-order Drucker Prager Yield equation
         *
         * F = ||dev ( stressMeasure ) || - \left( A^{\phi} \bar{c} - B^{\phi} \bar{p} \right) \leq 0
         * 
         * || dev ( stressMeasure ) || = \sqrt{ dev( referenceStressMeasure ) : dev( referenceStressMeasure ) }
         *  dev( referenceStressMeasure ) : dev( referenceStressMeasure ) = dev( referenceStressMeasure )_{IJ} dev( referenceStressMeasure )_{IJ}
         *  dev( referenceStressMeasure )_{IJ} = referenceStressMeasure_{IJ} - \bar{p} elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         *  Also compute the Jacobians
         * \frac{ \partial F }{ \partial stressMeasure_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure ) \frac{ \partial dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial stressMeasure_{IJ} }
         * \frac{ \partial F }{ \partial \bar{c} } = -A^{\phi}
         * \frac{ \partial F }{ \partial C_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure )_{AB} }{ \partial C_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial C_{IJ} }
         *
         *  The second deriatives of \frac{ \partial F }{ \partial \stressMeasure_{IJ}  } are
         *  \frac{ \partial^2 F }{ \partial stressMeasure_{IJ} \partial stressMeasure_{KL} } = \frac{ \partial^2 || dev( stressMeasure ) || }{ \partial dev( stressMeasure )_{AB} \partial dev( stressMeasure )_{CD} } \frac{ \partial dev( stressMeasure )_{AB} } { \partial stressMeasure_{IJ} } \frac{ \partial dev( stressMeasure )_{CD} } { \partial stressMeasure_{KL} } + \frac{ dev ( stressMeasure )_{AB} }{ || dev( stressMeasure ) || } \frac{ \partial^2 dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} \partial stressMeasure_{KL} }
         *  \frac{ \partial^2 F }{ \partial stressMeasure_{IJ} \partial RCG_{KL} } = \frac{ \partial^2 || dev( stressMeasure ) || }{ \partial dev( stressMeasure )_{AB} \partial dev( stressMeasure )_{CD} } \frac{ \partial dev( stressMeasure )_{AB} } { \partial stressMeasure_{IJ} } \frac{ \partial dev( stressMeasure )_{CD} } { \partial C_{KL} } + \frac{ dev ( stressMeasure )_{AB} }{ || dev( stressMeasure ) || } \frac{ \partial^2 dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} \partial C_{KL} } + B^{\phi} \frac{ \partial^2 \bar{p} }{ \partial stressMeasure_{IJ} \partial C_{KL} }
         *  
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration
         * :param const variableType &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableType &yieldValue: The yield value.
         * :param variableVector &dFdStress: The Jacobian of the yield surface w.r.t. the stress measure.
         * :param variableType &dFdc: The Jacobian of the yield surface w.r.t. the cohesion.
         * :param variableVector &dFdElasticRCG: The Jacobian of the yield surface w.r.t. the elastic 
         *     right Cauchy-Green deformation tensor.
         * :param variableMatrix &d2FdStress2: The second derivative of the flow direction w.r.t. the stress. This 
         *     is useful if one is using this expression as the flow potential and wants the jacobian of the flow direction \frac{ \partial G }{\partial Stress_{IJ} }
         * :param variableMatrix &d2FdStressdElasticRCG: The second derivative of the flow direction w.r.t. the stress.
         * :param double tol: The tolerance used to prevent nans in the Jacobians
         */
        //Assume 3D
        unsigned int dim = 3;

        parameterType AAngle, BAngle;
        errorOut error = computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation (second order jacobian)",
                                             "Error in computation of the Drucker-Prager internal parameters" );
            result->addNext( error );
            return result;
        }

        //Compute the decomposition of the stress
        variableType pressure;
        variableVector deviatoricReferenceStress;
        
        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableVector dPressuredStress, dPressuredRCG;

        variableMatrix d2DevStressdStressdRCG, d2PressuredStressdRCG;

        error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation (second order jacobian)",
                                             "Error in computation of second-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Evaluate the jacobians
        variableVector devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );

        dFdStress = tardigradeVectorTools::Tdot( dDevStressdStress, devStressDirection )
                  + BAngle * dPressuredStress;

        dFdc = - AAngle;

        dFdElasticRCG = tardigradeVectorTools::Tdot( dDevStressdRCG, devStressDirection )
                      + BAngle * dPressuredRCG;

        //Evaluate the second-order jacobians
        constantMatrix EYE = tardigradeVectorTools::eye< constantType >( dim * dim );
        variableMatrix dDevStressDirectiondDevStress = ( EYE - tardigradeVectorTools::dyadic( devStressDirection, devStressDirection ) ) / ( normDevStress + tol );

        d2FdStress2 = tardigradeVectorTools::Tdot( dDevStressdStress, tardigradeVectorTools::dot( dDevStressDirectiondDevStress, dDevStressdStress ) );

        d2FdStressdElasticRCG = tardigradeVectorTools::Tdot( tardigradeVectorTools::dot( dDevStressDirectiondDevStress, dDevStressdStress ), dDevStressdRCG )
                              + BAngle * d2PressuredStressdRCG;

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int A = 0; A < dim; A++ ){
                            for ( unsigned int B = 0; B < dim; B++ ){
                                d2FdStressdElasticRCG[ dim * I + J ][ dim * K + L ] += devStressDirection[ dim * A + B ] * d2DevStressdStressdRCG[ dim * A + B ][ dim * dim * dim * I + dim * dim * J + dim * K + L ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeHigherOrderDruckerPragerYieldEquation( const variableVector &referenceHigherOrderStress, 
                                                           const variableVector &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue ){
        /*!
         * Compute the higher-order Drucker Prager Yield equation
         *
         * F_K = ||dev ( M ) ||_K - \left( A^{\phi} \bar{c}_K - B^{\phi} \bar{p}_K \right) \leq 0
         * 
         * || dev ( stressMeasure ) ||_K = \sqrt{ dev( M )_{IJK} : dev( M )_{IJK} }
         * where the K's aren't summed.
         *  dev( M )_{IJK} = M_{IJK} - \bar{p}_K elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} M_{IJK}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         * :param const variableVector &referenceHigherOrderStress: The higher-order stress in the reference configuration
         * :param const variableVector &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableVector &yieldValue: The yield value.
         */

        parameterType AAngle, BAngle;
        errorOut error = computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation",
                                             "Error in computation of the Drucker-Prager internal parameters" );
            result->addNext( error );
            return result;
        }

        //Compute the decomposition of the stress
        variableVector pressure;
        variableVector deviatoricReferenceStress;
        
        error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( referenceHigherOrderStress,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation",
                                             "Error in computation of higher-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableVector normDevStress;
        error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation",
                                             "Error in computation of the deviatoric higher-order stress norm" );
            result->addNext( error );
            return result;
        }

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        return NULL;
    }

    errorOut computeHigherOrderDruckerPragerYieldEquation( const variableVector &referenceHigherOrderStress, 
                                                           const variableVector &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                           variableMatrix &dFdElasticRCG ){
        /*!
         * Compute the higher-order Drucker Prager Yield equation
         *
         * F_K = ||dev ( M ) ||_K - \left( A^{\phi} \bar{c}_K - B^{\phi} \bar{p}_K \right) \leq 0
         * 
         * || dev ( stressMeasure ) ||_K = \sqrt{ dev( M )_{IJK} : dev( M )_{IJK} }
         * where the K's aren't summed.
         *  dev( M )_{IJK} = M_{IJK} - \bar{p}_K elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} M_{IJK}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         *  Also computes the Jacobians
         *
         * :param const variableVector &referenceHigherOrderStress: The higher-order stress in the reference configuration
         * :param const variableVector &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableVector &yieldValue: The yield value.
         * :param variableMatrix &dFdStress: The Jacobian of the yield function w.r.t. the reference higher order stress.
         * :param variableMatrix &dFdc: The Jacobian of the yield function w.r.t. the cohesion.
         * :param variableMatrix &dFdElasticRCG: The Jacobian of the yield function w.r.t. the elastic right Cauchy-Green
         *     deformation tensor.
         */

        parameterType AAngle, BAngle;
        errorOut error = computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of the Drucker-Prager internal parameters" );
            result->addNext( error );
            return result;
        }

        //Compute the decomposition of the stress
        variableVector pressure;
        variableVector deviatoricReferenceStress;

        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableMatrix dPressuredStress, dPressuredRCG;
        
        error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( referenceHigherOrderStress,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of higher-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableVector normDevStress;
        variableMatrix dNormDevStressdDevStress;
        error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress, dNormDevStressdDevStress );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of the deviatoric higher-order stress norm" );
            result->addNext( error );
            return result;
        }

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Construct the Jacobians
        dFdStress = tardigradeVectorTools::dot( dNormDevStressdDevStress, dDevStressdStress )
                  + BAngle * dPressuredStress;

        dFdc = -AAngle * tardigradeVectorTools::eye< constantType >( cohesion.size() );

        dFdElasticRCG = tardigradeVectorTools::dot( dNormDevStressdDevStress, dDevStressdRCG )
                      + BAngle * dPressuredRCG;

        return NULL;
    }

    errorOut computeHigherOrderDruckerPragerYieldEquation( const variableVector &referenceHigherOrderStress, 
                                                           const variableVector &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                           variableMatrix &dFdElasticRCG, variableMatrix &d2FdStress2,
                                                           variableMatrix &d2FdStressdElasticRCG ){
        /*!
         * Compute the higher-order Drucker Prager Yield equation
         *
         * F_K = ||dev ( M ) ||_K - \left( A^{\phi} \bar{c}_K - B^{\phi} \bar{p}_K \right) \leq 0
         * 
         * || dev ( stressMeasure ) ||_K = \sqrt{ dev( M )_{IJK} : dev( M )_{IJK} }
         * where the K's aren't summed.
         *  dev( M )_{IJK} = M_{IJK} - \bar{p}_K elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} M_{IJK}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         *  Also computes the Jacobians
         *
         * :param const variableVector &referenceHigherOrderStress: The higher-order stress in the reference configuration
         * :param const variableVector &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableVector &yieldValue: The yield value.
         * :param variableMatrix &dFdStress: The Jacobian of the yield function w.r.t. the reference higher order stress.
         * :param variableMatrix &dFdc: The Jacobian of the yield function w.r.t. the cohesion.
         * :param variableMatrix &dFdElasticRCG: The Jacobian of the yield function w.r.t. the elastic right Cauchy-Green
         *     deformation tensor.
         * :param variableMatrix &d2FdStress2: The second order Jacobian of the yield function w.r.t. the reference 
         *     higher order stress.
         * :param variableMatrix &d2FdStressdElasticRCG: The second order Jacobian of the yield function w.r.t. the 
         *     reference higher order stress and the elastic right Cauchy-Green Deformation metric.
         */

        //Assume 3D
        unsigned int dim = 3;

        parameterType AAngle, BAngle;
        errorOut error = computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (second order jacobian)",
                                             "Error in computation of the Drucker-Prager internal parameters" );
            result->addNext( error );
            return result;
        }

        //Compute the decomposition of the stress
        variableVector pressure;
        variableVector deviatoricReferenceStress;

        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableMatrix dPressuredStress, dPressuredRCG;

        variableMatrix d2DevStressdStressdRCG, d2PressuredStressdRCG;
        
        error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( referenceHigherOrderStress,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (second order jacobian)",
                                             "Error in computation of higher-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableVector normDevStress;
        variableMatrix dNormDevStressdDevStress;
        variableMatrix d2NormDevStressdDevStress2;
        error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress,
                                                                 dNormDevStressdDevStress,
                                                                 d2NormDevStressdDevStress2 );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (second order jacobian)",
                                             "Error in computation of the deviatoric higher-order stress norm" );
            result->addNext( error );
            return result;
        }

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Construct the Jacobians
        dFdStress = tardigradeVectorTools::dot( dNormDevStressdDevStress, dDevStressdStress )
                  + BAngle * dPressuredStress;

        dFdc = -AAngle * tardigradeVectorTools::eye< constantType >( cohesion.size() );

        dFdElasticRCG = tardigradeVectorTools::dot( dNormDevStressdDevStress, dDevStressdRCG )
                      + BAngle * dPressuredRCG;

        //Construct the second-order jacobians
        d2FdStress2 = variableMatrix( dim, variableVector( dim * dim * dim * dim * dim * dim, 0 ) );
        d2FdStressdElasticRCG = tardigradeVectorTools::dot( dNormDevStressdDevStress, d2DevStressdStressdRCG )
                              + BAngle * d2PressuredStressdRCG;

        for ( unsigned int K = 0; K < 3; K++ ){
            for ( unsigned int L = 0; L < 3; L++ ){
                for ( unsigned int M = 0; M < 3; M++ ){
                    for ( unsigned int N = 0; N < 3; N++ ){
                        for ( unsigned int O = 0; O < 3; O++ ){
                            for ( unsigned int P = 0; P < 3; P++ ){
                                for ( unsigned int Q = 0; Q < 3; Q++ ){
                                    for ( unsigned int A = 0; A < 3; A++ ){
                                        for ( unsigned int B = 0; B < 3; B++ ){
                                            for ( unsigned int C = 0; C < 3; C++ ){
                                                for ( unsigned int D = 0; D < 3; D++ ){
                                                    for ( unsigned int E = 0; E < 3; E++ ){
                                                        d2FdStressdElasticRCG[ K ][ dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * N + dim * O + P ]
                                                            += d2NormDevStressdDevStress2[ K ][ dim * dim * dim * dim * dim * Q + dim * dim * dim * dim * A + dim * dim * dim * B + dim * dim * C + dim * D + E ]
                                                             * dDevStressdStress[ dim * dim * Q + dim * A + B ][ dim * dim * L + dim * M + N ]
                                                             * dDevStressdRCG[ dim * dim * C + dim * D + E ][ dim * O + P ];
                                                        for ( unsigned int F = 0; F < 3; F++ ){
                                                            d2FdStress2[ K ][ dim * dim * dim * dim * dim * L + dim * dim * dim * dim * M + dim * dim * dim * N + dim * dim * O + dim * P + Q ]
                                                                += d2NormDevStressdDevStress2[ K ][ dim * dim * dim * dim * dim * A + dim * dim * dim * dim * B + dim * dim * dim * C + dim * dim * D + dim * E + F ]
                                                                 * dDevStressdStress[ dim * dim * A + dim * B + C ][ dim * dim * L + dim * M + N ]
                                                                 * dDevStressdStress[ dim * dim * D + dim * E + F ][ dim * dim * O + dim * P + Q ];
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeElasticPartOfDeformation( const variableVector &deformationGradient, const variableVector &microDeformation,
                                              const variableVector &gradientMicroDeformation,
                                              const variableVector &plasticDeformationGradient,
                                              const variableVector &plasticMicroDeformation,
                                              const variableVector &plasticGradientMicroDeformation,
                                              variableVector &inversePlasticDeformationGradient,
                                              variableVector &inversePlasticMicroDeformation,
                                              variableVector &elasticDeformationGradient, variableVector &elasticMicroDeformation, 
                                              variableVector &elasticGradientMicroDeformation ){
        /*!
         * Compute the elastic parts of the various deformation measures.
         *
         * F_{i\bar{I}}^e = F_{iI} F_{I \bar{I}}^{p, -1}
         * \chi_{i\bar{I}}^e = \chi_{iI} \chi_{I \bar{I}}^{p, -1}
         * \chi_{ k\bar{ K },\bar{ L } } = \left( \chi_{ kK, L } F_{ L \bar{ L } }^{p, -1} - \chi_{ k \bar{ A } }^e \chi_{ \bar{ A } K,\bar{ L } }^{ p } \right) \chi_{ K \bar{ K } }^{p,-1}
         *
         * :param const variableVector &deformationGradient: The macroscale deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation tensor $\chi$
         * :param const variableVector &gradientMicroDeformation: The gradient of the micro-deformation tensor
         *     $\chi$ with respect to the reference configuration.
         * :param const variableVector &plasticDeformationGradient: The plastic part of the macroscale 
         *     deformation gradient.
         * :param const variableVector &plasticMicroDeformation: The plastic part of the micro-deformation tensor
         *     $\chi$.
         * :param const variableVector &plasticGradientMicroDeformation: The plastic part of the gradient of the 
         *     micro-deformation tensor $\chi$ with respect to the intermediate configuration.
         * :param variableVector &inversePlasticDeformationGradient: The inverse of the plastic part of the macro deformation gradient.
         * :param variableVector &inversePlasticMicroDeformation: The inverse of the plastic part of the micro deformation.
         * :param variableVector &elasticDeformationGradient: The elastic part of the macroscale deformation gradient.
         * :param variableVector &elasticMicroDeformation: The elastic part of the micro-deformation tensor $\chi$
         * :param variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the micro-deformation
         *     tensor $\chi$ w.r.t. the reference configuration.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The deformation gradient is not the correct size (3D)" );
        }

        if ( microDeformation.size() != dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The micro-deformation is not the correct size (3D)" );
        }

        if ( gradientMicroDeformation.size() != dim * dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The gradient of the micro-deformation is not the correct size (3D)" );
        }

        if ( plasticDeformationGradient.size() != dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The plastic deformation gradient is not the correct size (3D)" );
        }

        if ( plasticMicroDeformation.size() != dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The plastic micro-deformation is not the correct size (3D)" );
        }

        if ( plasticGradientMicroDeformation.size() != dim * dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The plastic gradient of the micro-deformation is not the correct size (3D)" );
        }

        //Compute the inverses of the plastic deformation gradient and micro-deformation
        inversePlasticDeformationGradient = tardigradeVectorTools::inverse( plasticDeformationGradient, dim, dim );

        inversePlasticMicroDeformation = tardigradeVectorTools::inverse( plasticMicroDeformation, dim, dim );

        //Assemble the elastic parts of the deformation measures
        elasticDeformationGradient = tardigradeVectorTools::matrixMultiply( deformationGradient, inversePlasticDeformationGradient,
                                                                  dim, dim, dim, dim );

        elasticMicroDeformation = tardigradeVectorTools::matrixMultiply( microDeformation, inversePlasticMicroDeformation,
                                                               dim, dim, dim, dim );

        elasticGradientMicroDeformation = variableVector( dim * dim * dim, 0 );
        for ( unsigned int k = 0; k < dim; k++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int Ab = 0; Ab < dim; Ab++ ){
                            elasticGradientMicroDeformation[ dim * dim * k + dim * Kb + Lb ]
                                += ( gradientMicroDeformation[ dim * dim * k + dim * K + Ab ]
                                 *   inversePlasticDeformationGradient[ dim * Ab + Lb ]
                                 -   elasticMicroDeformation[ dim * k + Ab ]
                                 *   plasticGradientMicroDeformation[ dim * dim * Ab + dim * K + Lb ] )
                                 * inversePlasticMicroDeformation[ dim * K + Kb ];
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeElasticPartOfDeformation( const variableVector &deformationGradient, const variableVector &microDeformation,
                                              const variableVector &gradientMicroDeformation,
                                              const variableVector &plasticDeformationGradient,
                                              const variableVector &plasticMicroDeformation,
                                              const variableVector &plasticGradientMicroDeformation,
                                              variableVector &elasticDeformationGradient, variableVector &elasticMicroDeformation, 
                                              variableVector &elasticGradientMicroDeformation ){
        /*!
         * Compute the elastic parts of the various deformation measures.
         * F_{i\bar{I}}^e = F_{iI} F_{I \bar{I}}^{p, -1}
         * \chi_{i\bar{I}}^e = \chi_{iI} \chi_{I \bar{I}}^{p, -1}
         * \chi_{ k\bar{ K },\bar{ L } } = \left( \chi_{ kK, L } F_{ L \bar{ L } }^{p, -1} - \chi_{ k \bar{ A } }^e \chi_{ \bar{ A } K,\bar{ L } }^{ p } \right) \chi_{ K \bar{ K } }^{p,-1}
         *
         * :param const variableVector &deformationGradient: The macroscale deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation tensor $\chi$
         * :param const variableVector &gradientMicroDeformation: The gradient of the micro-deformation tensor
         *     $\chi$ with respect to the reference configuration.
         * :param const variableVector &plasticDeformationGradient: The plastic part of the macroscale 
         *     deformation gradient.
         * :param const variableVector &plasticMicroDeformation: The plastic part of the micro-deformation tensor
         *     $\chi$.
         * :param const variableVector &plasticGradientMicroDeformation: The plastic part of the gradient of the 
         *     micro-deformation tensor $\chi$ with respect to the intermediate configuration.
         * :param variableVector &elasticDeformationGradient: The elastic part of the macroscale deformation gradient.
         * :param variableVector &elasticMicroDeformation: The elastic part of the micro-deformation tensor $\chi$
         * :param variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the micro-deformation
         *     tensor $\chi$ w.r.t. the reference configuration.
         */

        variableVector inversePlasticDeformationGradient, inversePlasticMicroDeformation;

        return computeElasticPartOfDeformation( deformationGradient, microDeformation, gradientMicroDeformation,
                                                plasticDeformationGradient, plasticMicroDeformation, plasticGradientMicroDeformation,
                                                inversePlasticDeformationGradient, inversePlasticMicroDeformation,
                                                elasticDeformationGradient, elasticMicroDeformation, elasticGradientMicroDeformation );
    }

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
                                              variableMatrix &dElasticGradChidPlasticChi ){
        /*!
         * Compute the elastic parts of the various deformation measures.
         * F_{i\bar{I}}^e = F_{iI} F_{I \bar{I}}^{p, -1}
         * \chi_{i\bar{I}}^e = \chi_{iI} \chi_{I \bar{I}}^{p, -1}
         * \chi_{ k\bar{ K },\bar{ L } } = \left( \chi_{ kK, L } F_{ L \bar{ L } }^{ p, -1}  - \chi_{ k \bar{ A } }^e \chi_{ \bar{ A } K,\bar{ L } }^{ p } \right) \chi_{ K \bar{ K } }^{p,-1}
         *
         * Also compute the Jacobians
         *
         * \frac{ \partial F_{i\bar{I}}^e }{ \partial F_{nN} } = \delta_{in} F_{N \bar{I} }^{p, -1}
         * \frac{ \partial F_{i\bar{I}}^e }{ \partial F_{\bar{N} N}^p } = -F_{iI} F_{I \bar{N}}^{p, -1} F_{N \bar{I} }^{p, -1}
         * \frac{ \partial \chi_{i\bar{I}}^e }{ \partial \chi_{nN} } = \delta_{in} \chi_{N \bar{I} }^{p, -1}
         * \frac{ \partial \chi_{i\bar{I}}^e }{ \partial \chi_{\bar{N} N}^p } = -\chi_{iI} \chi_{I \bar{N}}^{p, -1} \chi_{N \bar{I} }^{p, -1}
         * \frac{ \partial \chi_{k \bar{ K }, \bar{ L } } }{ \partial F_{ \bar{N} N }^{ p } } = \chi_{ kK, L } F_{ L \bar{ N } }^{p, -1} F_{ N \bar{ L } }^{p, -1}
         * \frac{ \partial \chi_{k\bar{ K },\bar{ L } } }{ \partial \chi_{ nN } } = - \frac{ \partial \chi_{k \bar{ A } }^e }{ \partial \chi_{ nN } } \chi_{ \bar{ A } K,\bar{ L } }^{ p } \chi_{ K \bar{ K } }^{p,-1}
         * \frac{ \partial \chi_{k\bar{ K },\bar{ L } } }{ \partial \chi_{ \bar{ N } N } } = \chi_{ k \bar{ A } }^e \chi_{ \bar{ A } K,\bar{ L } }^{ p } \chi_{ K \bar{ N } }^{p,-1} \chi_{ N \bar{ K } }^{p,-1}
         *
         * :param const variableVector &deformationGradient: The macroscale deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation tensor $\chi$
         * :param const variableVector &gradientMicroDeformation: The gradient of the micro-deformation tensor
         *     $\chi$ with respect to the reference configuration.
         * :param const variableVector &plasticDeformationGradient: The plastic part of the macroscale 
         *     deformation gradient.
         * :param const variableVector &plasticMicroDeformation: The plastic part of the micro-deformation tensor
         *     $\chi$.
         * :param const variableVector &plasticGradientMicroDeformation: The plastic part of the gradient of the 
         *     micro-deformation tensor $\chi$ with respect to the intermediate configuration.
         * :param variableVector &elasticDeformationGradient: The elastic part of the macroscale deformation gradient.
         * :param variableVector &elasticMicroDeformation: The elastic part of the micro-deformation tensor $\chi$
         * :param variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the micro-deformation
         *     tensor $\chi$ w.r.t. the intermediate configuration.
         * :param variableMatrix &dElasticFdF: The Jacobian of the elastic part of the deformation gradient w.r.t. the deformation 
         *     gradient.
         * :param variableMatrix &dElasticFdPlasticF: The Jacobian of the elastic part of the deformation gradient w.r.t. the 
         *     plastic part of the deformation gradient.
         * :param variableMatrix &dElasticChidChi: The Jacobian of the elastic part of the micro-deformation w.r.t. the
         *     micro deformation.
         * :param variableMatrix &dElasticChidPlasticChi: The Jacobian of the elastic part of the micro-deformation w.r.t.
         *     the plastic part of the micro deformation.
         * :param variableMatrix &dElasticGradChidGradChi: The Jacobian of the elastic part of the gradient of the micro-deformation
         *     w.r.t. the gradient of the micro-deformation.
         * :param variableMatrix &dElasticGradChidPlasticGradChi: The Jacobian of the elastic part of the gradient of the 
         *     micro-deformation w.r.t. the plastic part of the gradient of the micro-deformation.
         * :param variableMatrix &dElasticGradChidPlasticF: The Jacobian of the elastic part of the gradient of the micro-deformation
         *     w.r.t. the plastic deformation gradient.
         * :param variableMatrix &dElasticGradChidChi: The Jacobian of the elastic part of the gradient of the micro-deformation 
         *     w.r.t. the micro deformation.
         * :param variableMatrix &dElasticGradChidPlasticChi: The Jacobian of the elastic part of the gradient of the micro-deformation
         *     w.r.t. the plastic part of the micro-deformation.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Compute the required deformation measures
        variableVector inversePlasticDeformationGradient, inversePlasticMicroDeformation;

        errorOut error = computeElasticPartOfDeformation( deformationGradient, microDeformation, gradientMicroDeformation,
                                                          plasticDeformationGradient, plasticMicroDeformation, 
                                                          plasticGradientMicroDeformation, inversePlasticDeformationGradient,
                                                          inversePlasticMicroDeformation, elasticDeformationGradient,
                                                          elasticMicroDeformation, elasticGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "computeElasticPartOfDeformation (jacobian)",
                                             "Error in computation of the elastic part of the deformation measures" );
            result->addNext( error );
            return result;
        }

        //Assemble the Jacobians
        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        dElasticFdF = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dElasticFdPlasticF = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dElasticChidChi = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dElasticChidPlasticChi = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int Ib = 0; Ib < dim; Ib++ ){
                for ( unsigned int n = 0; n < dim; n++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        dElasticFdF[ dim * i + Ib ][ dim * n + N ] = eye[ dim * i + n ] * inversePlasticDeformationGradient[ dim * N + Ib ];
                        dElasticChidChi[ dim * i + Ib ][ dim * n + N ] = eye[ dim * i + n ] * inversePlasticMicroDeformation[ dim * N + Ib ];

                        for ( unsigned int I = 0; I < dim; I++ ){
                            dElasticFdPlasticF[ dim * i + Ib ][ dim * n + N ] -= deformationGradient[ dim * i + I ]
                                                                               * inversePlasticDeformationGradient[ dim * I + n ]
                                                                               * inversePlasticDeformationGradient[ dim * N + Ib ];
                            dElasticChidPlasticChi[ dim * i + Ib ][ dim * n + N ] -= microDeformation[ dim * i + I ]
                                                                                   * inversePlasticMicroDeformation[ dim * I + n ]
                                                                                   * inversePlasticMicroDeformation[ dim * N + Ib ];
                        }
                    }
                }
            }
        }

        dElasticGradChidPlasticF = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dElasticGradChidChi = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dElasticGradChidPlasticChi = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );

        dElasticGradChidGradChi = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dElasticGradChidPlasticGradChi = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int k = 0; k < dim; k++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    for ( unsigned int n = 0; n < dim; n++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            for ( unsigned int K = 0; K < dim; K++ ){

                                dElasticGradChidGradChi[ dim * dim * k + dim * Kb + Lb ][ dim * dim * n + dim * N + K ]
                                    = eye[ dim * k + n ] * inversePlasticDeformationGradient[ dim * K + Lb ]
                                    * inversePlasticMicroDeformation[ dim * N + Kb ];
                                
                                dElasticGradChidPlasticGradChi[ dim * dim * k + dim * Kb + Lb ][ dim * dim * n + dim * N + K ]
                                    = -elasticMicroDeformation[ dim * k + n ] * eye[ dim * Lb + K ]
                                    * inversePlasticMicroDeformation[ dim * N + Kb ];

                                for ( unsigned int Ab = 0; Ab < dim; Ab++ ){
                                    dElasticGradChidPlasticF[ dim * dim * k + dim * Kb + Lb ][ dim * n + N ]
                                        -= gradientMicroDeformation[ dim * dim * k + dim * K + Ab ]
                                         * inversePlasticDeformationGradient[ dim * Ab + n ]
                                         * inversePlasticDeformationGradient[ dim * N + Lb ]
                                         * inversePlasticMicroDeformation[ dim * K + Kb ];

                                    dElasticGradChidChi[ dim * dim * k + dim * Kb + Lb ][ dim * n + N ]
                                        -= dElasticChidChi[ dim * k + Ab ][ dim * n + N ]
                                        *  plasticGradientMicroDeformation[ dim * dim * Ab + dim * K + Lb ]
                                        *  inversePlasticMicroDeformation[ dim * K + Kb ];

                                    dElasticGradChidPlasticChi[ dim * dim * k + dim * Kb + Lb ][ dim * n + N ]
                                        -= dElasticChidPlasticChi[ dim * k + Ab ][ dim * n + N ]
                                        *  plasticGradientMicroDeformation[ dim * dim * Ab + dim * K + Lb ]
                                        *  inversePlasticMicroDeformation[ dim * K + Kb]
                                        +  ( gradientMicroDeformation[ dim * dim * k + dim * K + Ab ]
                                        *    inversePlasticDeformationGradient[ dim * Ab + Lb ]
                                        -    elasticMicroDeformation[ dim * k + Ab ]
                                        *    plasticGradientMicroDeformation[ dim * dim * Ab + dim * K + Lb ] )
                                        *    inversePlasticMicroDeformation[ dim * K + n ]
                                        *    inversePlasticMicroDeformation[ dim * N + Kb ] ;
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeElasticDeformationMeasures( const variableVector &elasticDeformationGradient,
                                                const variableVector &elasticMicroDeformation,
                                                const variableVector &elasticGradientMicroDeformation,
                                                variableVector &elasticRightCauchyGreen,
                                                variableVector &elasticMicroRightCauchyGreen,
                                                variableVector &elasticPsi, variableVector &elasticGamma ){
        /*!
         * Compute the elastic deformation measures required for the evolution of plasticity.
         *
         * :param const variableVector &elasticDeformationGradient: The elastic part of the deformation gradient.
         * :param const variableVector &elasticMicroDeformation: The elastic part of the micro-deformation.
         * :param const variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the
         *     micro-deformation.
         * :param variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation measure.
         * :param variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation 
         *     measure.
         * :param variableVector &elasticPsi: The elastic part of the micro-deformation measure Psi.
         * :param variableVector &elasticGamma: The elastic part of the higher order deformation measure.
         */

        errorOut error = tardigradeConstitutiveTools::computeRightCauchyGreen( elasticDeformationGradient, elasticRightCauchyGreen );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = tardigradeConstitutiveTools::computeRightCauchyGreen( elasticMicroDeformation, elasticMicroRightCauchyGreen );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic micro right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = tardigradeMicromorphicTools::computePsi( elasticDeformationGradient, elasticMicroDeformation, elasticPsi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic micro-deformation metric Psi" );
            result->addNext( error );
            return result;
        }

        error = tardigradeMicromorphicTools::computeGamma( elasticDeformationGradient, elasticGradientMicroDeformation, elasticGamma );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic higher order deformation metric Gamma" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computeElasticDeformationMeasures( const variableVector &elasticDeformationGradient,
                                                const variableVector &elasticMicroDeformation,
                                                const variableVector &elasticGradientMicroDeformation,
                                                variableVector &elasticRightCauchyGreen,
                                                variableVector &elasticMicroRightCauchyGreen,
                                                variableVector &elasticPsi, variableVector &elasticGamma,
                                                variableMatrix &dElasticRCGdElasticF, variableMatrix &dElasticMicroRCGdElasticChi,
                                                variableMatrix &dElasticPsidElasticF, variableMatrix &dElasticPsidElasticChi,
                                                variableMatrix &dElasticGammadElasticF, variableMatrix &dElasticGammadElasticGradChi ){
        /*!
         * Compute the elastic deformation measures required for the evolution of plasticity.
         *
         * :param const variableVector &elasticDeformationGradient: The elastic part of the deformation gradient.
         * :param const variableVector &elasticMicroDeformation: The elastic part of the micro-deformation.
         * :param const variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the
         *     micro-deformation.
         * :param variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation measure.
         * :param variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation 
         *     measure.
         * :param variableVector &elasticPsi: The elastic part of the micro-deformation measure Psi.
         * :param variableVector &elasticGamma: The elastic part of the higher order deformation measure.
         * :param variableMatrix &dElasticRCGdElasticF: The Jacobian of the right Cauchy-Green deformation measure
         *     w.r.t. the elastic deformation gradient.
         * :param variableMatrix &dElasticMicroRCGdElasticChi: The Jacobian of the micro right Cauchy-Green deformation 
         *     measure w.r.t. the elastic micro deformation
         * :param variableMatrix &dElasticPsidElasticF: The Jacobian of the elastic micro-deformation measure Psi
         *     w.r.t. the elastic deformation gradient.
         * :param variableMatrix &dElasticPsidElasticF: The Jacobian of the elastic micro-deformation measure Psi
         *     w.r.t. the elastic micro-deformation.
         * :param variableMatrix &dElasticGammadElasticF: The Jacobian of the higher order deformation measure Gamma w.r.t.
         *     the elastic deformation gradient.
         * :param variableMatrix &dElasticGammadElasticGradChi: The Jacobian of the higher order deformation measure Gamma
         *     w.r.t. the elastic part of the gradient of the micro-deformation.
         */

        errorOut error = tardigradeConstitutiveTools::computeRightCauchyGreen( elasticDeformationGradient, elasticRightCauchyGreen,
                                                                     dElasticRCGdElasticF );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = tardigradeConstitutiveTools::computeRightCauchyGreen( elasticMicroDeformation, elasticMicroRightCauchyGreen,
                                                            dElasticMicroRCGdElasticChi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic micro right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = tardigradeMicromorphicTools::computePsi( elasticDeformationGradient, elasticMicroDeformation, elasticPsi,
                                               dElasticPsidElasticF, dElasticPsidElasticChi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic micro-deformation metric Psi" );
            result->addNext( error );
            return result;
        }

        error = tardigradeMicromorphicTools::computeGamma( elasticDeformationGradient, elasticGradientMicroDeformation, elasticGamma,
                                                 dElasticGammadElasticF, dElasticGammadElasticGradChi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic higher order deformation metric Gamma" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient ){
        /*!
         * Compute the plastic macro velocity gradient in the intermediate configuration.
         *
         * \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &inverseElasticRightCauchyGreen: The inverse of the elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( inverseElasticRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticMacroVelocityGradient",
                                  "The inverse elastic right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( macroFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticMacroVelocityGradient",
                                  "The macro flow direction tensor must be 3D" );
        }

        if ( microFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticMacroVelocityGradient",
                                  "The micro flow direction tensor must be 3D" );
        }

        //Compute the macro-scale velocity gradient
        plasticMacroVelocityGradient = variableVector( dim * dim, 0 );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    plasticMacroVelocityGradient[ dim * Bb + Kb ]
                        += inverseElasticRightCauchyGreen[ dim * Bb + Lb ]
                         * ( macroGamma * macroFlowDirection[ dim * Kb + Lb ]
                         +   microGamma * microFlowDirection[ dim * Kb + Lb ] );
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma ){
        /*!
         * Compute the plastic macro velocity gradient in the intermediate configuration.
         *
         * \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &inverseElasticRightCauchyGreen: The inverse of the elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &dPlasticMacroLdMacroGamma: The Jacobian of the plastic velocity gradient w.r.t. the 
         *     macro plastic multiplier.
         * :param variableVector &dPlasticMacroLdMicroGamma: The Jacobian of the plastic velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMacroVelocityGradient (plastic multiplier jacobian)",
                                             "Error in computation of the plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        dPlasticMacroLdMacroGamma = variableVector( dim * dim, 0 );
        dPlasticMacroLdMicroGamma = variableVector( dim * dim, 0 );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    dPlasticMacroLdMacroGamma[ dim * Bb + Kb ] += inverseElasticRightCauchyGreen[ dim * Bb + Lb ]
                                                                * macroFlowDirection[ dim * Kb + Lb ];

                    dPlasticMacroLdMicroGamma[ dim * Bb + Kb ] += inverseElasticRightCauchyGreen[ dim * Bb + Lb ]
                                                                * microFlowDirection[ dim * Kb + Lb ];
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma, 
                                                  variableMatrix &dPlasticMacroLdElasticRCG, 
                                                  variableMatrix &dPlasticMacroLdMacroFlowDirection, 
                                                  variableMatrix &dPlasticMacroLdMicroFlowDirection ){
        /*!
         * Compute the plastic macro velocity gradient in the intermediate configuration.
         *
         * \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &inverseElasticRightCauchyGreen: The inverse of the elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &dPlasticMacroLdMacroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     macro plastic multiplier.
         * :param variableVector &dPlasticMacroLdMicroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableMatrix &dPlasticMacroLdElasticRCG: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     elastic right Cauchy-Green deformation tensor.
         * :param variableMatrix &dPlasticMacroLdMacroFlowDirection: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     macro flow direction tensor.
         * :param variableMatrix &dPlasticMacroLdMicroFlowDirection: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     micro flow direction tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient,
                                                              dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMacroVelocityGradient (full jacobian)",
                                             "Error in computation of the plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        dPlasticMacroLdElasticRCG = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMacroLdMacroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMacroLdMicroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Ob = 0; Ob < dim; Ob++ ){
                    for ( unsigned int Pb = 0; Pb < dim; Pb++ ){
                        dPlasticMacroLdElasticRCG[ dim * Bb + Kb ][ dim * Ob + Pb ]
                            -= inverseElasticRightCauchyGreen[ dim * Bb + Ob ]
                             * plasticMacroVelocityGradient[ dim * Pb + Kb ];
                        
                        dPlasticMacroLdMacroFlowDirection[ dim * Bb + Kb ][ dim * Ob + Pb ]
                            += macroGamma * inverseElasticRightCauchyGreen[ dim * Bb + Pb ] * eye[ dim * Kb + Ob ];

                        dPlasticMacroLdMicroFlowDirection[ dim * Bb + Kb ][ dim * Ob + Pb ]
                            += microGamma * inverseElasticRightCauchyGreen[ dim * Bb + Pb ] * eye[ dim * Kb + Ob ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient ){
        /*!
         * Compute the plastic micro velocity gradient
         *
         *  \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
         *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
         *
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
         * :param const variableVector &microFlowDirection: The micro plastic flow direction.
         * :param variableVector &plasticMicroVelocityGradient: The plastic micro velocity gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticMicroRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient",
                                  "The elastic micro right Cauchy-Green deformation tensor is not 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient",
                                  "The elastic micro deformation tensor Psi is not 3D" );
        }

        if ( inverseElasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient",
                                  "The inverse of the elastic micro deformation tensor Psi is not 3D" );
        }

        if ( microFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient",
                                  "The micro flow direction of the elastic micro plastic flow direction is not 3D" );
        }

        plasticMicroVelocityGradient = variableVector( dim * dim, 0 );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    for ( unsigned int Nb = 0; Nb < dim; Nb++ ){
                        for ( unsigned int Eb = 0; Eb < dim; Eb++ ){
                            plasticMicroVelocityGradient[ dim * Bb + Kb ]
                                += microGamma
                                 * inverseElasticPsi[ dim * Bb + Lb ]
                                 * microFlowDirection[ dim * Eb + Lb ]
                                 * inverseElasticPsi[ dim * Eb + Nb ]
                                 * elasticMicroRightCauchyGreen[ dim * Nb + Kb ];
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma ){
        /*!
         * Compute the plastic micro velocity gradient
         *
         *  \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
         *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
         *
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
         * :param const variableVector &microFlowDirection: The micro plastic flow direction.
         * :param variableVector &plasticMicroVelocityGradient: The plastic micro velocity gradient.
         * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro plastic multiplier.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticMicroRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient (jacobian)",
                                  "The elastic micro right Cauchy-Green deformation tensor is not 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient (jacobian)",
                                  "The elastic micro deformation tensor Psi is not 3D" );
        }

        if ( inverseElasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient (jacobian)",
                                  "The inverse of the elastic micro deformation tensor Psi is not 3D" );
        }

        if ( microFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient (jacobian)",
                                  "The micro flow direction of the elastic micro plastic flow direction is not 3D" );
        }

        plasticMicroVelocityGradient = variableVector( dim * dim, 0 );
        dPlasticMicroLdMicroGamma = variableVector( dim * dim, 0 );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    for ( unsigned int Nb = 0; Nb < dim; Nb++ ){
                        for ( unsigned int Eb = 0; Eb < dim; Eb++ ){
                            dPlasticMicroLdMicroGamma[ dim * Bb + Kb ]
                                += inverseElasticPsi[ dim * Bb + Lb ]
                                 * microFlowDirection[ dim * Eb + Lb ]
                                 * inverseElasticPsi[ dim * Eb + Nb ]
                                 * elasticMicroRightCauchyGreen[ dim * Nb + Kb ];
                        }
                    }
                }
                plasticMicroVelocityGradient[ dim * Bb + Kb ] = microGamma * dPlasticMicroLdMicroGamma[ dim * Bb + Kb ];
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma,
                                                  variableMatrix &dPlasticMicroLdElasticMicroRCG,
                                                  variableMatrix &dPlasticMicroLdElasticPsi,
                                                  variableMatrix &dPlasticMicroLdMicroFlowDirection ){
        /*!
         * Compute the plastic micro velocity gradient
         *
         *  \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
         *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
         *
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
         * :param const variableVector &microFlowDirection: The micro plastic flow direction.
         * :param variableVector &plasticMicroVelocityGradient: The plastic micro velocity gradient.
         * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro plastic multiplier.
         * :param variableMatrix &dPlasticMicroLdElasticMicroRCG: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro right Cauchy-Green deformation tensor.
         * :param variableMatrix &dPlasticMicroLdElasticPsi: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro deformation measure Psi.
         * :param variableMatrix &dPlasticMicroLdMicroFlowDirection: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro flow direction.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen,
                                                              elasticPsi, inverseElasticPsi, microFlowDirection,
                                                              plasticMicroVelocityGradient, dPlasticMicroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMicroVelocityGradient (full jacobian)",
                                             "Error in computation of the plastic micro velocity gradient" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        //Assemble the Jacobians
        dPlasticMicroLdElasticMicroRCG = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroLdElasticPsi = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroLdMicroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Ob = 0; Ob < dim; Ob++ ){
                    for ( unsigned int Pb = 0; Pb < dim; Pb++ ){
                        dPlasticMicroLdElasticPsi[ dim * Bb + Kb ][ dim * Ob + Pb ]
                            -= inverseElasticPsi[ dim * Bb + Ob ] * plasticMicroVelocityGradient[ dim * Pb + Kb ];

                        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                            dPlasticMicroLdMicroFlowDirection[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                += microGamma * inverseElasticPsi[ dim * Bb + Pb ]
                                 * inverseElasticPsi[ dim * Ob + Lb ]
                                 * elasticMicroRightCauchyGreen[ dim * Lb + Kb ];

                            for ( unsigned int Eb = 0; Eb < dim; Eb++ ){
                                dPlasticMicroLdElasticMicroRCG[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                    += microGamma * inverseElasticPsi[ dim * Bb + Lb ] * microFlowDirection[ dim * Eb + Lb ]
                                     * inverseElasticPsi[ dim * Eb + Ob ] * eye[ dim * Kb + Pb ];
                                for ( unsigned int Nb = 0; Nb < dim; Nb++ ){
                                    dPlasticMicroLdElasticPsi[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                        -= microGamma * inverseElasticPsi[ dim * Bb + Lb ] * microFlowDirection[ dim * Eb + Lb ]
                                         * inverseElasticPsi[ dim * Eb + Ob ] * inverseElasticPsi[ dim * Pb + Nb ]
                                         * elasticMicroRightCauchyGreen[ dim * Nb + Kb ];
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableVector &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         */

        variableVector skewTerm;
        return computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                            elasticGamma, microGradientFlowDirection,
                                                            plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient,
                                                            skewTerm );

    }
    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableVector &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         * :param variableVector &skewTerm: The skew term ( times 2 ) from the higher order computation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( microGradientGamma.size() != dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The micro gradient plastic multiplier must have a length of 3" );
        }

        if ( elasticPsi.size() != dim  * dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The elastic micro deformation measure Psi must be 3D" );
        }

        if ( inverseElasticPsi.size() != dim  * dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The inverse elastic micro deformation measure Psi must be 3D" );
        }

        if ( elasticGamma.size() != dim * dim * dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The elastic higher order deformation measure Gamma must be 3D" );
        }

        if ( microGradientFlowDirection.size() != dim * dim * dim * dim ){
            return new errorNode( "computePlasticMacroGradientVelocityGradient",
                                  "The micro gradient flow direction must be 3D" );
        }

        if ( plasticMicroVelocityGradient.size() != dim * dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The plastic micro velocity gradient must be 3D" );
        }

        //Assemble the 'skew' term
        skewTerm = variableVector( dim * dim * dim, 0 );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int Mb = 0; Mb < dim; Mb++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Cb = 0; Cb < dim; Cb++ ){
                        for ( unsigned int Fb = 0; Fb < dim; Fb++ ){
                            skewTerm[ dim * dim * Db + dim * Mb + Kb ]
                                += plasticMicroVelocityGradient[ dim * Db + Cb ]
                                 * inverseElasticPsi[ dim * Cb + Fb ]
                                 * elasticGamma[ dim * dim * Fb + dim * Mb + Kb ]
                                 - plasticMicroVelocityGradient[ dim * Cb + Mb ]
                                 * inverseElasticPsi[ dim * Db + Fb ]
                                 * elasticGamma[ dim * dim * Fb + dim * Cb + Kb ];
                        }
                    }
                }
            }
        }

        plasticMicroGradientVelocityGradient = variableVector( dim * dim * dim, 0 );

        for ( unsigned int Nb = 0; Nb < dim; Nb++ ){
            for ( unsigned int Mb = 0; Mb < dim; Mb++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                        for ( unsigned int Ib = 0; Ib < dim; Ib++ ){
                            plasticMicroGradientVelocityGradient[ dim * dim * Nb + dim * Mb + Kb ]
                                += inverseElasticPsi[ dim * Nb + Lb ]
                                 * ( microGradientGamma[ Ib ] * microGradientFlowDirection[ dim * dim * dim * Ib + dim * dim * Kb + dim * Lb + Mb ]
                                 +   elasticPsi[ dim * Lb + Ib ] * skewTerm[ dim * dim * Ib + dim * Mb + Kb ] );
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableVector &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient gamma.
         * :param variableMatrix &dPlasticMicroGradientLdPlasticMicroL: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the platic micro velocity gradient.
         */

        variableVector skewTerm;
        return computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi, elasticGamma,
                                                            microGradientFlowDirection, plasticMicroVelocityGradient,
                                                            plasticMicroGradientVelocityGradient, skewTerm,
                                                            dPlasticMicroGradientLdMicroGradientGamma,
                                                            dPlasticMicroGradientLdPlasticMicroL );
    }

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableVector &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         * :param variableVector &skewTerm: Two times the skew term.
         * :param variableMatrix &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient gamma.
         * :param variableMatrix &dPlasticMicroGradientLdPlasticMicroL: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the platic micro velocity gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                                      elasticGamma, microGradientFlowDirection,
                                                                      plasticMicroVelocityGradient, 
                                                                      plasticMicroGradientVelocityGradient, skewTerm );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMicroGradientVelocityGradient (jacobian)",
                                             "Error in computation of the plasticMicroVelocityGradient" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        dPlasticMicroGradientLdPlasticMicroL = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroGradientLdMicroGradientGamma = variableMatrix( dim * dim * dim, variableVector( dim, 0 ) );

        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
            for ( unsigned int Mb = 0; Mb < dim; Mb++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Ob = 0; Ob < dim; Ob++ ){
                        for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                            dPlasticMicroGradientLdMicroGradientGamma[ dim * dim * Lb + dim * Mb + Kb ][ Ob ]
                                += inverseElasticPsi[ dim * Lb + Pb ]
                                 * microGradientFlowDirection[ dim * dim * dim * Ob + dim * dim * Kb + dim * Pb + Mb ];

                            for ( unsigned int Qb = 0; Qb < dim; Qb++ ){
                                dPlasticMicroGradientLdPlasticMicroL[ dim * dim * Lb + dim * Mb + Kb ][ dim * Ob + Pb ]
                                    += eye[ dim * Lb + Ob ] * inverseElasticPsi[ dim * Pb + Qb ]
                                     * elasticGamma[ dim * dim * Qb + dim * Mb + Kb ]
                                     - eye[ dim * Mb + Pb ] * inverseElasticPsi[ dim * Lb + Qb ]
                                     * elasticGamma[ dim * dim * Qb + dim * Ob + Kb ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL,
                                                          variableMatrix &dPlasticMicroGradientLdElasticPsi,
                                                          variableMatrix &dPlasticMicroGradientLdElasticGamma,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableVector &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient gamma.
         * :param variableMatrix &dPlasticMicroGradientLdPlasticMicroL: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the platic micro velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientLdElasticPsi: The Jacobian of the plastic micro gradient
         *     velocity gradient w.r.t. the elastic micro deformation tensor Psi.
         * :param variableMatrix &dPlasticMicroGradientLdElasticGamma: The Jacobian of the plastic micro gradient
         *     velocity gradient w.r.t. the elastic higher ordrer deformation tensor Gamma.
         * :param variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection: The Jacobian of the plastic micro gradient
         *     velocity gradient w.r.t. the micro gradient flow direction.
         */

        //Assume 3D
        unsigned int dim = 3;

        variableVector skewTerm;
        errorOut error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                                      elasticGamma, microGradientFlowDirection,
                                                                      plasticMicroVelocityGradient, 
                                                                      plasticMicroGradientVelocityGradient, skewTerm,
                                                                      dPlasticMicroGradientLdMicroGradientGamma,
                                                                      dPlasticMicroGradientLdPlasticMicroL );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMicroGradientVelocityGradient (jacobian)",
                                             "Error in computation of the plasticMicroVelocityGradient" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        dPlasticMicroGradientLdElasticPsi = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroGradientLdElasticGamma = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dPlasticMicroGradientLdMicroGradientFlowDirection = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim * dim, 0 ) );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int Mb = 0; Mb < dim; Mb++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Ob = 0; Ob < dim; Ob++ ){
                        for ( unsigned int Pb = 0; Pb < dim; Pb++ ){
                            dPlasticMicroGradientLdElasticPsi[ dim * dim * Db + dim * Mb + Kb ][ dim * Ob + Pb ]
                                += -inverseElasticPsi[ dim * Db + Ob ]
                                 * plasticMicroGradientVelocityGradient[ dim * dim * Pb + dim * Mb + Kb ]
                                 + inverseElasticPsi[ dim * Db + Ob ] * skewTerm[ dim * dim * Pb + dim * Mb + Kb ];

                            for ( unsigned int Qb = 0; Qb < dim; Qb++ ){
                                dPlasticMicroGradientLdElasticGamma[ dim * dim * Db + dim * Mb + Kb ][ dim * dim * Ob + dim * Pb + Qb ]
                                    -= plasticMicroVelocityGradient[ dim * Pb + Mb ]
                                     * inverseElasticPsi[ dim * Db + Ob ] * eye[ dim * Kb + Qb ]; 

                                for ( unsigned int Rb = 0; Rb < dim; Rb++ ){
                                    dPlasticMicroGradientLdElasticPsi[ dim * dim * Db + dim * Mb + Kb ][ dim * Ob + Pb ]
                                        -= plasticMicroVelocityGradient[ dim * Db + Qb ]
                                         * inverseElasticPsi[ dim * Qb + Ob ] * inverseElasticPsi[ dim * Pb + Rb ]
                                         * elasticGamma[ dim * dim * Rb + dim * Mb + Kb ]
                                         - plasticMicroVelocityGradient[ dim * Qb + Mb ]
                                         * inverseElasticPsi[ dim * Db + Ob ] * inverseElasticPsi[ dim * Pb + Rb ]
                                         * elasticGamma[ dim * dim * Rb + dim * Qb + Kb ];
                                    
                                    dPlasticMicroGradientLdElasticGamma[ dim * dim * Db + dim * Mb + Kb ][ dim * dim * Ob + dim * Pb + Qb ]
                                        += plasticMicroVelocityGradient[ dim * Db + Rb ] * inverseElasticPsi[ dim * Rb + Ob ]
                                         * eye[ dim * Mb + Pb ] * eye[ dim * Kb + Qb ];

                                    dPlasticMicroGradientLdMicroGradientFlowDirection[ dim * dim * Db + dim * Mb + Kb ][ dim * dim * dim * Ob + dim * dim * Pb + dim * Qb + Rb ]
                                        += inverseElasticPsi[ dim * Db + Qb ] * microGradientGamma[ Ob ]
                                         * eye[ dim * Kb + Pb ] * eye[ dim * Mb + Rb ];
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma,
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableVector &microGradientFlowDirection,
                                              variableVector &plasticMacroVelocityGradient, variableVector &plasticMicroVelocityGradient,
                                              variableVector &plasticMicroGradientVelocityGradient ){
        /*!
         * Compute the plastic velocity gradients in the intermediate configuration.
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation metric Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation metric Gamma.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param const variableVector &microGradientFlowDirection: The flow direction of the micro gradient plasticity.
         *     Note: This is a matrix because it is computed as the gradient of the flow potential which is a vector.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &plasticMicroVelocityGradient: The plastic velocity gradient for the micro plastic deformation.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic velocity gradient for the micro gradient 
         *     plastic deformation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The elastic right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The elastic micro deformation metric Psi must be 3D" );
        }

        //Compute the required inverses of the deformation metrics
        variableVector inverseElasticRightCauchyGreen = tardigradeVectorTools::inverse( elasticRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = tardigradeVectorTools::inverse( elasticPsi, dim, dim );

        //Compute the macro velocity gradient
        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients",
                                             "Error in computation of plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro velocity gradient
        error = computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen, elasticPsi, inverseElasticPsi,
                                                     microFlowDirection, plasticMicroVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients",
                                             "Error in computation of plastic micro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro gradient velocity gradient
        error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                             elasticGamma, microGradientFlowDirection, 
                                                             plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients",
                                             "Error in computation of plastic micro gradient velocity gradient" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma,
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableVector &microGradientFlowDirection,
                                              variableVector &plasticMacroVelocityGradient, variableVector &plasticMicroVelocityGradient,
                                              variableVector &plasticMicroGradientVelocityGradient,
                                              variableVector &dPlasticMacroLdMacroGamma,
                                              variableVector &dPlasticMacroLdMicroGamma, variableVector &dPlasticMicroLdMicroGamma,
                                              variableVector &dPlasticMicroGradientLdMicroGamma,
                                              variableMatrix &dPlasticMicroGradientLdMicroGradientGamma ){
        /*!
         * Compute the plastic velocity gradients in the intermediate configuration.
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation metric Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation metric Gamma.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param const variableVector &microGradientFlowDirection: The flow direction of the micro gradient plasticity.
         *     Note: This is a matrix because it is computed as the gradient of the flow potential which is a vector.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &plasticMicroVelocityGradient: The plastic velocity gradient for the micro plastic deformation.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic velocity gradient for the micro gradient 
         *     plastic deformation.
         * :param variableVector &dPlasticMacroLdMacroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     macro plastic multiplier.
         * :param variableVector &dPlasticMacroLdMicroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableVector &dPlasticMicroGradientLdMicroGamma: The Jacobian of the plastic micro gradient velocity
         *     gradient w.r.t. the micro plastic multiplier.
         * :param variableVector &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient plastic multiplier.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients (jacobian)",
                                  "The elastic right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients (jacobian)",
                                  "The elastic micro deformation metric Psi must be 3D" );
        }

        //Compute the required inverses of the deformation metrics
        variableVector inverseElasticRightCauchyGreen = tardigradeVectorTools::inverse( elasticRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = tardigradeVectorTools::inverse( elasticPsi, dim, dim );

        //Compute the plastic macro velocity gradient
        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient,
                                                              dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (jacobian)",
                                             "Error in computation of plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro velocity gradient
        error = computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen, elasticPsi, inverseElasticPsi,
                                                     microFlowDirection, plasticMicroVelocityGradient,
                                                     dPlasticMicroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (jacobian)",
                                             "Error in computation of plastic micro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro gradient velocity gradient
        variableMatrix dPlasticMicroGradientLdPlasticMicroL;
        error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                             elasticGamma, microGradientFlowDirection, 
                                                             plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient,
                                                             dPlasticMicroGradientLdMicroGradientGamma,
                                                             dPlasticMicroGradientLdPlasticMicroL );

        dPlasticMicroGradientLdMicroGamma = tardigradeVectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                              dPlasticMicroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (jacobian)",
                                             "Error in computation of plastic micro gradient velocity gradient" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma,
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableVector &microGradientFlowDirection,
                                              variableVector &plasticMacroVelocityGradient, variableVector &plasticMicroVelocityGradient,
                                              variableVector &plasticMicroGradientVelocityGradient,
                                              variableVector &dPlasticMacroLdMacroGamma,
                                              variableVector &dPlasticMacroLdMicroGamma, variableVector &dPlasticMicroLdMicroGamma,
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
                                              variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection ){
        /*!
         * Compute the plastic velocity gradients in the intermediate configuration.
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation metric Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation metric Gamma.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param const variableVector &microGradientFlowDirection: The flow direction of the micro gradient plasticity.
         *     Note: This is a matrix because it is computed as the gradient of the flow potential which is a vector.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &plasticMicroVelocityGradient: The plastic velocity gradient for the micro plastic deformation.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic velocity gradient for the micro gradient 
         *     plastic deformation.
         * :param variableVector &dPlasticMacroLdMacroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     macro plastic multiplier.
         * :param variableVector &dPlasticMacroLdMicroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableVector &dPlasticMicroGradientLdMicroGamma: The Jacobian of the plastic micro gradient velocity
         *     gradient w.r.t. the micro plastic multiplier.
         * :param variableVector &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient plastic multiplier.
         * :param variableVector &dPlasticMacroLdElasticRCG: The Jacobian of the plastic macro velocity gradient w.r.t. 
         *     the elastic right Cauchy-Green deformation tensor.
         * :param variableVector &dPlasticMacroLdMacroFlowDirection: The Jacobian of the plastic macro velocity gradient
         *     w.r.t. the macro flow direction.
         * :param variableVector &dPlasticMacroLdMicroFlowDirection: The Jacobian of the plastic macro velocity gradient
         *     w.r.t. the micro flow direction.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                  "The elastic right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                  "The elastic micro deformation metric Psi must be 3D" );
        }

        //Compute the required inverses of the deformation metrics
        variableVector inverseElasticRightCauchyGreen = tardigradeVectorTools::inverse( elasticRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = tardigradeVectorTools::inverse( elasticPsi, dim, dim );

        //Compute the plastic macro velocity gradient
        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient,
                                                              dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma,
                                                              dPlasticMacroLdElasticRCG, dPlasticMacroLdMacroFlowDirection,
                                                              dPlasticMacroLdMicroFlowDirection );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                             "Error in computation of plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro velocity gradient
        error = computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen, elasticPsi, inverseElasticPsi,
                                                     microFlowDirection, plasticMicroVelocityGradient,
                                                     dPlasticMicroLdMicroGamma, dPlasticMicroLdElasticMicroRCG,
                                                     dPlasticMicroLdElasticPsi, dPlasticMicroLdMicroFlowDirection );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                             "Error in computation of plastic micro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro gradient velocity gradient
        variableMatrix dPlasticMicroGradientLdPlasticMicroL;
        error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                             elasticGamma, microGradientFlowDirection, 
                                                             plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient,
                                                             dPlasticMicroGradientLdMicroGradientGamma,
                                                             dPlasticMicroGradientLdPlasticMicroL,
                                                             dPlasticMicroGradientLdElasticPsi,
                                                             dPlasticMicroGradientLdElasticGamma,
                                                             dPlasticMicroGradientLdMicroGradientFlowDirection );

        dPlasticMicroGradientLdMicroGamma = tardigradeVectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                              dPlasticMicroLdMicroGamma );

        dPlasticMicroGradientLdElasticMicroRCG = tardigradeVectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                                   dPlasticMicroLdElasticMicroRCG );

        dPlasticMicroGradientLdElasticPsi += tardigradeVectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                               dPlasticMicroLdElasticPsi );

        dPlasticMicroGradientLdMicroFlowDirection = tardigradeVectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                                      dPlasticMicroLdMicroFlowDirection );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                             "Error in computation of plastic micro gradient velocity gradient" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

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
                                        const parameterType alpha ){
        /*!
         * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
         *
         * :param const variableType &Dt: The change in time.
         * :param const variableVector &currentPlasticMicroDeformation: The inverse of the current micro deformation.
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
         *     velocity gradient.
         * :param const variableVector &previousPlasticMicroDeformation: The plastic micro deformation 
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradient: The micro gradient deformation in the 
         *     intermediate configuation from the last converged increment.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
         *     velocity gradient from the last converged increment.
         * :param variableVector &currentPlasticMicroGradient: The current plastic micro gradient 
         *    deformation in the intermediate configuration.
         * :param parameterType alpha: The integration parameter.
         */

        variableMatrix LHS;
        return evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                          currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                          previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                          previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                          previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient, LHS,
                                          alpha );
    }

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
                                        const parameterType alpha ){
        /*!
         * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
         *
         * :param const variableType &Dt: The change in time.
         * :param const variableVector &currentPlasticMicroDeformation: The inverse of the current micro deformation.
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
         *     velocity gradient.
         * :param const variableVector &previousPlasticMicroDeformation: The the plastic micro deformation 
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradient: The micro gradient deformation in the 
         *     intermediate configuation from the last converged increment.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
         *     velocity gradient from the last converged increment.
         * :param variableVector &currentPlasticMicroGradient: The current plastic micro gradient 
         *    deformation in the intermediate configuration.
         * :param variableMatrix &LHS: The left-hand-side matrix.
         * :param parameterType alpha: The integration parameter.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( currentPlasticMicroDeformation.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The plastic micro-deformation must be 3D" );
        }

        if ( currentPlasticMacroVelocityGradient.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The plastic macro velocity gradient must be 3D" );
        }

        if ( currentPlasticMicroVelocityGradient.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The plastic micro velocity gradient must be 3D" );
        }

        if ( currentPlasticMicroGradientVelocityGradient.size() != dim * dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The plastic micro gradient velocity gradient must be 3D" );
        }

        if ( previousPlasticMicroDeformation.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic micro-deformation must be 3D" );
        }

        if ( previousPlasticMicroGradient.size() != dim * dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic micro gradient must be 3D" );
        }

        if ( previousPlasticMacroVelocityGradient.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic macro velocity gradient must be 3D" );
        }

        if ( previousPlasticMicroVelocityGradient.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic micro velocity gradient must be 3D" );
        }

        if ( previousPlasticMicroGradientVelocityGradient.size() != dim * dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic micro gradient velocity gradient must be 3D" );
        }

        //Compute the required identity terms
        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        //Assemble the A term ( forcing term ) and the fourth order A term
        variableVector DtAtilde( dim * dim * dim, 0 );
        variableVector previousFourthA( dim * dim * dim * dim, 0 );
        variableVector currentFourthA( dim * dim * dim * dim, 0 );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int B = 0; B < dim; B++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                        DtAtilde[ dim * dim * Db + dim * B + Kb ] += Dt 
                            * ( alpha * previousPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Lb + Kb ]
                                      *  previousPlasticMicroDeformation[ dim * Lb + B ]
                            + ( 1. - alpha ) * currentPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Lb + Kb ]
                                             * currentPlasticMicroDeformation[ dim * Lb + B ] );

                        previousFourthA[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Lb ]
                            = ( previousPlasticMicroVelocityGradient[ dim * Db + B ] * eye[ dim * Kb + Lb ]
                            -   previousPlasticMacroVelocityGradient[ dim * Lb + Kb ] * eye[ dim * Db + B ] );

                        currentFourthA[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Lb ]
                            = ( currentPlasticMicroVelocityGradient[ dim * Db + B ] * eye[ dim * Kb + Lb ]
                            -   currentPlasticMacroVelocityGradient[ dim * Lb + Kb ] * eye[ dim * Db + B ] );
                    }
                }
            }
        }

        //Assemble the right-hand side and left-hand side term
        variableVector RHS = DtAtilde;
        LHS = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int B = 0; B < dim; B++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                       for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
                          RHS[ dim * dim * Db + dim * B + Kb ]
                             += ( eye[ dim * Db + Bb ] * eye[ dim * Kb + Lb ] + Dt * alpha * previousFourthA[ dim * dim * dim * Db + dim * dim * Bb + dim * Kb + Lb ] )
                              * previousPlasticMicroGradient[ dim * dim * Bb + dim * B + Lb ];
                          for ( unsigned int Sb = 0; Sb < dim; Sb++ ){
                              LHS[ dim * dim * Db + dim * B + Kb ][ dim * dim * Lb + dim * Bb + Sb ]
                                  = ( eye[ dim * Db + Lb ] * eye[ dim * Kb + Sb ] - Dt * ( 1. - alpha ) * currentFourthA[ dim * dim * dim * Db + dim * dim * Lb + dim * Kb + Sb ] ) * eye[ dim * B + Bb ];
                          }
                       } 
                    }
                }
            }
        }

        //Solve for the current plastic micro gradient
        unsigned int rank;
        currentPlasticMicroGradient = tardigradeVectorTools::solveLinearSystem( LHS, RHS, rank );

        if ( rank != LHS.size() ){
            std::cout << "rank: " << rank << "\n";
//            const variableType &Dt,
//            const variableVector &currentPlasticMicroDeformation,
//            const variableVector &currentPlasticMacroVelocityGradient,
//            const variableVector &currentPlasticMicroVelocityGradient,
//            const variableVector &currentPlasticMicroGradientVelocityGradient,
//            const variableVector &previousPlasticMicroDeformation,
//            const variableVector &previousPlasticMicroGradient,
//            const variableVector &previousPlasticMacroVelocityGradient,
//            const variableVector &previousPlasticMicroVelocityGradient,
//            const variableVector &previousPlasticMicroGradientVelocityGradient,
//            variableVector &currentPlasticMicroGradient
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The left hand side matrix is not full rank" );
        }

        return NULL;
    }

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
                                        const parameterType alpha ){
        /*!
         * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
         *
         * :param const variableType &Dt: The change in time.
         * :param const variableVector &currentPlasticMicroDeformation: The inverse of the current micro deformation.
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
         *     velocity gradient.
         * :param const variableVector &previousPlasticMicroDeformation: The the plastic micro deformation 
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradient: The micro gradient deformation in the 
         *     intermediate configuation from the last converged increment.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
         *     velocity gradient from the last converged increment.
         * :param variableVector &currentPlasticMicroGradient: The current plastic micro gradient 
         *    deformation in the intermediate configuration.
         * :param variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroDeformation: The jacobian of the plastic 
         *     micro deformation w.r.t. the plastic micro deformation.
         * :param variableMatrix &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient: The jacobian of the plastic 
         *     micro deformation w.r.t. the plastic macro velocity gradient.
         * :param variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient: The jacobian of the plastic 
         *     micro deformation w.r.t. the plastic micro velocity gradient.
         * :param variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient: The jacobian of the plastic 
         *     micro deformation w.r.t. the plastic micro gradient velocity gradient.
         * :param parameterType alpha: The integration parameter.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Compute the required identity terms
        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        //Compute the new currentPlasticMicroGradient
        variableMatrix LHS;
        errorOut error = evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                                    currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                                    previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                                    previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                                    previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                                    LHS, alpha );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticMicroGradChi",
                                             "Error in the evolution of the plastic micro gradient of chi" );
            result->addNext( error );
            return result;
        }

        //Compute the negative partial derivatives w.r.t. currentFourthA and the current part of DtAtilde
        //We do this in vector form so that we can interface with Eigen easier
        variableVector negdRdCurrentDtAtilde( dim * dim * dim * dim * dim * dim, 0 );
        variableVector negdRdCurrentFourthA( dim * dim * dim * dim * dim * dim * dim, 0 );

        //Also assemble jacobians of the A terms
        variableMatrix dCurrentDTAtildedPlasticMicroDeformation( dim * dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dCurrentDTAtildedPlasticMicroGradientVelocityGradient( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        variableMatrix dCurrentFourthAdMacroVelocityGradient( dim * dim * dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dCurrentFourthAdMicroVelocityGradient( dim * dim * dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int B = 0; B < dim; B++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Rb = 0; Rb < dim; Rb++ ){
                        for ( unsigned int S = 0; S < dim; S++ ){
                            dCurrentDTAtildedPlasticMicroDeformation[ dim * dim * Db + dim * B + Kb ][ dim * Rb + S ]
                                += Dt * ( 1. - alpha ) * currentPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Rb + Kb ]
                                 * eye[ dim * B + S ];

                            for ( unsigned int Tb = 0; Tb < dim; Tb++ ){
                                negdRdCurrentDtAtilde[ dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * B + dim * dim * dim * Kb + dim * dim * Rb + dim * S + Tb ]
                                    += eye[ dim * Db + Rb ] * eye[ dim * B + S ] * eye[ dim * Kb + Tb ];

                                dCurrentDTAtildedPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * B + Kb ][ dim * dim * Rb + dim * S + Tb ]
                                    += Dt * ( 1. - alpha ) * eye[ dim * Db + Rb ] * eye[ dim * Kb + Tb ] * currentPlasticMicroDeformation[ dim * S + B ];

                                dCurrentFourthAdMacroVelocityGradient[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Rb ][ dim * S + Tb ]
                                    -= eye[ dim * Rb + S ] * eye[ dim * Kb + Tb ] * eye[ dim * Db + B ];

                                dCurrentFourthAdMicroVelocityGradient[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Rb ][ dim * S + Tb ]
                                    += eye[ dim * Db + S ] * eye[ dim * B + Tb ] * eye[ dim * Kb + Rb ];

                                for ( unsigned int Ub = 0; Ub < dim; Ub++ ){
                                    negdRdCurrentFourthA[ dim * dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * dim * B + dim * dim * dim * dim * Kb + dim * dim * dim * Rb + dim * dim * S + dim * Tb + Ub ]
                                        += Dt * ( 1. - alpha ) * eye[ dim * Db + Rb ] * eye[ dim * Kb + Tb ]
                                         * currentPlasticMicroGradient[ dim * dim * S + dim * B + Ub ];

                                }
                            }
                        }
                    }
                }
            }
        }

        //Solve for the Jacobians
        variableVector vecdCurrentPlasticMicroGradientdCurrentDTAtilde( dim * dim * dim * dim * dim * dim );
        variableVector vecdCurrentPlasticMicroGradientdCurrentFourthA( dim * dim * dim * dim * dim * dim * dim );

        variableVector floatLHS = tardigradeVectorTools::appendVectors( LHS );

        Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > LHSMat( floatLHS.data(), LHS.size(), LHS.size() );
        Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCDA( negdRdCurrentDtAtilde.data(), LHS.size(), dim * dim * dim );
        Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCFA( negdRdCurrentFourthA.data(), LHS.size(), dim * dim * dim * dim );

        Eigen::ColPivHouseholderQR< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > qrSolver( LHSMat );

        Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X1( vecdCurrentPlasticMicroGradientdCurrentDTAtilde.data(), LHS.size(), dim * dim * dim );
        Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X2( vecdCurrentPlasticMicroGradientdCurrentFourthA.data(), LHS.size(), dim * dim * dim * dim );

        X1 = qrSolver.solve( nDRDCDA );
        X2 = qrSolver.solve( nDRDCFA );

        variableMatrix dCurrentPlasticMicroGradientdCurrentDTAtilde = tardigradeVectorTools::inflate( vecdCurrentPlasticMicroGradientdCurrentDTAtilde, dim * dim * dim, dim * dim * dim );
        variableMatrix dCurrentPlasticMicroGradientdCurrentFourthA = tardigradeVectorTools::inflate( vecdCurrentPlasticMicroGradientdCurrentFourthA, dim * dim * dim, dim * dim * dim * dim );

        //Assemble the final terms of the deformation
        dCurrentPlasticMicroGradientdPlasticMicroDeformation = tardigradeVectorTools::dot( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                 dCurrentDTAtildedPlasticMicroDeformation );

        dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient = tardigradeVectorTools::dot( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                      dCurrentFourthAdMacroVelocityGradient );

        dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient = tardigradeVectorTools::dot( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                      dCurrentFourthAdMicroVelocityGradient );

        dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient = tardigradeVectorTools::dot( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                 dCurrentDTAtildedPlasticMicroGradientVelocityGradient );
        return NULL;
    }

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
                                       const parameterType alphaMacro,
                                       const parameterType alphaMicro,
                                       const parameterType alphaMicroGradient ){
        /*!
         * Evolve the plastic deformation
         *
         * :param const variableType &Dt: The timestep
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient 
         *     velocity gradient.
         * :param const variableVector &previousPlasticDeformationGradient: The plastic deformation gradient at the end of the last 
         *     converged timestep.
         * :param const variableVector &previousPlasticMicroDeformation: The plastic micro deformation at the end of the last converged 
         *     timestep.
         * :param const variableVector &previousPlasticMicroGradient: The plastic micro gradient at the end of the last converged 
         *     timestep.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient at the end of the 
         *     last converged timestep.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient at the end of the 
         *     last converged timestep.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient 
         *     at the end of the last converged timestep.
         * :param variableVector &currentPlasticDeformationGradient: The current value of the plastic deformation gradient.
         * :param variableVector &currentPlasticMicroDeformation: The current value of the plastic micro deformation.
         * :param variableVector &currentPlasticMicroGradient: The current value of the plastic micro gradient.
         * :param parameterType alphaMacro: The integration parameter for the macro plasticity. Defaults to 0.5.
         * :param parameterType alphaMicro: The integration parameter for the micro plasticity. Defaults to 0.5.
         * :param parameterType alphaMicroGradient: The integration parameter for the micro gradient plasticity. Defaults to 0.5.
         */

        errorOut error = tardigradeConstitutiveTools::evolveF( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
                                                     currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
                                                     alphaMacro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation",
                                             "Error in computation of the plastic macro deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = tardigradeConstitutiveTools::evolveF( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            alphaMicro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation",
                                             "Error in computation of the plastic micro deformation" );
            result->addNext( error );
            return result;
        }

        error = evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                           previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                           previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                           previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                           alphaMicroGradient );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation",
                                             "Error in computation of the plastic micro gradient" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

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
                                       const parameterType alphaMacro,
                                       const parameterType alphaMicro,
                                       const parameterType alphaMicroGradient ){
        /*!
         * Evolve the plastic deformation
         *
         * :param const variableType &Dt: The timestep
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient 
         *     velocity gradient.
         * :param const variableVector &previousPlasticDeformationGradient: The plastic deformation gradient at the end of the last 
         *     converged timestep.
         * :param const variableVector &previousPlasticMicroDeformation: The plastic micro deformation at the end of the last converged 
         *     timestep.
         * :param const variableVector &previousPlasticMicroGradient: The plastic micro gradient at the end of the last converged 
         *     timestep.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient at the end of the 
         *     last converged timestep.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient at the end of the 
         *     last converged timestep.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient 
         *     at the end of the last converged timestep.
         * :param variableVector &currentPlasticDeformationGradient: The current value of the plastic deformation gradient.
         * :param variableVector &currentPlasticMicroDeformation: The current value of the plastic micro deformation.
         * :param variableVector &currentPlasticMicroGradient: The current value of the plastic micro gradient.
         * :param variableMatrix &dPlasticFdPlasticMacroL: The Jacobian of the plastic deformation gradient w.r.t. the plastic 
         *     macro velocity gradient.
         * :param variableMatrix &dPlasticMicroDeformationdPlasticMicroL: The Jacobian of the plastic micro-deformation w.r.t. 
         *     the plastic micro velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientdPlasticMacroL: The Jacobian of the plastic micro gradient deformation 
         *     w.r.t. the plastic macro velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientdPlasticMicroL: The Jacobian of the plastic micro gradient deformation
         *     w.r.t. the plastic micro velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientdPlasticMicroGradientL: The Jacobian of the plastic micro gradient deformation
         *     w.r.t. the plastic micro gradient velocity gradient.
         * :param parameterType alphaMacro: The integration parameter for the macro plasticity. Defaults to 0.5.
         * :param parameterType alphaMicro: The integration parameter for the micro plasticity. Defaults to 0.5.
         * :param parameterType alphaMicroGradient: The integration parameter for the micro gradient plasticity. Defaults to 0.5.
         */

        errorOut error = tardigradeConstitutiveTools::evolveF( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
                                                     currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
                                                     dPlasticFdPlasticMacroL, alphaMacro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation (jacobian)",
                                             "Error in computation of the plastic macro deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = tardigradeConstitutiveTools::evolveF( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            dPlasticMicroDeformationdPlasticMicroL, alphaMicro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation (jacobian)",
                                             "Error in computation of the plastic micro deformation" );
            result->addNext( error );
            return result;
        }

        variableMatrix dPlasticMicroGradientdPlasticMicroDeformation;
        error = evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                           previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                           previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                           previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                           dPlasticMicroGradientdPlasticMicroDeformation,
                                           dPlasticMicroGradientdPlasticMacroL, dPlasticMicroGradientdPlasticMicroL,
                                           dPlasticMicroGradientdPlasticMicroGradientL, alphaMicroGradient );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation (jacobian)",
                                             "Error in computation of the plastic micro gradient" );
            result->addNext( error );
            return result;
        }

        dPlasticMicroGradientdPlasticMicroL += tardigradeVectorTools::dot( dPlasticMicroGradientdPlasticMicroDeformation,
                                                                 dPlasticMicroDeformationdPlasticMicroL );

        return NULL;
    }

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
                                         const parameterType alphaMacro,
                                         const parameterType alphaMicro,
                                         const parameterType alphaMicroGradient ){
        /*!
         * Evolve the strain-like state variables.
         *
         * :param const constantType &Dt: The timestep
         * :param const variableType &currentMacroGamma: The current macro-scale plastic multiplier.
         * :param const variableType &currentMicroGamma: The current micro-scale plastic multiplier.
         * :param const variableVector &currentMicroGradientGamma: The current micro-scale gradient plastic multiplier.
         * :param const variableType &currentdMacroGdMacroC: The current Jacobian of the macro plastic potential 
         *     w.r.t. the macro cohesion.
         * :param const variableType &currentdMicroGdMicroC: The current Jacobian of the micro plastic potential 
         *     w.r.t. the micro cohesion.
         * :param const variableVector &currentdMicroGradientGdMicroGradientC: The current Jacobian of the micro gradient 
         *     plastic potential w.r.t. the micro gradient cohesion.
         * :param variableType &previousMacroStrainISV: The previous value of the macro strain like ISV.
         * :param variableType &previousMicroStrainISV: The previous value of the micro strain like ISV.
         * :param variableVector &previousMicroGradientStrainISV: The previous value of the micro gradient strain like ISV.
         * :param const variableType &previousMacroGamma: The previous macro-scale plastic multiplier.
         * :param const variableType &previousMicroGamma: The previous micro-scale plastic multiplier.
         * :param const variableVector &previousMicroGradientGamma: The previous micro-scale gradient plastic multiplier.
         * :param const variableType &previousdMacroGdMacroC: The previous Jacobian of the macro plastic potential 
         *     w.r.t. the macro cohesion.
         * :param const variableType &previousdMicroGdMicroC: The previous Jacobian of the micro plastic potential 
         *     w.r.t. the micro cohesion.
         * :param const variableVector &previousdMicroGradientGdMicroGradientC: The previous Jacobian of the micro gradient 
         *     plastic potential w.r.t. the micro gradient cohesion.
         * :param variableType &currentMacroStrainISV: The new value of the macro strain ISV
         * :param variableType &currentMicroStrainISV: The new value of the micro strain ISV
         * :param variableVector &currentMicroGradientStrainISV: The new value of the micro gradient strain ISV.
         * :param const variableType alphaMacro: The integration parameter for the macro strain-like ISV
         * :param const variableType alphaMicro: The integration parameter for the micro strain-like ISV
         * :param const variableType alphaMicroGradient: The integration parameter for the micro Gradient strain-like ISV
         */

        if ( currentMicroGradientGamma.size() != 3 ){
            return new errorNode( "evolveStrainStateVariables",
                                  "The current micro gradient gamma must be 3 dimensional" );
        }

        if ( previousMicroGradientGamma.size() != 3 ){
            return new errorNode( "evolveStrainStateVariables",
                                  "The previous micro gradient gamma must be 3 dimensional" );
        }

        if ( currentMicroGradientGamma.size() != currentdMicroGradientGdMicroGradientC.size() ){
            return new errorNode( "evolveStrainStateVariables",
                                  "The current micro gradient gamma and the current derivative of the plastic potential function w.r.t. the micro gradient cohesion are not consistent" );
        }

        if ( previousMicroGradientGamma.size() != previousdMicroGradientGdMicroGradientC.size() ){
            return new errorNode( "evolveStrainStateVariables",
                                  "The previous micro gradient gamma and the previous derivative of the plastic potential function w.r.t. the micro gradient cohesion are not consistent" );
        }

        for ( unsigned int i = 0; i < currentdMicroGradientGdMicroGradientC.size(); i++ ){
            if ( currentdMicroGradientGdMicroGradientC[ i ].size() != previousMicroGradientStrainISV.size() ){
                return new errorNode( "evolveStrainStateVariables",
                                      "The current derivative of the plastic potential function w.r.t. the micro gradient cohesion is not a square matrix" ); 
            }
        }

        for ( unsigned int i = 0; i < previousdMicroGradientGdMicroGradientC.size(); i++ ){
            if ( previousdMicroGradientGdMicroGradientC[ i ].size() != previousMicroGradientStrainISV.size() ){
                return new errorNode( "evolveStrainStateVariables",
                                      "The previous derivative of the plastic potential function w.r.t. the micro gradient cohesion is not a square matrix" ); 
            }
        }
        
        if ( std::isnan( currentMacroGamma ) ){

            return new errorNode( __func__, "The current macro plastic multiplier (gamma) is nan" );

        }

        if ( std::isnan( currentMicroGamma ) ){

            return new errorNode( __func__, "The current micro plastic multiplier (gamma) is nan" );

        }

        for ( auto cMGG = currentMicroGradientGamma.begin( ); cMGG != currentMicroGradientGamma.end( ); cMGG++ ){

            if ( std::isnan( *cMGG ) ){
    
                return new errorNode( __func__, "The " + std::to_string( cMGG - currentMicroGradientGamma.begin( ) ) + "th index of the current micro gradient plastic multiplier (gamma) is nan" );
    
            }

        }

        if ( std::isnan( currentdMacroGdMacroC ) ){

            return new errorNode( __func__, "The current gradient of the macro plastic potential w.r.t. the macro cohesion is nan" );

        }

        if ( std::isnan( currentdMicroGdMicroC ) ){

            return new errorNode( __func__, "The current gradient of the micro plastic potential w.r.t. the micro cohesion is nan" );

        }

        for ( auto cMGG = currentdMicroGradientGdMicroGradientC.begin( ); cMGG != currentdMicroGradientGdMicroGradientC.end( ); cMGG++ ){

            for ( auto cMGG2 = cMGG->begin( ); cMGG2 != cMGG->end( ); cMGG2++ ){

                if ( std::isnan( *cMGG2 ) ){
        
                    return new errorNode( __func__, "The " + std::to_string( cMGG - currentdMicroGradientGdMicroGradientC.begin( ) ) + "," + std::to_string( cMGG2 - cMGG->begin( ) ) + "th index of the current micro gradient plastic potential w.r.t. the micro-gradient cohesion is nan" );
        
                }

            }

        }

        //Evolve the macro-scale internal strain-like state variable
        currentMacroStrainISV = previousMacroStrainISV + Dt * (         - alphaMacro * previousMacroGamma * previousdMacroGdMacroC
                                                                - ( 1 - alphaMacro ) *  currentMacroGamma * currentdMacroGdMacroC );

        //Evolve the micro-scale internal strain-like state variable
        currentMicroStrainISV = previousMicroStrainISV + Dt * (         - alphaMicro * previousMicroGamma * previousdMicroGdMicroC
                                                                - ( 1 - alphaMicro ) *  currentMicroGamma * currentdMicroGdMicroC );

        //Evolve the micro gradient internal strain-like state variables
        currentMicroGradientStrainISV = previousMicroGradientStrainISV
        + Dt * ( 
            - alphaMicroGradient * tardigradeVectorTools::Tdot( previousdMicroGradientGdMicroGradientC, previousMicroGradientGamma )
            - ( 1 - alphaMicroGradient ) * tardigradeVectorTools::Tdot( currentdMicroGradientGdMicroGradientC, currentMicroGradientGamma )
        );

        return NULL;
    }

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
                                         const parameterType alphaMacro,
                                         const parameterType alphaMicro,
                                         const parameterType alphaMicroGradient){
        /*!
         * Evolve the strain-like state variables.
         *
         * :param const constantType &Dt: The timestep
         * :param const variableType &currentMacroGamma: The current macro-scale plastic multiplier.
         * :param const variableType &currentMicroGamma: The current micro-scale plastic multiplier.
         * :param const variableVector &currentMicroGradientGamma: The current micro-scale gradient plastic multiplier.
         * :param const variableType &currentdMacroGdMacroC: The current Jacobian of the macro plastic potential 
         *     w.r.t. the macro cohesion.
         * :param const variableType &currentdMicroGdMicroC: The current Jacobian of the micro plastic potential 
         *     w.r.t. the micro cohesion.
         * :param const variableVector &currentdMicroGradientGdMicroGradientC: The current Jacobian of the micro gradient 
         *     plastic potential w.r.t. the micro gradient cohesion.
         * :param variableType &previousMacroStrainISV: The previous value of the macro strain like ISV.
         * :param variableType &previousMicroStrainISV: The previous value of the micro strain like ISV.
         * :param variableVector &previousMicroGradientStrainISV: The previous value of the micro gradient strain like ISV.
         * :param const variableType &previousMacroGamma: The previous macro-scale plastic multiplier.
         * :param const variableType &previousMicroGamma: The previous micro-scale plastic multiplier.
         * :param const variableVector &previousMicroGradientGamma: The previous micro-scale gradient plastic multiplier.
         * :param const variableType &previousdMacroGdMacroC: The previous Jacobian of the macro plastic potential 
         *     w.r.t. the macro cohesion.
         * :param const variableType &previousdMicroGdMicroC: The previous Jacobian of the micro plastic potential 
         *     w.r.t. the micro cohesion.
         * :param const variableVector &previousdMicroGradientGdMicroGradientC: The previous Jacobian of the micro gradient 
         *     plastic potential w.r.t. the micro gradient cohesion.
         * :param variableType &dCurrentMacroISVdCurrentMacroGamma: The Jacobian of the 
         * :param variableType &currentMacroStrainISV: The new value of the macro strain ISV
         * :param variableType &currentMicroStrainISV: The new value of the micro strain ISV
         * :param variableVector &currentMicroGradientStrainISV: The new value of the micro gradient strain ISV.
         * :param variableType &dCurrentMacroISVdCurrentMacroGamma: The Jacobian of the current macro strain ISV w.r.t. the 
         *     macro plastic multiplier.
         * :param variableType &dCurrentMacroISVddMacroGdMacroC: The Jacobian of the current macro strain ISV w.r.t. the 
         *     derivative of the macro plastic potential w.r.t. the macro cohesion.
         * :param variableType &dCurrentMicroISVdCurrentMicroGamma: The Jacobian of the current micro strain ISV w.r.t. the 
         *     micro plastic multiplier.
         * :param variableType &dCurrentMicroISVddMicroGdMicroC: The Jacobian of the current micro strain ISV w.r.t. the 
         *     derivative of the micro plastic potential w.r.t. the micro cohesion.
         * :param variableMatrix &dCurrentMicroGradientISVdCurrentMicroGradientGamma: The Jacobian of the current micro 
         *     strain ISV w.r.t. the micro plastic multiplier.
         * :param variableMatrix &dCurrentMicroGradientISVddMicroGradientGdMicroGradientC: The Jacobian of the current 
         *     micro gradient strain ISV w.r.t. the derivative of the micro gradient plastic potential w.r.t. the micro 
         *     gradient cohesion.
         * :param const variableType alphaMacro: The integration parameter for the macro strain-like ISV
         * :param const variableType alphaMicro: The integration parameter for the micro strain-like ISV
         * :param const variableType alphaMicroGradient: The integration parameter for the micro Gradient strain-like ISV
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma, currentMicroGradientGamma,
                                                     currentdMacroGdMacroC, currentdMicroGdMicroC,
                                                     currentdMicroGradientGdMicroGradientC,
                                                     previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
                                                     previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
                                                     previousdMacroGdMacroC, previousdMicroGdMicroC,
                                                     previousdMicroGradientGdMicroGradientC,
                                                     currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                                     alphaMacro, alphaMicro, alphaMicroGradient );

        if ( error ){
            errorOut result = new errorNode( "evolveStrainStateVariables (jacobian)",
                                             "Error in the evolution of the strain-like ISVs" );
            result->addNext( error );
            return result;
        }

        //Compute the Jacobians of the macro ISV evolution w.r.t. the current macro gamma and the derivative of the plastic 
        //potential function w.r.t. the macro cohesion
        dCurrentMacroISVdCurrentMacroGamma = - Dt * ( 1 - alphaMacro ) * currentdMacroGdMacroC;
        dCurrentMacroISVddMacroGdMacroC = - Dt * ( 1 - alphaMacro ) * currentMacroGamma;

        //Compute the Jacobians of the micro ISV evolution w.r.t. the current micro gamma and the derivative of the plastic 
        //potential function w.r.t. the micro cohesion
        dCurrentMicroISVdCurrentMicroGamma = - Dt * ( 1 - alphaMicro ) * currentdMicroGdMicroC;
        dCurrentMicroISVddMicroGdMicroC = - Dt * ( 1 - alphaMicro ) * currentMicroGamma;

        //Compute the Jacobians of the micro gradient ISV evolution w.r.t. the current micro gradient gamma and the derivative
        //of the plastic potential function w.r.t. the micro gradient cohesion.
        dCurrentMicroGradISVdCurrentMicroGradGamma = variableMatrix( dim, variableVector( dim, 0 ) );
        dCurrentMicroGradISVddMicroGradGdMicroGradC = variableMatrix( dim, variableVector( dim * dim, 0 ) );

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                dCurrentMicroGradISVdCurrentMicroGradGamma[ i ][ j ] -= Dt * ( 1 - alphaMicroGradient )
                                                                      * currentdMicroGradientGdMicroGradientC[ j ][ i ];
                for ( unsigned int k = 0; k < dim; k++ ){
                    dCurrentMicroGradISVddMicroGradGdMicroGradC[ i ][ dim * j + k ]
                        -= Dt * ( 1 - alphaMicroGradient ) * currentMicroGradientGamma[ j ] * eye[ dim * i + k ];
                }
            }
        }

        return NULL;
    }

    errorOut computeFlowDirections( const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                    const variableVector &referenceHigherOrderStress, const variableType &macroCohesion,
                                    const variableType &microCohesion, const variableVector &microGradientCohesion,
                                    const variableVector &elasticRightCauchyGreen, const parameterVector &macroFlowParameters,
                                    const parameterVector &microFlowParameters, const parameterVector &microGradientFlowParameters,
                                    variableVector &macroFlowDirection, variableVector &microFlowDirection,
                                    variableVector &microGradientFlowDirection, variableType &dGdMacroCohesion, 
                                    variableType &dGdMicroCohesion, variableMatrix &dGdMicroGradientCohesion ){
        /*!
         * Compute all of the flow directions.
         *
         * :param const variableVector &PK2Stress: The second Piola-Kirchhoff stress.
         * :param const variableVector &referenceMicroStress: The reference symmetric micro stress.
         * :param const variableVector &referenceHigherOrderStress: The reference higher order stress.
         * :param const variableType &macroCohesion: The macro cohesion value.
         * :param const variableType &microCohesion: The micro cohesion value.
         * :param const variableVector &microGradientCohesion: The micro gradient cohesion value.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right cauchy green deformation.
         * :param const parameterVector &macroFlowParameters: The macro plastic parameters.
         *     [ friction angle, beta ]
         * :param const parameterVector &microFlowParameters: The micro plastic parameters.
         *     [ friction angle, beta ]
         * :param const parameterVector &microGradientFlowParameters: The micro gradient plastic parameters.
         *     [ friction angle, beta ]
         * :param variableVector &macroFlowDirection: The flow direction for the macro scale plasticity.
         * :param variableVector &microFlowDirection: The flow direction for the micro scale plasticity.
         * :param variableVector &microGradientFlowDirection: The flow direction for the micro gradient 
         *     plasticity.
         * :param variableType &dGdMacroCohesion: The Jacobian of the macro plastic potential w.r.t. the 
         *     macro cohesion.
         * :param variableType &dGdMicroCohesion: The Jacobian of the micro plastic potential w.r.t. the 
         *     micro cohesion.
         * :param variableVector &dGdMicroGradientCohesion: The Jacobian of the micro gradient plastic 
         *     potential w.r.t. the micro gradient cohesion.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Error handling
        if ( macroFlowParameters.size() != 2 ){
            return new errorNode( "computeFlowDirections",
                                  "The number of macro flow parameters must be 2" );
        }

        if ( microFlowParameters.size() != 2 ){
            return new errorNode( "computeFlowDirections",
                                  "The number of micro flow parameters must be 2" );
        }

        if ( microGradientFlowParameters.size() != 2 ){
            return new errorNode( "computeFlowDirections",
                                  "The number of micro gradient flow parameters must be 2" );
        }

        //Set temporary variables
        variableVector tmpVec;
        variableMatrix tmpMat;

        parameterType macroFrictionAngle = macroFlowParameters[ 0 ];
        parameterType macroBeta          = macroFlowParameters[ 1 ];
        parameterType macroFlowPotential;

        errorOut error;

        error = computeSecondOrderDruckerPragerYieldEquation( PK2Stress, macroCohesion, elasticRightCauchyGreen,
                                                              macroFrictionAngle, macroBeta, macroFlowPotential, 
                                                              macroFlowDirection, dGdMacroCohesion, tmpVec );

        if ( error ){
            errorOut result = new errorNode( "computeFlowDirections",
                                             "Error in the computation of the macro flow direction" );
            result->addNext( error );
            return result;
        }

        parameterType microFrictionAngle = microFlowParameters[ 0 ];
        parameterType microBeta          = microFlowParameters[ 1 ];
        parameterType microFlowPotential;

        error = computeSecondOrderDruckerPragerYieldEquation( referenceMicroStress, microCohesion, elasticRightCauchyGreen,
                                                              microFrictionAngle, microBeta, microFlowPotential, 
                                                              microFlowDirection, dGdMicroCohesion, tmpVec );

        if ( error ){
            errorOut result = new errorNode( "computeFlowDirections",
                                             "Error in the computation of the macro flow direction" );
            result->addNext( error );
            return result;
        }

        parameterType microGradientFrictionAngle = microGradientFlowParameters[ 0 ];
        parameterType microGradientBeta          = microGradientFlowParameters[ 1 ];
        parameterVector microGradientFlowPotential;

        variableMatrix _microGradientFlowDirection;
        error = computeHigherOrderDruckerPragerYieldEquation( referenceHigherOrderStress, microGradientCohesion, elasticRightCauchyGreen,
                                                              microGradientFrictionAngle, microGradientBeta, microGradientFlowPotential,
                                                              _microGradientFlowDirection, dGdMicroGradientCohesion, tmpMat ); 

        microGradientFlowDirection = variableVector( dim * dim * dim * dim, 0 );
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        microGradientFlowDirection[ dim * dim * dim * I + dim * dim * K + dim * L + M ]
                            = _microGradientFlowDirection[ I ][ dim * dim * K + dim * L + M ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeFlowDirections( const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                    const variableVector &referenceHigherOrderStress, const variableType &macroCohesion,
                                    const variableType &microCohesion, const variableVector &microGradientCohesion,
                                    const variableVector &elasticRightCauchyGreen, const parameterVector &macroFlowParameters,
                                    const parameterVector &microFlowParameters, const parameterVector &microGradientFlowParameters,
                                    variableVector &macroFlowDirection, variableVector &microFlowDirection,
                                    variableVector &microGradientFlowDirection, variableType &dGdMacroCohesion,
                                    variableType &dGdMicroCohesion, variableMatrix &dGdMicroGradientCohesion,
                                    variableMatrix &dMacroFlowDirectiondPK2Stress,
                                    variableMatrix &dMacroFlowDirectiondElasticRCG,
                                    variableMatrix &dMicroFlowDirectiondReferenceMicroStress,
                                    variableMatrix &dMicroFlowDirectiondElasticRCG,
                                    variableMatrix &dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                    variableMatrix &dMicroGradientFlowDirectiondElasticRCG ){
        /*!
         * Compute all of the flow directions.
         *
         * :param const variableVector &PK2Stress: The second Piola-Kirchhoff stress.
         * :param const variableVector &referenceMicroStress: The reference symmetric micro stress.
         * :param const variableVector &referenceHigherOrderStress: The reference higher order stress.
         * :param const variableType &macroCohesion: The macro cohesion value.
         * :param const variableType &microCohesion: The micro cohesion value.
         * :param const variableVector &microGradientCohesion: The micro gradient cohesion value.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right cauchy green deformation.
         * :param const parameterVector &macroFlowParameters: The macro plastic parameters.
         *     [ friction angle, beta ]
         * :param const parameterVector &microFlowParameters: The micro plastic parameters.
         *     [ friction angle, beta ]
         * :param const parameterVector &microGradientFlowParameters: The micro gradient plastic parameters.
         *     [ friction angle, beta ]
         * :param variableVector &macroFlowDirection: The flow direction for the macro scale plasticity.
         * :param variableVector &microFlowDirection: The flow direction for the micro scale plasticity.
         * :param variableVector &microGradientFlowDirection: The flow direction for the micro gradient 
         *     plasticity.
         * :param variableType &dGdMacroCohesion: The Jacobian of the macro plastic potential w.r.t. the 
         *     macro cohesion.
         * :param variableType &dGdMicroCohesion: The Jacobian of the micro plastic potential w.r.t. the 
         *     micro cohesion.
         * :param variableVector &dGdMicroGradientCohesion: The Jacobian of the micro gradient plastic 
         *     potential w.r.t. the micro gradient cohesion.
         * :param variableMatrix &dMacroFlowDirectiondPK2Stress: The Jacobian of the macro flow direction w.r.t.
         *     the PK2 stress.
         * :param variableMatrix &dMacroFlowDirectiondElasticRCG: The Jacobian of the macro flow direction w.r.t.
         *     the elastic right Green-Lagrange deformation tensor.
         * :param variableMatrix &dMicroFlowDirectiondReferenceMicroStress: The Jacobian of the micro flow direction w.r.t.
         *     the reference symmetric micro stress.
         * :param variableMatrix &dMicroFlowDirectiondElasticRCG: The Jacobian of the micro flow direction w.r.t.
         *     the elastic right Green-Lagrange deformation tensor.
         * :param variableMatrix &dMicroGradientFlowDirectiondReferenceHigherOrderStress: The Jacobian of the micro gradient flow 
         *     direction w.r.t.the reference higher order stress.
         * :param variableMatrix &dMicroGradientFlowDirectiondElasticRCG: The Jacobian of the micro gradient flow direction w.r.t.
         *     the elastic right Green-Lagrange deformation tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Error handling
        if ( macroFlowParameters.size() != 2 ){
            return new errorNode( "computeFlowDirections",
                                  "The number of macro flow parameters must be 2" );
        }

        if ( microFlowParameters.size() != 2 ){
            return new errorNode( "computeFlowDirections",
                                  "The number of micro flow parameters must be 2" );
        }

        if ( microGradientFlowParameters.size() != 2 ){
            return new errorNode( "computeFlowDirections",
                                  "The number of micro gradient flow parameters must be 2" );
        }

        variableVector tmpVec;
        variableMatrix tmpMat;

        parameterType macroFrictionAngle = macroFlowParameters[ 0 ];
        parameterType macroBeta          = macroFlowParameters[ 1 ];
        parameterType macroFlowPotential;

        errorOut error;

        error = computeSecondOrderDruckerPragerYieldEquation( PK2Stress, macroCohesion, elasticRightCauchyGreen,
                                                              macroFrictionAngle, macroBeta, macroFlowPotential, 
                                                              macroFlowDirection, dGdMacroCohesion, tmpVec,
                                                              dMacroFlowDirectiondPK2Stress, dMacroFlowDirectiondElasticRCG );

        if ( error ){
            errorOut result = new errorNode( "computeFlowDirections",
                                             "Error in the computation of the macro flow direction" );
            result->addNext( error );
            return result;
        }

        parameterType microFrictionAngle = microFlowParameters[ 0 ];
        parameterType microBeta          = microFlowParameters[ 1 ];
        parameterType microFlowPotential;

        error = computeSecondOrderDruckerPragerYieldEquation( referenceMicroStress, microCohesion, elasticRightCauchyGreen,
                                                              microFrictionAngle, microBeta, microFlowPotential, 
                                                              microFlowDirection, dGdMicroCohesion, tmpVec,
                                                              dMicroFlowDirectiondReferenceMicroStress,
                                                              dMicroFlowDirectiondElasticRCG );

        if ( error ){
            errorOut result = new errorNode( "computeFlowDirections",
                                             "Error in the computation of the macro flow direction" );
            result->addNext( error );
            return result;
        }

        parameterType microGradientFrictionAngle = microGradientFlowParameters[ 0 ];
        parameterType microGradientBeta          = microGradientFlowParameters[ 1 ];
        parameterVector microGradientFlowPotential;

        variableMatrix tmp1, tmp2;

        variableMatrix _microGradientFlowDirection;

        error = computeHigherOrderDruckerPragerYieldEquation( referenceHigherOrderStress, microGradientCohesion, elasticRightCauchyGreen,
                                                              microGradientFrictionAngle, microGradientBeta, microGradientFlowPotential,
                                                              _microGradientFlowDirection, dGdMicroGradientCohesion, tmpMat,
                                                              tmp1, tmp2 );

        microGradientFlowDirection = variableVector( dim * dim * dim * dim, 0 );
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        microGradientFlowDirection[ dim * dim * dim * I + dim * dim * K + dim * L + M ]
                            = _microGradientFlowDirection[ I ][ dim * dim * K + dim * L + M ];
                    }
                }
            }
        }

        //Reform the jacobians to a matrix form which enables matrix multiplication of the jacobians
        dMicroGradientFlowDirectiondReferenceHigherOrderStress = variableMatrix( dim * dim * dim * dim,
                                                                 variableVector( dim * dim * dim, 0 ) );
        dMicroGradientFlowDirectiondElasticRCG = variableMatrix( dim * dim * dim * dim,
                                                 variableVector( dim * dim, 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        for ( unsigned int m = 0; m < dim; m++ ){
                            for ( unsigned int n = 0; n < dim; n++ ){
                                dMicroGradientFlowDirectiondElasticRCG[ dim * dim * dim * i + dim * dim * j + dim * k + l ][ dim * m + n ]
                                    = tmp2[ i ][ dim * dim * dim * dim * j + dim * dim * dim * k + dim * dim * l + dim * m + n ];
                                for ( unsigned int o = 0; o < dim; o++ ){
                                    dMicroGradientFlowDirectiondReferenceHigherOrderStress[ dim * dim * dim * i + dim * dim * j + dim * k + l ][ dim * dim * m + dim * n + o ]
                                        = tmp1[ i ][ dim * dim * dim * dim * dim * j + dim * dim * dim * dim * k + dim * dim * dim * l + dim * dim * m + dim * n + o ];
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticDeformationResidual( const tardigradeSolverTools::floatVector &x, const tardigradeSolverTools::floatMatrix &floatArgs,
                                                const tardigradeSolverTools::intMatrix &intArgs, tardigradeSolverTools::floatVector &residual,
                                                tardigradeSolverTools::floatMatrix &jacobian, tardigradeSolverTools::floatMatrix &floatOuts,
                                                tardigradeSolverTools::intMatrix &intOuts
#ifdef DEBUG_MODE
                                                , tardigradeSolverTools::debugMap &DEBUG
#endif
                                              ){
        /*!
         * Compute the residual on the amount of plastic deformation.
         * 
         * :param tardigradeSolverTools::floatVector &x: The unknown vector. Organized as
         *     [ plasticDeformationGradient, plasticMicroDeformation, plasticGradientMicroDeformation,
         *       plasticMacroStrainISV, plasticMicroStrainISV, plasticMicroGradientStrainISV,
         *       currentMacroGamma, currentMicroGamma, currentMicroGradientGamma ]
         * :param const tardigradeSolverTools::floatMatrix &floatArgs: The floating point arguments which do not vary
         *     during the solve.
         * :param const tardigradeSolverTools::intMatrix &intArgs: The integer arguments which do not vary during 
         *     the solve.
         * :param tardigradeSolverTools::floatVector &residual: The value of the residual. This will be the 
         *     the values in the x vector - the estimated amount of plastic deformation
         * :param tardigradeSolverTools::floatMatrix &jacobian: The jacobian matrix
         * :param tardigradeSolverTools::floatMatrix &floatOuts: The floating point values that do change during the solve.
         * :param tardigradeSolverTools::intMatrix &intOuts: The integer values that do change during the solve.
         * :param std::map< std::string, tardigradeSolverTools::floatVector > &DEBUG: The debug map. Only available if
         *     DEBUG_MODE is defined.
         *
         * Ordering of floatArgs
         * floatArgs[  0 ] = Dt, The change in time
         * floatArgs[  1 ] = currentDeformationGradient, The current value of the deformation gradient
         * floatArgs[  2 ] = currentMicroDeformation, The current value of the micro-deformation
         * floatArgs[  3 ] = currentGradientMicroDeformation, The current value of the gradient of the
         *     micro deformation in the reference configuration.
         * floatArgs[  4 ] = currentMacroGamma, The current value of the macro plastic multiplier
         * floatArgs[  5 ] = currentMicroGamma, The current value of the micro plastic multiplier
         * floatArgs[  6 ] = currentMicroGradientGamma, The current values of the micro plastic
         *     multipliers
         * floatArgs[  7 ] = previousMacroGamma, The previous value of the macro plastic multiplier
         * floatArgs[  8 ] = previousMicroGamma, The previous value of the micro plastic multiplier
         * floatArgs[  9 ] = previousMicroGradientGamma, The previous values of the micro plastic
         *     multipliers
         * floatArgs[ 10 ] = previousPlasticDeformationGradient, The previous value of the plastic
         *     deformation gradient.
         * floatArgs[ 11 ] = previousPlasticMicroDeformation, The previous value of the plastic micro
         *     deforamtion.
         * floatArgs[ 12 ] = previousPlasticGradientMicroDeformation, The previous value of the plastic
         *     intermediate configuration gradient of the micro deformation.
         * floatArgs[ 13 ] = previousMacroStrainISV, The previous value of the macro strain-like ISV.
         * floatArgs[ 14 ] = previousMicroStrainISV, The previous value of the micro strain-like ISV.
         * floatArgs[ 15 ] = previousMicroGradientStrainISV, The previous values of the micro gradient
         *     strain-like ISVs.
         * floatArgs[ 16 ] = previousdMacroGdMicroCohesion, The previous value of the Jacobian of the
         *     macro flow direction w.r.t. the macro cohesion value.
         * floatArgs[ 17 ] = previousdMicroGdMicroCohesion, The previous value of the Jacobian of the
         *     micro flow direction w.r.t. the micro cohesion value.
         * floatArgs[ 18 ] = previousdMicroGradientGdMicroGradientCohesion, The previous value of the
         *     Jacobian of the micro gradient flow direction w.r.t. the micro gradient cohesion value.
         * floatArgs[ 19 ] = previousPlasticMacroVelocityGradient, The previous plastic macro velocity
         *     gradient.
         * floatArgs[ 20 ] = previousPlasticMicroVelocityGradient, The previous plastic micro velocity
         *     gradient.
         * floatArgs[ 21 ] = previousPlasticMicroGradientVelocityGradient, The previous plastic micro
         *     gradient velocity gradient.
         * floatArgs[ 22 ] = macroFlowParameters, The macro flow parameters.
         * floatArgs[ 23 ] = microFlowParameters, The micro flow parameters.
         * floatArgs[ 24 ] = microGradientFlowParameters, The micro gradient flow parameters.
         * floatArgs[ 25 ] = macroHardeningParameters, The macro hardening parameters.
         * floatArgs[ 26 ] = microHardeningParameters, The micro hardening parameters.
         * floatArgs[ 27 ] = microGradientHardeningParameters, The micro gradient hardening parameters.
         * floatArgs[ 28 ] = Amatrix, The A stiffness tensor.
         * floatArgs[ 29 ] = Bmatrix, The B stiffness tensor.
         * floatArgs[ 30 ] = Cmatrix, The C stiffness tensor.
         * floatArgs[ 31 ] = Dmatrix, The D stiffness tensor.
         * floatArgs[ 32 ] = alphaMacro, The macro integration parameter.
         * floatArgs[ 33 ] = alphaMicro, The micro integration parameter.
         * floatArgs[ 34 ] = alphaMicroGradient, The micro gradient integration parameter.
         *
         * Ordering of intArgs
         * intArgs[ 0 ] = computeFullDerivatives, Flag which indicates if the full derivatives should be
         *     computed.
         * 
         * Ordering of floatOuts
         * floatOuts[  0 ] = currentPK2Stress, The current value of the second Piola-Kirchoff stress
         * floatOuts[  1 ] = currentReferenceMicroStress, The current value of the reference micro
         *     stress.
         * floatOuts[  2 ] = currentReferenceHigherOrderStress, The current value of the reference
         *     higher order stress.
         * floatOuts[  3 ] = currentMacroStrainISV, The current value of the macro strain-like ISV
         * floatOuts[  4 ] = currentMacroStrainISV, The current value of the micro strain-like ISV
         * floatOuts[  5 ] = currentMacroStrainISV, The current value of the micro gradient strain-like ISV
         * floatOuts[  6 ] = currentMacroCohesion, The current value of the macro cohesion
         * floatOuts[  7 ] = currentMicroCohesion, The current value of the micro cohesion
         * floatOuts[  8 ] = currentMicroGradientCohesion, The current value of the micro gradient cohesion
         * floatOuts[  9 ] = currentElasticRightCauchyGreen, The current value of the right Cauchy green
         * floatOuts[ 10 ] = dElasticRightCauchyGreendPlasticDeformationGradient, The partial derivative
         *     of the stresses w.r.t. the solution vector.
         * floatOuts[ 11 ] = dStressdx, The partial derivative of the stresses w.r.t. the 
         *     solution vector.
         * floatOuts[ 12 ] = dStressdDeformation, The partial derivative of the stresses w.r.t.
         *     the fundamental deformation measures.
         * floatOuts[ 13 ] = dResidualdDeformation, The partial derivative of the residual w.r.t.
         *     the fundamental deformation measures.
         * floatOuts[ 14 ] = dElasticRightCauchyGreendDeformation, The partial derivative of the
         *     elastic right cauchy green w.r.t. the deformation gradient.
         * floatOuts[ 15 ] = dResidualdGammas, The partial derivative of the residual w.r.t. the 
         *     fundamental plastic multipliers.
         * floatOuts[ 16 ] = dCohesionsdGammas, The partial derivative of the residual w.r.t. the
         *     plastic multipliers
         */

        if ( x.size() != 45 ){
            return new errorNode( "computePlasticDeformationResidual",
                                  "The x vector must have a length of 45" );
        }

        if ( floatArgs.size() != 35 ){
            return new errorNode( "computePlasticDeformationResidual",
                                  "The floating point argument matrix floatArgs must have a length of 35" );
        }

        if ( intArgs.size() != 1 ){
            return new errorNode( "computePlasticDeformationResidual",
                                  "The integer argument matrix intArgs must have a length of 1" );
        }

        if ( intArgs[ 0 ].size() != 1 ){
            return new errorNode( "computePlasticDeformationResidual",
                                  "The deformation measure evaluation flag must have a size of 1" );
        }

        bool evaluateFullDerivatives = false;
        if ( intArgs[ 0 ][ 0 ] > 0 ){
            evaluateFullDerivatives = true;
        }

        if ( evaluateFullDerivatives ){
            if ( floatOuts.size() != 17 ){
                return new errorNode( "computePlasticDeformationResidual",
                                      "The floating point output matrix floatOuts must have a length of 17" );
            }
        }
        else{ 
            if ( floatOuts.size() != 10 ){
                return new errorNode( "computePlasticDeformationResidual",
                                      "The floating point output matrix floatOuts must have a length of 10" );
            }
        }

        if ( intOuts.size() != 0 ){
            return new errorNode( "computePlasticDeformationResidual",
                                  "The integer output matrix intOuts must have a length of 0" );
        }

        /*=============================
        | Extract the incoming values |
        =============================*/

        const variableVector currentPlasticDeformationGradient( x.begin(), x.begin() + 9 );
        const variableVector currentPlasticMicroDeformation( x.begin() + 9, x.begin() + 18 );
        const variableVector currentPlasticGradientMicroDeformation( x.begin() + 18, x.begin() + 45 );

        unsigned int ii = 0;
        const constantType    *Dt                                            = &floatArgs[ ii++ ][ 0 ];
        const variableVector  *currentDeformationGradient                    = &floatArgs[ ii++ ];
        const variableVector  *currentMicroDeformation                       = &floatArgs[ ii++ ];
        const variableVector  *currentGradientMicroDeformation               = &floatArgs[ ii++ ];
        const variableType    *currentMacroGamma                             = &floatArgs[ ii++ ][ 0 ];
        const variableType    *currentMicroGamma                             = &floatArgs[ ii++ ][ 0 ];
        const variableVector  *currentMicroGradientGamma                     = &floatArgs[ ii++ ];
        const variableType    *previousMacroGamma                            = &floatArgs[ ii++ ][ 0 ];
        const variableType    *previousMicroGamma                            = &floatArgs[ ii++ ][ 0 ];
        const variableVector  *previousMicroGradientGamma                    = &floatArgs[ ii++ ];
        const variableVector  *previousPlasticDeformationGradient            = &floatArgs[ ii++ ];
        const variableVector  *previousPlasticMicroDeformation               = &floatArgs[ ii++ ];
        const variableVector  *previousPlasticGradientMicroDeformation       = &floatArgs[ ii++ ];
        const variableType    *previousMacroStrainISV                        = &floatArgs[ ii++ ][ 0 ];
        const variableType    *previousMicroStrainISV                        = &floatArgs[ ii++ ][ 0 ];
        const variableVector  *previousMicroGradientStrainISV                = &floatArgs[ ii++ ];
        const variableType    *previousdMacroGdMacroCohesion                 = &floatArgs[ ii++ ][ 0 ];
        const variableType    *previousdMicroGdMicroCohesion                 = &floatArgs[ ii++ ][ 0 ];
        const variableMatrix   previousdMicroGradientGdMicroGradientCohesion = tardigradeVectorTools::inflate( floatArgs[ ii++ ], 3, 3 );
        const variableVector  *previousPlasticMacroVelocityGradient          = &floatArgs[ ii++ ];
        const variableVector  *previousPlasticMicroVelocityGradient          = &floatArgs[ ii++ ];
        const variableVector  *previousPlasticMicroGradientVelocityGradient  = &floatArgs[ ii++ ];
        const parameterVector *macroFlowParameters                           = &floatArgs[ ii++ ];
        const parameterVector *microFlowParameters                           = &floatArgs[ ii++ ];
        const parameterVector *microGradientFlowParameters                   = &floatArgs[ ii++ ];
        const parameterVector *macroHardeningParameters                      = &floatArgs[ ii++ ];
        const parameterVector *microHardeningParameters                      = &floatArgs[ ii++ ];
        const parameterVector *microGradientHardeningParameters              = &floatArgs[ ii++ ];
//        const parameterVector *macroYieldParameters                          = &floatArgs[ ii++ ];
//        const parameterVector *microYieldParameters                          = &floatArgs[ ii++ ];
//        const parameterVector *microGradientYieldParameters                  = &floatArgs[ ii++ ];
        const parameterVector *Amatrix                                       = &floatArgs[ ii++ ];
        const parameterVector *Bmatrix                                       = &floatArgs[ ii++ ];
        const parameterVector *Cmatrix                                       = &floatArgs[ ii++ ];
        const parameterVector *Dmatrix                                       = &floatArgs[ ii++ ];
        const parameterType   *alphaMacro                                    = &floatArgs[ ii++ ][ 0 ];
        const parameterType   *alphaMicro                                    = &floatArgs[ ii++ ][ 0 ];
        const parameterType   *alphaMicroGradient                            = &floatArgs[ ii++ ][ 0 ];

#ifdef DEBUG_MODE
        variableVector temp = { *currentMacroGamma };
        DEBUG.emplace( "currentMacroGamma", temp );

        temp = { *currentMicroGamma };
        DEBUG.emplace( "currentMicroGamma", temp );

        DEBUG.emplace( "currentMicroGradientGamma", *currentMicroGradientGamma );
#endif

        /*=================================
        | Compute the Elastic Deformation |
        =================================*/

        variableVector currentElasticDeformationGradient, currentElasticMicroDeformation, currentElasticGradientMicroDeformation;

        variableMatrix dElasticDeformationGradientdDeformationGradient, dElasticDeformationGradientdPlasticDeformationGradient,
                       dElasticMicroDeformationdMicroDeformation, dElasticMicroDeformationdPlasticMicroDeformation,
                       dElasticGradientMicroDeformationdGradientMicroDeformation,
                       dElasticGradientMicroDeformationdPlasticGradientMicroDeformation,
                       dElasticGradientMicroDeformationdPlasticDeformationGradient,
                       dElasticGradientMicroDeformationdMicroDeformation,
                       dElasticGradientMicroDeformationdPlasticMicroDeformation;

        errorOut error = computeElasticPartOfDeformation( *currentDeformationGradient, *currentMicroDeformation,
                                                          *currentGradientMicroDeformation,
                                                           currentPlasticDeformationGradient, currentPlasticMicroDeformation,
                                                           currentPlasticGradientMicroDeformation,
                                                           currentElasticDeformationGradient, currentElasticMicroDeformation,
                                                           currentElasticGradientMicroDeformation,
                                                           dElasticDeformationGradientdDeformationGradient,
                                                           dElasticDeformationGradientdPlasticDeformationGradient,
                                                           dElasticMicroDeformationdMicroDeformation,
                                                           dElasticMicroDeformationdPlasticMicroDeformation,
                                                           dElasticGradientMicroDeformationdGradientMicroDeformation,
                                                           dElasticGradientMicroDeformationdPlasticGradientMicroDeformation,
                                                           dElasticGradientMicroDeformationdPlasticDeformationGradient,
                                                           dElasticGradientMicroDeformationdMicroDeformation,
                                                           dElasticGradientMicroDeformationdPlasticMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the elastic part of the deformation" );
            result->addNext( error );
            return result;
        }

//        std::cout << "IN PLASTIC DEFORMATION RESIDUAL\n";
//        std::cout << "\n  gammas:\n";
//        std::cout << "    " << currentMacroGamma << "\n";
//        std::cout << "    " << currentMicroGamma << "\n";
//        std::cout << "    "; tardigradeVectorTools::print( currentMicroGradientGamma );
//        std::cout << "\n  ISVS:\n";
//        std::cout << "    " << currentMacroStrainISV << "\n";
//        std::cout << "    " << currentMicroStrainISV << "\n";
//        std::cout << "    "; tardigradeVectorTools::print( currentMicroGradientStrainISV ); 
//        std::cout << "  elastic fundamental deformation measures\n";
//        std::cout << "    "; tardigradeVectorTools::print( currentElasticDeformationGradient );
//        std::cout << "    "; tardigradeVectorTools::print( currentElasticMicroDeformation );
//        std::cout << "    "; tardigradeVectorTools::print( currentElasticGradientMicroDeformation );

#ifdef DEBUG_MODE

        //Save the elastic fundamental deformation measures
        DEBUG.emplace( "currentElasticDeformationGradient", currentElasticDeformationGradient );
        DEBUG.emplace( "currentElasticMicroDeformation", currentElasticMicroDeformation );
        DEBUG.emplace( "currentElasticGradientMicroDeformation", currentElasticGradientMicroDeformation );

        //Save the flattened jacobians
        DEBUG.emplace( "dElasticDeformationGradientdDeformationGradient",
                       tardigradeVectorTools::appendVectors( dElasticDeformationGradientdDeformationGradient )  );
        DEBUG.emplace( "dElasticDeformationGradientdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dElasticDeformationGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticMicroDeformationdMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticMicroDeformationdMicroDeformation ) );
        DEBUG.emplace( "dElasticMicroDeformationdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticMicroDeformationdPlasticMicroDeformation ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticGradientMicroDeformationdGradientMicroDeformation ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticGradientMicroDeformationdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dElasticGradientMicroDeformationdPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticGradientMicroDeformationdMicroDeformation ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticGradientMicroDeformationdPlasticMicroDeformation ) );

#endif

        /*!=================================================
        | Compute the elastic derived deformation measures |
        ==================================================*/

        variableVector currentElasticRightCauchyGreen, currentElasticMicroRightCauchyGreen, currentElasticPsi, currentElasticGamma;
        variableMatrix dElasticRightCauchyGreendElasticDeformationGradient,
                       dElasticMicroRightCauchyGreendElasticMicroDeformation,
                       dElasticPsidElasticDeformationGradient,
                       dElasticPsidElasticMicroDeformation,
                       dElasticGammadElasticDeformationGradient,
                       dElasticGammadElasticGradientMicroDeformation;

        error = computeElasticDeformationMeasures( currentElasticDeformationGradient, currentElasticMicroDeformation,
                                                   currentElasticGradientMicroDeformation,
                                                   currentElasticRightCauchyGreen, currentElasticMicroRightCauchyGreen,
                                                   currentElasticPsi, currentElasticGamma,
                                                   dElasticRightCauchyGreendElasticDeformationGradient,
                                                   dElasticMicroRightCauchyGreendElasticMicroDeformation,
                                                   dElasticPsidElasticDeformationGradient,
                                                   dElasticPsidElasticMicroDeformation,
                                                   dElasticGammadElasticDeformationGradient,
                                                   dElasticGammadElasticGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the derived elastic deformation measures" );
            result->addNext( error );
            return result;
        }

        //Save the current Elastic Right Cauchy Green Deformation Tensor
        floatOuts[ 9 ] = currentElasticRightCauchyGreen;

//        std::cout << "\n  elastic derived deformation measures\n";
//        std::cout << "    "; tardigradeVectorTools::print( currentElasticRightCauchyGreen );
//        std::cout << "    "; tardigradeVectorTools::print( currentElasticMicroRightCauchyGreen );
//        std::cout << "    "; tardigradeVectorTools::print( currentElasticPsi );
//        std::cout << "    "; tardigradeVectorTools::print( currentElasticGamma );

        /*!==================================================================
        | Assemble the Jacobian of the elastic derived deformation measures |
        ===================================================================*/

        //Compute the Jacobians w.r.t. the plastic deformation
        variableMatrix dElasticRightCauchyGreendPlasticDeformationGradient
            = tardigradeVectorTools::dot( dElasticRightCauchyGreendElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient );

        variableMatrix dElasticMicroRightCauchyGreendPlasticMicroDeformation
            = tardigradeVectorTools::dot( dElasticMicroRightCauchyGreendElasticMicroDeformation,
                                dElasticMicroDeformationdPlasticMicroDeformation );

        variableMatrix dElasticPsidPlasticDeformationGradient
            = tardigradeVectorTools::dot( dElasticPsidElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient );

        variableMatrix dElasticPsidPlasticMicroDeformation
            = tardigradeVectorTools::dot( dElasticPsidElasticMicroDeformation,
                                dElasticMicroDeformationdPlasticMicroDeformation );

        variableMatrix dElasticGammadPlasticDeformationGradient
            = tardigradeVectorTools::dot( dElasticGammadElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticDeformationGradient );

        variableMatrix dElasticGammadPlasticMicroDeformation
            = tardigradeVectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticMicroDeformation );

        variableMatrix dElasticGammadPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticGradientMicroDeformation );

        //Construct the full derivatives
        variableMatrix dElasticRightCauchyGreendDeformationGradient, dElasticMicroRightCauchyGreendMicroDeformation,
                       dElasticPsidDeformationGradient, dElasticPsidMicroDeformation,
                       dElasticGammadDeformationGradient, dElasticGammadMicroDeformation, dElasticGammadGradientMicroDeformation;
        if ( evaluateFullDerivatives ){
            dElasticRightCauchyGreendDeformationGradient = tardigradeVectorTools::dot( dElasticRightCauchyGreendElasticDeformationGradient,
                                                                             dElasticDeformationGradientdDeformationGradient );
            dElasticMicroRightCauchyGreendMicroDeformation = tardigradeVectorTools::dot( dElasticMicroRightCauchyGreendElasticMicroDeformation,
                                                                               dElasticMicroDeformationdMicroDeformation );
            dElasticPsidDeformationGradient = tardigradeVectorTools::dot( dElasticPsidElasticDeformationGradient,
                                                                dElasticDeformationGradientdDeformationGradient );
            dElasticPsidMicroDeformation = tardigradeVectorTools::dot( dElasticPsidElasticMicroDeformation,
                                                             dElasticMicroDeformationdMicroDeformation );
            dElasticGammadDeformationGradient = tardigradeVectorTools::dot( dElasticGammadElasticDeformationGradient,
                                                                  dElasticDeformationGradientdDeformationGradient );
            dElasticGammadMicroDeformation = tardigradeVectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                                               dElasticGradientMicroDeformationdMicroDeformation );
            dElasticGammadGradientMicroDeformation = tardigradeVectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                                                       dElasticGradientMicroDeformationdGradientMicroDeformation );

            floatOuts[ 10 ] = tardigradeVectorTools::appendVectors( dElasticRightCauchyGreendPlasticDeformationGradient );
            floatOuts[ 14 ] = tardigradeVectorTools::appendVectors( dElasticRightCauchyGreendDeformationGradient );
        }

#ifdef DEBUG_MODE

        //Save the elastic derived deformation measures
        DEBUG.emplace( "currentElasticRightCauchyGreen", currentElasticRightCauchyGreen );
        DEBUG.emplace( "currentElasticMicroRightCauchyGreen", currentElasticMicroRightCauchyGreen );
        DEBUG.emplace( "currentElasticPsi", currentElasticPsi );
        DEBUG.emplace( "currentElasticGamma", currentElasticGamma );

        //Save the Jacobians
        DEBUG.emplace( "dElasticRightCauchyGreendPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dElasticRightCauchyGreendPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticMicroRightCauchyGreendPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticMicroRightCauchyGreendPlasticMicroDeformation ) );
        DEBUG.emplace( "dElasticPsidPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dElasticPsidPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticPsidPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticPsidPlasticMicroDeformation ) );
        DEBUG.emplace( "dElasticGammadPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dElasticGammadPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticGammadPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticGammadPlasticMicroDeformation ) );
        DEBUG.emplace( "dElasticGammadPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dElasticGammadPlasticGradientMicroDeformation ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dElasticRightCauchyGreendDeformationGradient",
                            tardigradeVectorTools::appendVectors( dElasticRightCauchyGreendDeformationGradient ) );
            DEBUG.emplace( "dElasticMicroRightCauchyGreendMicroDeformation",
                            tardigradeVectorTools::appendVectors( dElasticMicroRightCauchyGreendMicroDeformation ) );
            DEBUG.emplace( "dElasticPsidDeformationGradient",
                            tardigradeVectorTools::appendVectors( dElasticPsidDeformationGradient ) );
            DEBUG.emplace( "dElasticPsidMicroDeformation",
                            tardigradeVectorTools::appendVectors( dElasticPsidMicroDeformation ) );
            DEBUG.emplace( "dElasticGammadDeformationGradient",
                            tardigradeVectorTools::appendVectors( dElasticGammadDeformationGradient ) );
            DEBUG.emplace( "dElasticGammadMicroDeformation",
                            tardigradeVectorTools::appendVectors( dElasticGammadMicroDeformation ) );
            DEBUG.emplace( "dElasticGammadGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dElasticGammadGradientMicroDeformation ) );
        }

#endif

        /*===========================
        | Compute the Stress values |
        ===========================*/

        variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

        variableMatrix dPK2StressdElasticDeformationGradient, dPK2StressdElasticMicroDeformation,
                       dPK2StressdElasticGradientMicroDeformation,
                       dReferenceMicroStressdElasticDeformationGradient, dReferenceMicroStressdElasticMicroDeformation,
                       dReferenceMicroStressdElasticGradientMicroDeformation,
                       dReferenceHigherOrderStressdElasticDeformationGradient,
                       dReferenceHigherOrderStressdElasticGradientMicroDeformation;

        error = tardigradeMicromorphicLinearElasticity::linearElasticityReference(  currentElasticDeformationGradient,
                                                                          currentElasticMicroDeformation,
                                                                          currentElasticGradientMicroDeformation,
                                                                         *Amatrix, *Bmatrix, *Cmatrix, *Dmatrix,
                                                                          currentPK2Stress, currentReferenceMicroStress,
                                                                          currentReferenceHigherOrderStress,
                                                                          dPK2StressdElasticDeformationGradient,
                                                                          dPK2StressdElasticMicroDeformation,
                                                                          dPK2StressdElasticGradientMicroDeformation,
                                                                          dReferenceMicroStressdElasticDeformationGradient,
                                                                          dReferenceMicroStressdElasticMicroDeformation,
                                                                          dReferenceMicroStressdElasticGradientMicroDeformation,
                                                                          dReferenceHigherOrderStressdElasticDeformationGradient,
                                                                          dReferenceHigherOrderStressdElasticGradientMicroDeformation );
        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the stresses" );
            result->addNext( error );
            return result;
        }

//        std::cout << "\n  current stress measures\n";
//        std::cout << "    "; tardigradeVectorTools::print( currentPK2Stress );
//        std::cout << "    "; tardigradeVectorTools::print( currentReferenceMicroStress );
//        std::cout << "    "; tardigradeVectorTools::print( currentReferenceHigherOrderStress );

        /*===============================
        | Assemble the stress Jacobians |
        ===============================*/

        //Jacobians w.r.t. the plastic deformation
        variableMatrix dPK2StressdPlasticDeformationGradient
            = tardigradeVectorTools::dot( dPK2StressdElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticDeformationGradient );

        variableMatrix dPK2StressdPlasticMicroDeformation
            = tardigradeVectorTools::dot( dPK2StressdElasticMicroDeformation,
                                dElasticMicroDeformationdPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticMicroDeformation );

        variableMatrix dPK2StressdPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticGradientMicroDeformation );

        variableMatrix dReferenceMicroStressdPlasticDeformationGradient
            = tardigradeVectorTools::dot( dReferenceMicroStressdElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticDeformationGradient );

        variableMatrix dReferenceMicroStressdPlasticMicroDeformation
            = tardigradeVectorTools::dot( dReferenceMicroStressdElasticMicroDeformation,
                                dElasticMicroDeformationdPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticMicroDeformation );

        variableMatrix dReferenceMicroStressdPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticGradientMicroDeformation );

        variableMatrix dReferenceHigherOrderStressdPlasticDeformationGradient
            = tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticDeformationGradient );

        variableMatrix dReferenceHigherOrderStressdPlasticMicroDeformation
            = tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticMicroDeformation );

        variableMatrix dReferenceHigherOrderStressdPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticGradientMicroDeformation );

        //Jacobians w.r.t. the fundamental deformation measures
        variableMatrix dPK2StressdDeformationGradient, dPK2StressdMicroDeformation, dPK2StressdGradientMicroDeformation,
                       dReferenceMicroStressdDeformationGradient, dReferenceMicroStressdMicroDeformation,
                       dReferenceMicroStressdGradientMicroDeformation,
                       dReferenceHigherOrderStressdDeformationGradient, dReferenceHigherOrderStressdMicroDeformation,
                       dReferenceHigherOrderStressdGradientMicroDeformation;

        if ( evaluateFullDerivatives ){
            dPK2StressdDeformationGradient = tardigradeVectorTools::dot( dPK2StressdElasticDeformationGradient,
                                                               dElasticDeformationGradientdDeformationGradient );
            dPK2StressdMicroDeformation = tardigradeVectorTools::dot( dPK2StressdElasticMicroDeformation,
                                                            dElasticMicroDeformationdMicroDeformation )
                                        + tardigradeVectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                                            dElasticGradientMicroDeformationdMicroDeformation );
            dPK2StressdGradientMicroDeformation = tardigradeVectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            dReferenceMicroStressdDeformationGradient = tardigradeVectorTools::dot( dReferenceMicroStressdElasticDeformationGradient,
                                                                          dElasticDeformationGradientdDeformationGradient );
            dReferenceMicroStressdMicroDeformation = tardigradeVectorTools::dot( dReferenceMicroStressdElasticMicroDeformation,
                                                                       dElasticMicroDeformationdMicroDeformation )
                                                   + tardigradeVectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                                                       dElasticGradientMicroDeformationdMicroDeformation );
            dReferenceMicroStressdGradientMicroDeformation
                = tardigradeVectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            dReferenceHigherOrderStressdDeformationGradient
                = tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticDeformationGradient,
                                    dElasticDeformationGradientdDeformationGradient );

            dReferenceHigherOrderStressdMicroDeformation
                = tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdMicroDeformation );

            dReferenceHigherOrderStressdGradientMicroDeformation
                = tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            //Assemble the jacobians into the output vector
            
            floatOuts[ 11 ] = tardigradeSolverTools::floatVector( 45 * 45, 0 ); //Jacobians w.r.t. the solution vector x
            floatOuts[ 12 ] = tardigradeSolverTools::floatVector( 45 * 45, 0 ); //Jacobians w.r.t. the fundamental deformation measures

            //Save the Jacobians of the PK2 and reference symmetric micro stresses
            for ( unsigned int i = 0; i < 9; i++ ){
                for ( unsigned int j = 0; j < 9; j++ ){
                    floatOuts[ 11 ][ 45 * i + j ]     = dPK2StressdPlasticDeformationGradient[ i ][ j ];
                    floatOuts[ 12 ][ 45 * i + j ]     = dPK2StressdDeformationGradient[ i ][ j ];
                    floatOuts[ 11 ][ 45 * i + j + 9 ] = dPK2StressdPlasticMicroDeformation[ i ][ j ];
                    floatOuts[ 12 ][ 45 * i + j + 9 ] = dPK2StressdMicroDeformation[ i ][ j ];

                    floatOuts[ 11 ][ 45 * ( i + 9 ) + j ] = dReferenceMicroStressdPlasticDeformationGradient[ i ][ j ];
                    floatOuts[ 12 ][ 45 * ( i + 9 ) + j ] = dReferenceMicroStressdDeformationGradient[ i ][ j ];
                    floatOuts[ 11 ][ 45 * ( i + 9 ) + j + 9 ] = dReferenceMicroStressdPlasticMicroDeformation[ i ][ j ];
                    floatOuts[ 12 ][ 45 * ( i + 9 ) + j + 9 ] = dReferenceMicroStressdMicroDeformation[ i ][ j ];
                }

                for ( unsigned int j = 0; j < 27; j++ ){
                    floatOuts[ 11 ][ 45 * i + j + 18 ] = dPK2StressdPlasticGradientMicroDeformation[ i ][ j ];
                    floatOuts[ 12 ][ 45 * i + j + 18 ] = dPK2StressdGradientMicroDeformation[ i ][ j ];

                    floatOuts[ 11 ][ 45 * ( i + 9 ) + j + 18 ] = dReferenceMicroStressdPlasticGradientMicroDeformation[ i ][ j ];
                    floatOuts[ 12 ][ 45 * ( i + 9 ) + j + 18 ] = dReferenceMicroStressdGradientMicroDeformation[ i ][ j ];
                }
            }

            //Save the Jacobians of the reference higher order stress
            for ( unsigned int i = 0; i < 27; i++ ){
                for ( unsigned int j = 0; j < 9; j++ ){
                    floatOuts[ 11 ][ 45 * ( i + 18 ) + j ] = dReferenceHigherOrderStressdPlasticDeformationGradient[ i ][ j ];
                    floatOuts[ 12 ][ 45 * ( i + 18 ) + j ] = dReferenceHigherOrderStressdDeformationGradient[ i ][ j ];
                    floatOuts[ 11 ][ 45 * ( i + 18 ) + j + 9 ] = dReferenceHigherOrderStressdPlasticMicroDeformation[ i ][ j ];
                    floatOuts[ 12 ][ 45 * ( i + 18 ) + j + 9 ] = dReferenceHigherOrderStressdMicroDeformation[ i ][ j ];
                }

                for ( unsigned int j = 0; j < 27; j++ ){
                    floatOuts[ 11 ][ 45 * ( i + 18 ) + j + 18 ] = dReferenceHigherOrderStressdPlasticGradientMicroDeformation[ i ][ j ];
                    floatOuts[ 12 ][ 45 * ( i + 18 ) + j + 18 ] = dReferenceHigherOrderStressdGradientMicroDeformation[ i ][ j ];
                }
            }
        }

#ifdef DEBUG_MODE

        //Save the stress values
        DEBUG.emplace( "currentPK2Stress", currentPK2Stress );
        DEBUG.emplace( "currentReferenceMicroStress", currentReferenceMicroStress );
        DEBUG.emplace( "currentReferenceHigherOrderStress", currentReferenceHigherOrderStress );

        //Save the Jacobians

        //Save the Jacobians w.r.t. the plastic deformations
        DEBUG.emplace( "dPK2StressdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dPK2StressdPlasticDeformationGradient ) );
        DEBUG.emplace( "dPK2StressdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dPK2StressdPlasticMicroDeformation ) );
        DEBUG.emplace( "dPK2StressdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dPK2StressdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dReferenceMicroStressdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dReferenceMicroStressdPlasticDeformationGradient ) );
        DEBUG.emplace( "dReferenceMicroStressdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dReferenceMicroStressdPlasticMicroDeformation ) );
        DEBUG.emplace( "dReferenceMicroStressdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dReferenceMicroStressdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dReferenceHigherOrderStressdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dReferenceHigherOrderStressdPlasticDeformationGradient ) );
        DEBUG.emplace( "dReferenceHigherOrderStressdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dReferenceHigherOrderStressdPlasticMicroDeformation ) );
        DEBUG.emplace( "dReferenceHigherOrderStressdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dReferenceHigherOrderStressdPlasticGradientMicroDeformation ) );

        //Save the Jacobians w.r.t. the fundamental deformation measures
        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dPK2StressdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dPK2StressdDeformationGradient ) );
            DEBUG.emplace( "dPK2StressdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPK2StressdMicroDeformation ) );
            DEBUG.emplace( "dPK2StressdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPK2StressdGradientMicroDeformation ) );
    
            DEBUG.emplace( "dReferenceMicroStressdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dReferenceMicroStressdDeformationGradient ) );
            DEBUG.emplace( "dReferenceMicroStressdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dReferenceMicroStressdMicroDeformation ) );
            DEBUG.emplace( "dReferenceMicroStressdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dReferenceMicroStressdGradientMicroDeformation ) );
    
            DEBUG.emplace( "dReferenceHigherOrderStressdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dReferenceHigherOrderStressdDeformationGradient ) );
            DEBUG.emplace( "dReferenceHigherOrderStressdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dReferenceHigherOrderStressdMicroDeformation ) );
            DEBUG.emplace( "dReferenceHigherOrderStressdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dReferenceHigherOrderStressdGradientMicroDeformation ) );
        }

#endif

        /*!=================================================================
        | Compute the derivative of the flow potential w.r.t. the cohesion |
        ==================================================================*/

        //TODO: This is currently assumed to be constant and is only valid for the Drucker-Prager model
        //      This could be expanded to account for non-linear or non-Drucker-Prager behavior in the
        //      future. This would require a non-linear solve for the cohesion.
        //
        //TODO: The computation of the cohesion values should be removed from this function as it has
        //      no interaction with the rest of the function and only with the evolution of Gamma.

        parameterType currentdMacroGdMacroCohesion, currentdMicroGdMicroCohesion,
                      magnitudeCurrentdMicroGradientGdMicroGradientCohesion, tempFloat;
        parameterMatrix currentdMicroGradientGdMicroGradientCohesion = tardigradeVectorTools::eye< parameterType >( 3 );

        error = computeDruckerPragerInternalParameters( (*macroFlowParameters)[ 0 ], (*macroFlowParameters)[ 1 ],
                                                        currentdMacroGdMacroCohesion, tempFloat );

        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the derivative of the macro flow potential w.r.t. the cohesion" );
            result->addNext( error );
            return result;
        }
        currentdMacroGdMacroCohesion *= -1;

        error = computeDruckerPragerInternalParameters( (*microFlowParameters)[ 0 ], (*microFlowParameters)[ 1 ],
                                                        currentdMicroGdMicroCohesion, tempFloat );

        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the derivative of the micro flow potential w.r.t. the cohesion" );
            result->addNext( error );
            return result;
        }
        currentdMicroGdMicroCohesion *= -1;

        error = computeDruckerPragerInternalParameters( (*microGradientFlowParameters)[ 0 ], (*microGradientFlowParameters)[ 1 ],
                                                        magnitudeCurrentdMicroGradientGdMicroGradientCohesion, tempFloat );

        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the derivative of the micro flow potential w.r.t. the cohesion" );
            result->addNext( error );
            return result;
        }

        currentdMicroGradientGdMicroGradientCohesion *= -magnitudeCurrentdMicroGradientGdMicroGradientCohesion;

        /*!=====================================
        | Compute the current strain-like ISVs |
        ======================================*/

        variableType currentMacroStrainISV, currentMicroStrainISV;
        variableVector currentMicroGradientStrainISV;

        variableType dCurrentMacroStrainISVdMacroGamma, dCurrentMacroStrainISVddMacroGdMacroCohesion;
        variableType dCurrentMicroStrainISVdMicroGamma, dCurrentMicroStrainISVddMicroGdMicroCohesion;

        variableMatrix dCurrentMicroGradientStrainISVdMicroGradientGamma,
                       dCurrentMicroGradientStrainISVddMicroGradientGdMicroGradientCohesion;

        error = evolveStrainStateVariables( *Dt,
                                            *currentMacroGamma, *currentMicroGamma, *currentMicroGradientGamma,
                                             currentdMacroGdMacroCohesion, currentdMicroGdMicroCohesion,
                                             currentdMicroGradientGdMicroGradientCohesion,
                                            *previousMacroStrainISV, *previousMicroStrainISV, *previousMicroGradientStrainISV,
                                            *previousMacroGamma, *previousMicroGamma, *previousMicroGradientGamma,
                                            *previousdMacroGdMacroCohesion, *previousdMicroGdMicroCohesion,
                                             previousdMicroGradientGdMicroGradientCohesion,
                                             currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                             dCurrentMacroStrainISVdMacroGamma, dCurrentMacroStrainISVddMacroGdMacroCohesion,
                                             dCurrentMicroStrainISVdMicroGamma, dCurrentMicroStrainISVddMicroGdMicroCohesion,
                                             dCurrentMicroGradientStrainISVdMicroGradientGamma,
                                             dCurrentMicroGradientStrainISVddMicroGradientGdMicroGradientCohesion,
                                            *alphaMacro, *alphaMicro, *alphaMicroGradient );
        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the expected strain-like ISVs" );
            result->addNext( error );
            return result;
        }

#ifdef DEBUG_MODE

        //Save the expected macro-strain values
        temp = { currentMacroStrainISV };
        DEBUG.emplace( "currentMacroStrainISV", temp );

        temp = { currentMicroStrainISV };
        DEBUG.emplace( "currentMicroStrainISV", temp );

        DEBUG.emplace( "currentMicroGradientStrainISV", currentMicroGradientStrainISV );

        //Save the jacobians
        temp = { dCurrentMacroStrainISVdMacroGamma };
        DEBUG.emplace( "dCurrentMacroStrainISVdMacroGamma", temp );

        temp = { dCurrentMacroStrainISVddMacroGdMacroCohesion };
        DEBUG.emplace( "dCurrentMacroStrainISVddMacroGdMacroCohesion", temp );

        temp = { dCurrentMicroStrainISVdMicroGamma };
        DEBUG.emplace( "dCurrentMicroStrainISVdMicroGamma", temp );

        temp = { dCurrentMicroStrainISVddMicroGdMicroCohesion };
        DEBUG.emplace( "dCurrentMicroStrainISVddMicroGdMicroCohesion", temp );

        DEBUG.emplace( "dCurrentMicroGradientStrainISVdMicroGradientGamma",
                        tardigradeVectorTools::appendVectors( dCurrentMicroGradientStrainISVdMicroGradientGamma ) );
        DEBUG.emplace( "dCurrentMicroGradientStrainISVddMicroGradientGdMicroGradientCohesion",
                        tardigradeVectorTools::appendVectors( dCurrentMicroGradientStrainISVddMicroGradientGdMicroGradientCohesion ) );

#endif

        /*!============================
        | Compute the cohesion values |
        =============================*/

        variableType   currentMacroCohesion, currentMicroCohesion;
        variableVector currentMicroGradientCohesion;

        variableType   dCurrentMacroCohesiondMacroStrainISV, dCurrentMicroCohesiondMicroStrainISV;
        variableMatrix dCurrentMicroGradientCohesiondMicroGradientStrainISV;

        error = computeCohesion(  currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                 *macroHardeningParameters, *microHardeningParameters, *microGradientHardeningParameters,
                                  currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion,
                                  dCurrentMacroCohesiondMacroStrainISV, dCurrentMicroCohesiondMicroStrainISV,
                                  dCurrentMicroGradientCohesiondMicroGradientStrainISV );
        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the cohesion values" );
            result->addNext( error );
            return result;
        }

        variableType dCurrentMacroCohesiondMacroGamma = dCurrentMacroCohesiondMacroStrainISV * dCurrentMacroStrainISVdMacroGamma;
        variableType dCurrentMicroCohesiondMicroGamma = dCurrentMicroCohesiondMicroStrainISV * dCurrentMicroStrainISVdMicroGamma;
        variableMatrix dCurrentMicroGradientCohesiondMicroGradientGamma
            = tardigradeVectorTools::dot( dCurrentMicroGradientCohesiondMicroGradientStrainISV,
                                dCurrentMicroGradientStrainISVdMicroGradientGamma );

//        std::cout << "\n  current cohesion values\n";
//        std::cout << "    " << currentMacroCohesion << "\n";
//        std::cout << "    " << currentMicroCohesion << "\n";
//        std::cout << "    "; tardigradeVectorTools::print( currentMicroGradientCohesion );

#ifdef DEBUG_MODE

        //Save the cohesion values
        temp = { currentMacroCohesion };

        DEBUG.emplace( "currentMacroCohesion", temp );

        temp = { currentMicroCohesion };

        DEBUG.emplace( "currentMicroCohesion", temp );

        DEBUG.emplace( "currentMicroGradientCohesion", currentMicroGradientCohesion );

        //Save the cohesion Jacobians
        temp = { dCurrentMacroCohesiondMacroGamma };

        DEBUG.emplace( "dCurrentMacroCohesiondMacroGamma", temp );

        temp = { dCurrentMicroCohesiondMicroGamma };

        DEBUG.emplace( "dCurrentMicroCohesiondMicroGamma", temp );

        DEBUG.emplace( "dCurrentMicroGradientCohesiondMicroGradientGamma",
                       tardigradeVectorTools::appendVectors( dCurrentMicroGradientCohesiondMicroGradientGamma ) );

#endif

        if ( evaluateFullDerivatives ){
            floatOuts[ 16 ].resize( 25, 0 );
            floatOuts[ 16 ][ 0 ] = dCurrentMacroCohesiondMacroGamma;
            floatOuts[ 16 ][ 6 ] = dCurrentMicroCohesiondMicroGamma;
            for ( unsigned int i = 0; i < currentMicroGradientCohesion.size(); i++ ){
                for ( unsigned int j = 0; j < currentMicroGradientGamma->size(); j++ ){
                    floatOuts[ 16 ][ 5 * ( i + 2 ) + j + 2 ] = dCurrentMicroGradientCohesiondMicroGradientGamma[ i ][ j ];
                }
            }
        }

        /*!============================
        | Compute the Flow Directions |
        =============================*/

        variableVector currentMacroFlowDirection, currentMicroFlowDirection, currentMicroGradientFlowDirection;
        variableType _currentdMacroGdMacroCohesion, _currentdMicroGdMicroCohesion; //TODO: This is here for development purposes.
        variableMatrix _currentdMicroGradientGdMicroGradientCohesion;              //TODO: This is here for development purposes.

        variableMatrix dMacroFlowDirectiondPK2Stress, dMacroFlowDirectiondElasticRightCauchyGreen,
                       dMicroFlowDirectiondReferenceMicroStress, dMicroFlowDirectiondElasticRightCauchyGreen,
                       dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                       dMicroGradientFlowDirectiondElasticRightCauchyGreen;

        error = computeFlowDirections(  currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress,
                                        currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion,
                                        currentElasticRightCauchyGreen,
                                       *macroFlowParameters, *microFlowParameters, *microGradientFlowParameters,
                                        currentMacroFlowDirection, currentMicroFlowDirection, currentMicroGradientFlowDirection,
                                        _currentdMacroGdMacroCohesion, _currentdMicroGdMicroCohesion,
                                        _currentdMicroGradientGdMicroGradientCohesion,
                                        dMacroFlowDirectiondPK2Stress, dMacroFlowDirectiondElasticRightCauchyGreen,
                                        dMicroFlowDirectiondReferenceMicroStress, dMicroFlowDirectiondElasticRightCauchyGreen,
                                        dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                        dMicroGradientFlowDirectiondElasticRightCauchyGreen );
        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the flow directions" );
            result->addNext( error );
            return result;
        }

        //TODO: This error section exists for debugging. May not be required once the code is in production.
        if ( !tardigradeVectorTools::fuzzyEquals( _currentdMacroGdMacroCohesion, currentdMacroGdMacroCohesion ) ){
            std::string message     = "The derivative of the macro plastic flow potential w.r.t. the cohesion are not consistent. Report this to Nathan Miller";
            std::string errorValues = "_currentdMacroGdMacroCohesion: " + std::to_string( _currentdMacroGdMacroCohesion ) + "\n"
                                    + " currentdMacroGdMacroCohesion: " + std::to_string( currentdMacroGdMacroCohesion );
            return new errorNode( __func__, message + errorValues );
        }

        if ( !tardigradeVectorTools::fuzzyEquals( _currentdMicroGdMicroCohesion, currentdMicroGdMicroCohesion ) ){
            std::string message     = "The derivative of the micro plastic flow potential w.r.t. the cohesion are not consistent. Report this to Nathan Miller";
            std::string errorValues = "_currentdMicroGdMicroCohesion: " + std::to_string( _currentdMicroGdMicroCohesion ) + "\n"
                                    + " currentdMicroGdMicroCohesion: " + std::to_string( currentdMicroGdMicroCohesion );
            return new errorNode( __func__, message + errorValues );
        }

        if ( !tardigradeVectorTools::fuzzyEquals( _currentdMicroGradientGdMicroGradientCohesion, currentdMicroGradientGdMicroGradientCohesion ) ){
            std::string message     = "The derivative of the micro gradient plastic flow potential w.r.t. the cohesion are not consistent. Report this to Nathan Miller";
            std::string errorValues = "    _currentdMicroGradientGdMicroGradientCohesion:\n";
            for ( auto iter = _currentdMicroGradientGdMicroGradientCohesion.begin( ); iter != _currentdMicroGradientGdMicroGradientCohesion.end( ); iter++ ){
                errorValues += "        ";
                for ( auto iter2 = iter->begin( ); iter2 != iter->end( ); iter2++ ){
                    errorValues += std::to_string( *iter2 ) + " ";
                }
                errorValues += "\n";
            }
            errorValues += "\n";
            errorValues += "     currentdMicroGradientGdMicroGradientCohesion:";
            for ( auto iter = currentdMicroGradientGdMicroGradientCohesion.begin( ); iter != currentdMicroGradientGdMicroGradientCohesion.end( ); iter++ ){
                errorValues += "        ";
                for ( auto iter2 = iter->begin( ); iter2 != iter->end( ); iter2++ ){
                    errorValues += std::to_string( *iter2 ) + " ";
                }
                errorValues += "\n";
            }
            errorValues += "\n";
            return new errorNode( __func__, message + errorValues );
        }
        //TODO: End of debugging section

//        std::cout << "\n  current flow directions\n";
//        std::cout << "    "; tardigradeVectorTools::print( currentMacroFlowDirection );
//        std::cout << "    "; tardigradeVectorTools::print( currentMicroFlowDirection );
//        std::cout << "    "; tardigradeVectorTools::print( currentMicroGradientFlowDirection );

        /*!============================
        | Assemble the flow Jacobians |
        =============================*/

        //Assemble the Jacobians w.r.t. the plastic deformation
        variableMatrix dMacroFlowDirectiondPlasticDeformationGradient
            = tardigradeVectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dMacroFlowDirectiondElasticRightCauchyGreen, dElasticRightCauchyGreendPlasticDeformationGradient );

        variableMatrix dMacroFlowDirectiondPlasticMicroDeformation
            = tardigradeVectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdPlasticMicroDeformation );

        variableMatrix dMacroFlowDirectiondPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdPlasticGradientMicroDeformation );

        variableMatrix dMicroFlowDirectiondPlasticDeformationGradient
            = tardigradeVectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dMicroFlowDirectiondElasticRightCauchyGreen, dElasticRightCauchyGreendPlasticDeformationGradient );

        variableMatrix dMicroFlowDirectiondPlasticMicroDeformation
            = tardigradeVectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdPlasticMicroDeformation );

        variableMatrix dMicroFlowDirectiondPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdPlasticGradientMicroDeformation );

        variableMatrix dMicroGradientFlowDirectiondPlasticDeformationGradient
            = tardigradeVectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                dReferenceHigherOrderStressdPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dMicroGradientFlowDirectiondElasticRightCauchyGreen,
                                dElasticRightCauchyGreendPlasticDeformationGradient );

        variableMatrix dMicroGradientFlowDirectiondPlasticMicroDeformation
            = tardigradeVectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                dReferenceHigherOrderStressdPlasticMicroDeformation );

        variableMatrix dMicroGradientFlowDirectiondPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                dReferenceHigherOrderStressdPlasticGradientMicroDeformation );

        //Assemble the Jacobians w.r.t. the fundamental deformation measures
        variableMatrix dMacroFlowDirectiondDeformationGradient, dMacroFlowDirectiondMicroDeformation,
                       dMacroFlowDirectiondGradientMicroDeformation,
                       dMicroFlowDirectiondDeformationGradient, dMicroFlowDirectiondMicroDeformation,
                       dMicroFlowDirectiondGradientMicroDeformation,
                       dMicroGradientFlowDirectiondDeformationGradient, dMicroGradientFlowDirectiondMicroDeformation,
                       dMicroGradientFlowDirectiondGradientMicroDeformation;

        if ( evaluateFullDerivatives ){
            dMacroFlowDirectiondDeformationGradient
                = tardigradeVectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdDeformationGradient )
                + tardigradeVectorTools::dot( dMacroFlowDirectiondElasticRightCauchyGreen, dElasticRightCauchyGreendDeformationGradient );

            dMacroFlowDirectiondMicroDeformation
                = tardigradeVectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdMicroDeformation );

            dMacroFlowDirectiondGradientMicroDeformation
                = tardigradeVectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdGradientMicroDeformation );

            dMicroFlowDirectiondDeformationGradient
                = tardigradeVectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdDeformationGradient )
                + tardigradeVectorTools::dot( dMicroFlowDirectiondElasticRightCauchyGreen, dElasticRightCauchyGreendDeformationGradient );

            dMicroFlowDirectiondMicroDeformation
                = tardigradeVectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdMicroDeformation );

            dMicroFlowDirectiondGradientMicroDeformation
                = tardigradeVectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdGradientMicroDeformation );

            dMicroGradientFlowDirectiondDeformationGradient
                = tardigradeVectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                    dReferenceHigherOrderStressdDeformationGradient )
                + tardigradeVectorTools::dot( dMicroGradientFlowDirectiondElasticRightCauchyGreen,
                                    dElasticRightCauchyGreendDeformationGradient );

            dMicroGradientFlowDirectiondMicroDeformation
                = tardigradeVectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                    dReferenceHigherOrderStressdMicroDeformation );

            dMicroGradientFlowDirectiondGradientMicroDeformation
                = tardigradeVectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                    dReferenceHigherOrderStressdGradientMicroDeformation );
        }

#ifdef DEBUG_MODE

        //Save the flow directions
        DEBUG.emplace( "currentMacroFlowDirection", currentMacroFlowDirection );
        DEBUG.emplace( "currentMicroFlowDirection", currentMicroFlowDirection );
        DEBUG.emplace( "currentMicroGradientFlowDirection", currentMicroGradientFlowDirection );

        //Save the Jacobians of the flow directions
        DEBUG.emplace( "dMacroFlowDirectiondPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dMacroFlowDirectiondPlasticDeformationGradient ) );
        DEBUG.emplace( "dMacroFlowDirectiondPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dMacroFlowDirectiondPlasticMicroDeformation ) );
        DEBUG.emplace( "dMacroFlowDirectiondPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dMacroFlowDirectiondPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dMicroFlowDirectiondPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dMicroFlowDirectiondPlasticDeformationGradient ) );
        DEBUG.emplace( "dMicroFlowDirectiondPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dMicroFlowDirectiondPlasticMicroDeformation ) );
        DEBUG.emplace( "dMicroFlowDirectiondPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dMicroFlowDirectiondPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dMicroGradientFlowDirectiondPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dMicroGradientFlowDirectiondPlasticDeformationGradient ) );
        DEBUG.emplace( "dMicroGradientFlowDirectiondPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dMicroGradientFlowDirectiondPlasticMicroDeformation ) );
        DEBUG.emplace( "dMicroGradientFlowDirectiondPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dMicroGradientFlowDirectiondPlasticGradientMicroDeformation  ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dMacroFlowDirectiondDeformationGradient",
                            tardigradeVectorTools::appendVectors( dMacroFlowDirectiondDeformationGradient ) );
            DEBUG.emplace( "dMacroFlowDirectiondMicroDeformation",
                            tardigradeVectorTools::appendVectors( dMacroFlowDirectiondMicroDeformation ) );
            DEBUG.emplace( "dMacroFlowDirectiondGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dMacroFlowDirectiondGradientMicroDeformation ) );
            DEBUG.emplace( "dMicroFlowDirectiondDeformationGradient",
                            tardigradeVectorTools::appendVectors( dMicroFlowDirectiondDeformationGradient ) );
            DEBUG.emplace( "dMicroFlowDirectiondMicroDeformation",
                            tardigradeVectorTools::appendVectors( dMicroFlowDirectiondMicroDeformation ) );
            DEBUG.emplace( "dMicroFlowDirectiondGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dMicroFlowDirectiondGradientMicroDeformation ) );
            DEBUG.emplace( "dMicroGradientFlowDirectiondDeformationGradient",
                            tardigradeVectorTools::appendVectors( dMicroGradientFlowDirectiondDeformationGradient ) );
            DEBUG.emplace( "dMicroGradientFlowDirectiondMicroDeformation",
                            tardigradeVectorTools::appendVectors( dMicroGradientFlowDirectiondMicroDeformation ) );
            DEBUG.emplace( "dMicroGradientFlowDirectiondGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dMicroGradientFlowDirectiondGradientMicroDeformation ) );
        }

#endif

        /*!=======================================
        | Compute the plastic velocity gradients |
        ========================================*/

        variableVector currentPlasticMacroVelocityGradient, currentPlasticMicroVelocityGradient,
                       currentPlasticMicroGradientVelocityGradient;

        variableVector dPlasticMacroVelocityGradientdMacroGamma, dPlasticMacroVelocityGradientdMicroGamma,
                       dPlasticMicroVelocityGradientdMicroGamma, dPlasticMicroGradientVelocityGradientdMicroGamma;

        variableMatrix dPlasticMicroGradientVelocityGradientdMicroGradientGamma,
                       dPlasticMacroVelocityGradientdElasticRightCauchyGreen,
                       dPlasticMacroVelocityGradientdMacroFlowDirection,
                       dPlasticMacroVelocityGradientdMicroFlowDirection,
                       dPlasticMicroVelocityGradientdElasticMicroRightCauchyGreen,
                       dPlasticMicroVelocityGradientdElasticPsi,
                       dPlasticMicroVelocityGradientdMicroFlowDirection,
                       dPlasticMicroGradientVelocityGradientdElasticMicroRightCauchyGreen,
                       dPlasticMicroGradientVelocityGradientdElasticPsi,
                       dPlasticMicroGradientVelocityGradientdElasticGamma,
                       dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                       dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection;

        error = computePlasticVelocityGradients( *currentMacroGamma, *currentMicroGamma, *currentMicroGradientGamma, 
                                                  currentElasticRightCauchyGreen, currentElasticMicroRightCauchyGreen,
                                                  currentElasticPsi, currentElasticGamma, currentMacroFlowDirection,
                                                  currentMicroFlowDirection, currentMicroGradientFlowDirection,
                                                  currentPlasticMacroVelocityGradient, currentPlasticMicroVelocityGradient,
                                                  currentPlasticMicroGradientVelocityGradient,
                                                  dPlasticMacroVelocityGradientdMacroGamma,
                                                  dPlasticMacroVelocityGradientdMicroGamma,
                                                  dPlasticMicroVelocityGradientdMicroGamma,
                                                  dPlasticMicroGradientVelocityGradientdMicroGamma,
                                                  dPlasticMicroGradientVelocityGradientdMicroGradientGamma,
                                                  dPlasticMacroVelocityGradientdElasticRightCauchyGreen,
                                                  dPlasticMacroVelocityGradientdMacroFlowDirection,
                                                  dPlasticMacroVelocityGradientdMicroFlowDirection,
                                                  dPlasticMicroVelocityGradientdElasticMicroRightCauchyGreen,
                                                  dPlasticMicroVelocityGradientdElasticPsi,
                                                  dPlasticMicroVelocityGradientdMicroFlowDirection,
                                                  dPlasticMicroGradientVelocityGradientdElasticMicroRightCauchyGreen,
                                                  dPlasticMicroGradientVelocityGradientdElasticPsi,
                                                  dPlasticMicroGradientVelocityGradientdElasticGamma,
                                                  dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                                  dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection );
        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the plastic velocity gradients" );
            result->addNext( error );
            return result;
        }

        /*!=========================================================
        | Assemble the jacobians of the plastic velocity gradients |
        ==========================================================*/

        //Compute the Jacobians w.r.t. the plastic deformation
        variableMatrix dPlasticMacroVelocityGradientdPlasticDeformationGradient
            = tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdElasticRightCauchyGreen,
                                dElasticRightCauchyGreendPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                dMacroFlowDirectiondPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticDeformationGradient );

        variableMatrix dPlasticMacroVelocityGradientdPlasticMicroDeformation
            = tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                dMacroFlowDirectiondPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticMicroDeformation );

        variableMatrix dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                dMacroFlowDirectiondPlasticGradientMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticGradientMicroDeformation );

        variableMatrix dPlasticMicroVelocityGradientdPlasticDeformationGradient
            = tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdElasticPsi,
                                dElasticPsidPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticDeformationGradient );

        variableMatrix dPlasticMicroVelocityGradientdPlasticMicroDeformation
            = tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdElasticMicroRightCauchyGreen,
                                dElasticMicroRightCauchyGreendPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdElasticPsi,
                                dElasticPsidPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticMicroDeformation );

        variableMatrix dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticGradientMicroDeformation );

        variableMatrix dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient
            = tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticPsi,
                                dElasticPsidPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                dElasticGammadPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                dMicroGradientFlowDirectiondPlasticDeformationGradient );

        variableMatrix dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation
            = tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticMicroRightCauchyGreen,
                                dElasticMicroRightCauchyGreendPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticPsi,
                                dElasticPsidPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                dElasticGammadPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                dMicroGradientFlowDirectiondPlasticMicroDeformation );

        variableMatrix dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                dElasticGammadPlasticGradientMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticGradientMicroDeformation )
            + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                dMicroGradientFlowDirectiondPlasticGradientMicroDeformation );

        variableMatrix dPlasticMacroVelocityGradientdDeformationGradient, dPlasticMacroVelocityGradientdMicroDeformation,
                       dPlasticMacroVelocityGradientdGradientMicroDeformation,
                       dPlasticMicroVelocityGradientdDeformationGradient, dPlasticMicroVelocityGradientdMicroDeformation,
                       dPlasticMicroVelocityGradientdGradientMicroDeformation,
                       dPlasticGradientMicroVelocityGradientdDeformationGradient, dPlasticGradientMicroVelocityGradientdMicroDeformation,
                       dPlasticGradientMicroVelocityGradientdGradientMicroDeformation;

        if ( evaluateFullDerivatives ){
            dPlasticMacroVelocityGradientdDeformationGradient
                = tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdElasticRightCauchyGreen,
                                    dElasticRightCauchyGreendDeformationGradient )
                + tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                    dMacroFlowDirectiondDeformationGradient )
                + tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondDeformationGradient );
    
            dPlasticMacroVelocityGradientdMicroDeformation
                = tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                    dMacroFlowDirectiondMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondMicroDeformation );
    
            dPlasticMacroVelocityGradientdGradientMicroDeformation
                = tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                    dMacroFlowDirectiondGradientMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondGradientMicroDeformation );
    
            dPlasticMicroVelocityGradientdDeformationGradient
                = tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdElasticPsi,
                                    dElasticPsidDeformationGradient )
                + tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondDeformationGradient );
    
            dPlasticMicroVelocityGradientdMicroDeformation
                = tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdElasticMicroRightCauchyGreen,
                                    dElasticMicroRightCauchyGreendMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdElasticPsi,
                                    dElasticPsidMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondMicroDeformation );
    
            dPlasticMicroVelocityGradientdGradientMicroDeformation
                = tardigradeVectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondGradientMicroDeformation );
    
            dPlasticGradientMicroVelocityGradientdDeformationGradient
                = tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticPsi,
                                    dElasticPsidDeformationGradient )
                + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                    dElasticGammadDeformationGradient )
                + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondDeformationGradient )
                + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                    dMicroGradientFlowDirectiondDeformationGradient );
    
            dPlasticGradientMicroVelocityGradientdMicroDeformation
                = tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticMicroRightCauchyGreen,
                                    dElasticMicroRightCauchyGreendMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticPsi,
                                    dElasticPsidMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                    dElasticGammadMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                    dMicroGradientFlowDirectiondMicroDeformation );
    
            dPlasticGradientMicroVelocityGradientdGradientMicroDeformation
                = tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                    dElasticGammadGradientMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondGradientMicroDeformation )
                + tardigradeVectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                    dMicroGradientFlowDirectiondGradientMicroDeformation );
        }


#ifdef DEBUG_MODE

        //Save the velocity gradients
        DEBUG.emplace( "currentPlasticMacroVelocityGradient", currentPlasticMacroVelocityGradient );
        DEBUG.emplace( "currentPlasticMicroVelocityGradient", currentPlasticMicroVelocityGradient );
        DEBUG.emplace( "currentPlasticMicroGradientVelocityGradient", currentPlasticMicroGradientVelocityGradient );

        //Save the Jacobians w.r.t. the plastic multipliers
        DEBUG.emplace( "dPlasticMacroVelocityGradientdMacroGamma",
                       dPlasticMacroVelocityGradientdMacroGamma );
        DEBUG.emplace( "dPlasticMacroVelocityGradientdMicroGamma",
                       dPlasticMacroVelocityGradientdMicroGamma );
        DEBUG.emplace( "dPlasticMicroVelocityGradientdMicroGamma",
                       dPlasticMicroVelocityGradientdMicroGamma );
        DEBUG.emplace( "dPlasticMicroGradientVelocityGradientdMicroGamma",
                       dPlasticMicroGradientVelocityGradientdMicroGamma );
        DEBUG.emplace( "dPlasticMicroGradientVelocityGradientdMicroGradientGamma",
                       tardigradeVectorTools::appendVectors( dPlasticMicroGradientVelocityGradientdMicroGradientGamma ) );

        //Save the Jacobians of the velocity gradients
        DEBUG.emplace( "dPlasticMacroVelocityGradientdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dPlasticMacroVelocityGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dPlasticMacroVelocityGradientdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dPlasticMacroVelocityGradientdPlasticMicroDeformation ) );
        DEBUG.emplace( "dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dPlasticMicroVelocityGradientdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dPlasticMicroVelocityGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dPlasticMicroVelocityGradientdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dPlasticMicroVelocityGradientdPlasticMicroDeformation ) );
        DEBUG.emplace( "dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation ) );
        DEBUG.emplace( "dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dPlasticMacroVelocityGradientdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dPlasticMacroVelocityGradientdDeformationGradient ) );
            DEBUG.emplace( "dPlasticMacroVelocityGradientdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticMacroVelocityGradientdMicroDeformation ) );
            DEBUG.emplace( "dPlasticMacroVelocityGradientdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticMacroVelocityGradientdGradientMicroDeformation ) );
            DEBUG.emplace( "dPlasticMicroVelocityGradientdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dPlasticMicroVelocityGradientdDeformationGradient ) );
            DEBUG.emplace( "dPlasticMicroVelocityGradientdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticMicroVelocityGradientdMicroDeformation ) );
            DEBUG.emplace( "dPlasticMicroVelocityGradientdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticMicroVelocityGradientdGradientMicroDeformation ) );
            DEBUG.emplace( "dPlasticMicroGradientVelocityGradientdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdDeformationGradient ) );
            DEBUG.emplace( "dPlasticMicroGradientVelocityGradientdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdMicroDeformation ) );
            DEBUG.emplace( "dPlasticMicroGradientVelocityGradientdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdGradientMicroDeformation ) );
        }

#endif

        /*!===============================
        | Evolve the plastic deformation |
        ================================*/

        variableVector expectedPlasticDeformationGradient, expectedPlasticMicroDeformation, expectedPlasticGradientMicroDeformation;

        variableMatrix dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                       dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                       dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                       dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                       dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient;

        error = evolvePlasticDeformation( *Dt,
                                           currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient,
                                           currentPlasticMicroGradientVelocityGradient,
                                          *previousPlasticDeformationGradient,
                                          *previousPlasticMicroDeformation,
                                          *previousPlasticGradientMicroDeformation,
                                          *previousPlasticMacroVelocityGradient,
                                          *previousPlasticMicroVelocityGradient,
                                          *previousPlasticMicroGradientVelocityGradient,
                                           expectedPlasticDeformationGradient,
                                           expectedPlasticMicroDeformation,
                                           expectedPlasticGradientMicroDeformation,
                                           dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                           dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                           dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                           dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                           dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                          *alphaMacro, *alphaMicro, *alphaMicroGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the plastic deformation" );
            result->addNext( error );
            return result;
        }

        /*!==========================================
        | Assemble the plastic deformation Jacobian |
        ===========================================*/

        //Compute the Jacobians w.r.t. the plastic multipliers
        variableVector dExpectedPlasticDeformationGradientdMacroGamma
            = tardigradeVectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdMacroGamma );

        variableVector dExpectedPlasticDeformationGradientdMicroGamma
            = tardigradeVectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdMicroGamma );

        variableVector dExpectedPlasticMicroDeformationdMicroGamma
            = tardigradeVectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdMicroGamma );

        variableVector dExpectedPlasticGradientMicroDeformationdMacroGamma
            = tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdMacroGamma );

        variableVector dExpectedPlasticGradientMicroDeformationdMicroGamma
            = tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdMicroGamma )
            + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdMicroGamma )
            + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                dPlasticMicroGradientVelocityGradientdMicroGamma );

        variableMatrix dExpectedPlasticGradientMicroDeformationdMicroGradientGamma
            = tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                dPlasticMicroGradientVelocityGradientdMicroGradientGamma );

        //Compute the Jacobians w.r.t. the plastic deformation
        variableMatrix dExpectedPlasticDeformationGradientdPlasticDeformationGradient
            = tardigradeVectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticDeformationGradient );

        variableMatrix dExpectedPlasticDeformationGradientdPlasticMicroDeformation
            = tardigradeVectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticMicroDeformation );

        variableMatrix dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation );

        variableMatrix dExpectedPlasticMicroDeformationdPlasticDeformationGradient
            = tardigradeVectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticDeformationGradient );

        variableMatrix dExpectedPlasticMicroDeformationdPlasticMicroDeformation
            = tardigradeVectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticMicroDeformation );

        variableMatrix dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation );

        variableMatrix dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient
            = tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticDeformationGradient )
            + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient );

        variableMatrix dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation
            = tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticMicroDeformation )
            + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation );

        variableMatrix dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation
            = tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation )
            + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation )
            + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation );

        variableMatrix dExpectedPlasticDeformationGradientdDeformationGradient, dExpectedPlasticDeformationGradientdMicroDeformation,
                       dExpectedPlasticDeformationGradientdGradientMicroDeformation,
                       dExpectedPlasticMicroDeformationdDeformationGradient, dExpectedPlasticMicroDeformationdMicroDeformation,
                       dExpectedPlasticMicroDeformationdGradientMicroDeformation,
                       dExpectedPlasticGradientMicroDeformationdDeformationGradient,
                       dExpectedPlasticGradientMicroDeformationdMicroDeformation,
                       dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation;

        if ( evaluateFullDerivatives ){
            dExpectedPlasticDeformationGradientdDeformationGradient
                = tardigradeVectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdDeformationGradient );
    
            dExpectedPlasticDeformationGradientdMicroDeformation
                = tardigradeVectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdMicroDeformation );
    
            dExpectedPlasticDeformationGradientdGradientMicroDeformation
                = tardigradeVectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdGradientMicroDeformation );
    
            dExpectedPlasticMicroDeformationdDeformationGradient
                = tardigradeVectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdDeformationGradient );
    
            dExpectedPlasticMicroDeformationdMicroDeformation
                = tardigradeVectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdMicroDeformation );
    
            dExpectedPlasticMicroDeformationdGradientMicroDeformation
                = tardigradeVectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdGradientMicroDeformation );
    
            dExpectedPlasticGradientMicroDeformationdDeformationGradient
                = tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdDeformationGradient )
                + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdDeformationGradient )
                + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                    dPlasticGradientMicroVelocityGradientdDeformationGradient );
    
            dExpectedPlasticGradientMicroDeformationdMicroDeformation
                = tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdMicroDeformation )
                + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdMicroDeformation )
                + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                    dPlasticGradientMicroVelocityGradientdMicroDeformation );
    
            dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation
                = tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdGradientMicroDeformation )
                + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdGradientMicroDeformation )
                + tardigradeVectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                    dPlasticGradientMicroVelocityGradientdGradientMicroDeformation );
        }

#ifdef DEBUG_MODE

        //Save the velocity gradients
        DEBUG.emplace( "expectedPlasticDeformationGradient", expectedPlasticDeformationGradient );
        DEBUG.emplace( "expectedPlasticMicroDeformation", expectedPlasticMicroDeformation );
        DEBUG.emplace( "expectedPlasticGradientMicroDeformation", expectedPlasticGradientMicroDeformation );

        //Save the Jacobians w.r.t. the plastic multipliers
        DEBUG.emplace( "dExpectedPlasticDeformationGradientdMacroGamma",
                        dExpectedPlasticDeformationGradientdMacroGamma );
        DEBUG.emplace( "dExpectedPlasticDeformationGradientdMicroGamma",
                        dExpectedPlasticDeformationGradientdMicroGamma );
        DEBUG.emplace( "dExpectedPlasticMicroDeformationdMicroGamma",
                        dExpectedPlasticMicroDeformationdMicroGamma );
        DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdMacroGamma",
                        dExpectedPlasticGradientMicroDeformationdMacroGamma );
        DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdMicroGamma",
                        dExpectedPlasticGradientMicroDeformationdMicroGamma );
        DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdMicroGradientGamma",
                        tardigradeVectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdMicroGradientGamma ) );

        //Save the Jacobians of the velocity gradients
        DEBUG.emplace( "dExpectedPlasticDeformationGradientdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dExpectedPlasticDeformationGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dExpectedPlasticDeformationGradientdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dExpectedPlasticDeformationGradientdPlasticMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticMicroDeformationdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dExpectedPlasticMicroDeformationdPlasticDeformationGradient ) );
        DEBUG.emplace( "dExpectedPlasticMicroDeformationdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dExpectedPlasticMicroDeformationdPlasticMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient",
                       tardigradeVectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient ) );
        DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation",
                       tardigradeVectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation",
                       tardigradeVectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dExpectedPlasticDeformationGradientdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dExpectedPlasticDeformationGradientdDeformationGradient ) );
            DEBUG.emplace( "dExpectedPlasticDeformationGradientdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dExpectedPlasticDeformationGradientdMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticDeformationGradientdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dExpectedPlasticDeformationGradientdGradientMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticMicroDeformationdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dExpectedPlasticMicroDeformationdDeformationGradient ) );
            DEBUG.emplace( "dExpectedPlasticMicroDeformationdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dExpectedPlasticMicroDeformationdMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticMicroDeformationdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dExpectedPlasticMicroDeformationdGradientMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdDeformationGradient ) );
            DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation ) );
        }

#endif

//        /*!=============================
//        | Evaluate the yield equations |
//        ==============================*/
//
//        variableVector yieldFunctionValues;
//        
//        variableType   dMacroYielddMacroCohesion, dMicroYielddMicroCohesion;
//        variableVector dMacroYielddPK2Stress, dMacroYielddElasticRightCauchyGreen;
//        variableVector dMicroYielddReferenceMicroStress, dMicroYielddElasticRightCauchyGreen;
//
//        variableMatrix dMicroGradientYielddReferenceHigherOrderStress, dMicroGradientYielddMicroGradientCohesion,
//                       dMicroGradientYielddElasticRightCauchyGreen;
//
//        error = evaluateYieldFunctions( currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress,
//                                        currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion,
//                                        currentElasticRightCauchyGreen,
//                                       *macroYieldParameters, *microYieldParameters, *microGradientYieldParameters,
//                                        yieldFunctionValues,
//                                        dMacroYielddPK2Stress, dMacroYielddMacroCohesion, dMacroYielddElasticRightCauchyGreen,
//                                        dMicroYielddReferenceMicroStress, dMicroYielddMicroCohesion, dMicroYielddElasticRightCauchyGreen,
//                                        dMicroGradientYielddReferenceHigherOrderStress, dMicroGradientYielddMicroGradientCohesion,
//                                        dMicroGradientYielddElasticRightCauchyGreen
//#ifdef DEBUG_MODE
//                                     , DEBUG
//#endif
//                                    );
//
//        if ( error ){
//            errorOut result = new errorNode( "computePlasticDeformationResidual",
//                                             "Error in the computation of the yield functions\n" );
//            result->addNext( error );
//            return result;
//        }
//
////        std::cout << "\n  expected yield values\n";
////        std::cout << "    "; tardigradeVectorTools::print( yieldFunctionValues );
//
//        /*!===============================================
//        | Construct the Jacobians of the yield equations |
//        ================================================*/
//
//        //Construct the Jacobians w.r.t. the plastic deformation measures
//        variableVector dMacroYielddPlasticDeformationGradient
//            = tardigradeVectorTools::Tdot( dPK2StressdPlasticDeformationGradient, dMacroYielddPK2Stress )
//            + tardigradeVectorTools::Tdot( dElasticRightCauchyGreendPlasticDeformationGradient, dMacroYielddElasticRightCauchyGreen );
//
//        variableVector dMacroYielddPlasticMicroDeformation
//            = tardigradeVectorTools::Tdot( dPK2StressdPlasticMicroDeformation, dMacroYielddPK2Stress );
//
//        variableVector dMacroYielddPlasticGradientMicroDeformation
//            = tardigradeVectorTools::Tdot( dPK2StressdPlasticGradientMicroDeformation, dMacroYielddPK2Stress );
//
//        variableVector dMicroYielddPlasticDeformationGradient
//            = tardigradeVectorTools::Tdot( dReferenceMicroStressdPlasticDeformationGradient, dMicroYielddReferenceMicroStress )
//            + tardigradeVectorTools::Tdot( dElasticRightCauchyGreendPlasticDeformationGradient, dMicroYielddElasticRightCauchyGreen );
//
//        variableVector dMicroYielddPlasticMicroDeformation
//            = tardigradeVectorTools::Tdot( dReferenceMicroStressdPlasticMicroDeformation, dMicroYielddReferenceMicroStress );
//
//        variableVector dMicroYielddPlasticGradientMicroDeformation
//            = tardigradeVectorTools::Tdot( dReferenceMicroStressdPlasticGradientMicroDeformation, dMicroYielddReferenceMicroStress );
//
//        variableMatrix dMicroGradientYielddPlasticDeformationGradient
//            = tardigradeVectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress, dReferenceHigherOrderStressdPlasticDeformationGradient )
//            + tardigradeVectorTools::dot( dMicroGradientYielddElasticRightCauchyGreen, dElasticRightCauchyGreendPlasticDeformationGradient );
//
//        variableMatrix dMicroGradientYielddPlasticMicroDeformation
//            = tardigradeVectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress, dReferenceHigherOrderStressdPlasticMicroDeformation );
//
//        variableMatrix dMicroGradientYielddPlasticGradientMicroDeformation
//            = tardigradeVectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress,
//                                dReferenceHigherOrderStressdPlasticGradientMicroDeformation );
//
//        //Construct the Jacobians w.r.t. the plastic multipliers
//        variableType dMacroYielddMacroGamma = dMacroYielddMacroCohesion * dCurrentMacroCohesiondMacroGamma;
//        variableType dMicroYielddMicroGamma = dMicroYielddMicroCohesion * dCurrentMicroCohesiondMicroGamma;
//        variableMatrix dMicroGradientYielddMicroGradientGamma = tardigradeVectorTools::dot( dMicroGradientYielddMicroGradientCohesion,
//                                                                                  dCurrentMicroGradientCohesiondMicroGradientGamma );
//
//        variableVector dMacroYielddDeformationGradient, dMacroYielddMicroDeformation, dMacroYielddGradientMicroDeformation,
//                       dMicroYielddDeformationGradient, dMicroYielddMicroDeformation, dMicroYielddGradientMicroDeformation;
//        variableMatrix dMicroGradientYielddDeformationGradient, dMicroGradientYielddMicroDeformation,
//                       dMicroGradientYielddGradientMicroDeformation;
//
//        if ( evaluateFullDerivatives ){
//            dMacroYielddDeformationGradient
//                = tardigradeVectorTools::Tdot( dPK2StressdDeformationGradient, dMacroYielddPK2Stress )
//                + tardigradeVectorTools::Tdot( dElasticRightCauchyGreendDeformationGradient, dMacroYielddElasticRightCauchyGreen );
//    
//            dMacroYielddMicroDeformation
//                = tardigradeVectorTools::Tdot( dPK2StressdMicroDeformation, dMacroYielddPK2Stress );
//    
//            dMacroYielddGradientMicroDeformation
//                = tardigradeVectorTools::Tdot( dPK2StressdGradientMicroDeformation, dMacroYielddPK2Stress );
//    
//            dMicroYielddDeformationGradient
//                = tardigradeVectorTools::Tdot( dReferenceMicroStressdDeformationGradient, dMicroYielddReferenceMicroStress )
//                + tardigradeVectorTools::Tdot( dElasticRightCauchyGreendDeformationGradient, dMicroYielddElasticRightCauchyGreen );
//    
//            dMicroYielddMicroDeformation
//                = tardigradeVectorTools::Tdot( dReferenceMicroStressdMicroDeformation, dMicroYielddReferenceMicroStress );
//    
//            dMicroYielddGradientMicroDeformation
//                = tardigradeVectorTools::Tdot( dReferenceMicroStressdGradientMicroDeformation, dMicroYielddReferenceMicroStress );
//    
//            dMicroGradientYielddDeformationGradient
//                = tardigradeVectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress, dReferenceHigherOrderStressdDeformationGradient )
//                + tardigradeVectorTools::dot( dMicroGradientYielddElasticRightCauchyGreen, dElasticRightCauchyGreendDeformationGradient );
//    
//            dMicroGradientYielddMicroDeformation
//                = tardigradeVectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress, dReferenceHigherOrderStressdMicroDeformation );
//    
//            dMicroGradientYielddGradientMicroDeformation
//                = tardigradeVectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress,
//                                    dReferenceHigherOrderStressdGradientMicroDeformation );
//        }
//
//#ifdef DEBUG_MODE
//
//        //Save the values to the debug map
//        temp = { yieldFunctionValues[ 0 ] };
//        DEBUG.emplace( "macroYieldFunction", temp );
// 
//        temp = { yieldFunctionValues[ 1 ] };
//        DEBUG.emplace( "microYieldFunction", temp );
// 
//        DEBUG.emplace( "microGradientYieldFunction",
//                        variableVector( yieldFunctionValues.begin() + 2, yieldFunctionValues.begin() + 5 ) );
// 
//        //Save the jacobians w.r.t. the plastic deformation
//        DEBUG.emplace( "dMacroYielddPlasticDeformationGradient", dMacroYielddPlasticDeformationGradient );
//        DEBUG.emplace( "dMacroYielddPlasticMicroDeformation", dMacroYielddPlasticMicroDeformation );
//        DEBUG.emplace( "dMacroYielddPlasticGradientMicroDeformation", dMacroYielddPlasticGradientMicroDeformation );
// 
//        DEBUG.emplace( "dMicroYielddPlasticDeformationGradient", dMicroYielddPlasticDeformationGradient );
//        DEBUG.emplace( "dMicroYielddPlasticMicroDeformation", dMicroYielddPlasticMicroDeformation );
//        DEBUG.emplace( "dMicroYielddPlasticGradientMicroDeformation", dMicroYielddPlasticGradientMicroDeformation );
// 
//        DEBUG.emplace( "dMicroGradientYielddPlasticDeformationGradient",
//                        tardigradeVectorTools::appendVectors( dMicroGradientYielddPlasticDeformationGradient ) );
//        DEBUG.emplace( "dMicroGradientYielddPlasticMicroDeformation",
//                        tardigradeVectorTools::appendVectors( dMicroGradientYielddPlasticMicroDeformation ) );
//        DEBUG.emplace( "dMicroGradientYielddPlasticGradientMicroDeformation",
//                        tardigradeVectorTools::appendVectors( dMicroGradientYielddPlasticGradientMicroDeformation ) );
// 
//        //Save the Jacobians w.r.t. the strain-like ISVs
//        temp = { dMacroYielddMacroGamma };
//        DEBUG.emplace( "dMacroYielddMacroGamma", temp );
//        temp = { dMicroYielddMicroGamma };
//        DEBUG.emplace( "dMicroYielddMicroGamma", temp );
//        DEBUG.emplace( "dMicroGradientYielddMicroGradientGamma",
//                        tardigradeVectorTools::appendVectors( dMicroGradientYielddMicroGradientGamma ) );
//
//        if ( evaluateFullDerivatives ){
//            DEBUG.emplace( "dMacroYielddDeformationGradient", dMacroYielddDeformationGradient );
//            DEBUG.emplace( "dMacroYielddMicroDeformation", dMacroYielddMicroDeformation );
//            DEBUG.emplace( "dMacroYielddGradientMicroDeformation", dMacroYielddGradientMicroDeformation );
//            DEBUG.emplace( "dMicroYielddDeformationGradient", dMicroYielddDeformationGradient );
//            DEBUG.emplace( "dMicroYielddMicroDeformation", dMicroYielddMicroDeformation );
//            DEBUG.emplace( "dMicroYielddGradientMicroDeformation", dMicroYielddGradientMicroDeformation );
//            DEBUG.emplace( "dMicroGradientYielddDeformationGradient",
//                            tardigradeVectorTools::appendVectors( dMicroGradientYielddDeformationGradient ) );
//            DEBUG.emplace( "dMicroGradientYielddMicroDeformation",
//                            tardigradeVectorTools::appendVectors( dMicroGradientYielddMicroDeformation ) );
//            DEBUG.emplace( "dMicroGradientYielddGradientMicroDeformation",
//                            tardigradeVectorTools::appendVectors( dMicroGradientYielddGradientMicroDeformation ) );
//        }
//
//#endif

        
        /*!===============================================
        | Compute the residual equation and the Jacobian |
        ================================================*/

        residual = tardigradeSolverTools::floatVector( x.size(), 0 );
        jacobian = tardigradeVectorTools::eye< tardigradeSolverTools::floatType >( x.size() );

        if ( evaluateFullDerivatives ){
            floatOuts[ 13 ] = tardigradeSolverTools::floatVector( x.size() * 45, 0 );
            floatOuts[ 15 ] = tardigradeSolverTools::floatVector( x.size() * 5, 0 );
        }

        //Compute the residuals and the jacobians for the plastic deformations
        for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
            residual[ i ] = currentPlasticDeformationGradient[ i ] - expectedPlasticDeformationGradient[ i ];

            //The plastic deformation residuals
            for ( unsigned int j = 0; j < currentPlasticDeformationGradient.size(); j++ ){
                jacobian[ i ][ j ] -= dExpectedPlasticDeformationGradientdPlasticDeformationGradient[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 13 ][ 45 * i + j ] -= dExpectedPlasticDeformationGradientdDeformationGradient[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticMicroDeformation.size(); j++ ){
                jacobian[ i ][ j + 9 ] -= dExpectedPlasticDeformationGradientdPlasticMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 13 ][ 45 * i + j + 9 ] -= dExpectedPlasticDeformationGradientdMicroDeformation[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticGradientMicroDeformation.size(); j++ ){
                jacobian[ i ][ j + 18 ] -= dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 13 ][ 45 * i + j + 18 ] -= dExpectedPlasticDeformationGradientdGradientMicroDeformation[ i ][ j ];
                }
            }
            
            //The plastic multiplier residuals
            if ( evaluateFullDerivatives ){
                floatOuts[ 15 ][ 5 * i + 0 ] -= dExpectedPlasticDeformationGradientdMacroGamma[ i ];
                floatOuts[ 15 ][ 5 * i + 1 ] -= dExpectedPlasticDeformationGradientdMicroGamma[ i ];
            }
        }

        for ( unsigned int i = 0; i < currentPlasticMicroDeformation.size(); i++ ){
            residual[ i + 9 ] = currentPlasticMicroDeformation[ i ] - expectedPlasticMicroDeformation[ i ];

            for ( unsigned int j = 0; j < currentPlasticDeformationGradient.size(); j++ ){
                jacobian[ i + 9 ][ j ] -= dExpectedPlasticMicroDeformationdPlasticDeformationGradient[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 13 ][ 45 * ( i + 9 ) + j ] -= dExpectedPlasticMicroDeformationdDeformationGradient[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticMicroDeformation.size(); j++ ){
                jacobian[ i + 9 ][ j + 9 ] -= dExpectedPlasticMicroDeformationdPlasticMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 13 ][ 45 * ( i + 9 ) + j + 9 ] -= dExpectedPlasticMicroDeformationdMicroDeformation[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticGradientMicroDeformation.size(); j++ ){
                jacobian[ i + 9 ][ j + 18 ] -= dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 13 ][ 45 * ( i + 9 ) + j + 18 ] -= dExpectedPlasticMicroDeformationdGradientMicroDeformation[ i ][ j ];
                }
            }
            
            //The plastic multiplier residuals
            if ( evaluateFullDerivatives ){
                floatOuts[ 15 ][ 5 * ( i + 9 ) + 1 ] -= dExpectedPlasticMicroDeformationdMicroGamma[ i ];
            }
        }

        for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
            residual[ i + 18 ] = currentPlasticGradientMicroDeformation[ i ] - expectedPlasticGradientMicroDeformation[ i ];

            for ( unsigned int j = 0; j < currentPlasticDeformationGradient.size(); j++ ){
                jacobian[ i + 18 ][ j ] -= dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 13 ][ 45 * ( i + 18 ) + j ] -= dExpectedPlasticGradientMicroDeformationdDeformationGradient[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticMicroDeformation.size(); j++ ){
                jacobian[ i + 18 ][ j + 9 ] -= dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 13 ][ 45 * ( i + 18 ) + j + 9 ] -= dExpectedPlasticGradientMicroDeformationdMicroDeformation[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticGradientMicroDeformation.size(); j++ ){
                jacobian[ i + 18 ][ j + 18 ] -= dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 13 ][ 45 * ( i + 18 ) + j + 18 ] -= dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation[ i ][ j ];
                }
            }

            //The plastic multiplier residuals
            if ( evaluateFullDerivatives ){
                floatOuts[ 15 ][ 5 * ( i + 18 ) + 0 ] -= dExpectedPlasticGradientMicroDeformationdMacroGamma[ i ];
                floatOuts[ 15 ][ 5 * ( i + 18 ) + 1 ] -= dExpectedPlasticGradientMicroDeformationdMicroGamma[ i ]; 
                for ( unsigned int j = 0; j < currentMicroGradientGamma->size(); j++ ){
                    floatOuts[ 15 ][ 5 * ( i + 18 ) + 2 + j ] -= dExpectedPlasticGradientMicroDeformationdMicroGradientGamma[ i ][ j ];
                }
            }
        }

//        //Compute the residuals and the Jacobians for the Kuhn-Tucker and yield conditions
//
//        //Set the residuals and Jacobians for the macro-scale terms
//        variableType dMacMacroYielddMacroYield;
//        residual[ 45 ] = currentMacroGamma * yieldFunctionValues[ 0 ]
//                       + tardigradeConstitutiveTools::mac( yieldFunctionValues[ 0 ], dMacMacroYielddMacroYield );
//
//        //The Jacobian terms w.r.t. the plastic deformation
//        for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
//            jacobian[ 45 ][ i ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddPlasticDeformationGradient[ i ];
//
//            if ( evaluateFullDerivatives ){
//                floatOuts[ 8 ][ 45 * 45 + i ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddDeformationGradient[ i ];
//            }
//        }
//
//        for ( unsigned int i = 0; i < currentPlasticMicroDeformation.size(); i++ ){
//            jacobian[ 45 ][ i + 9 ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddPlasticMicroDeformation[ i ];
//
//            if ( evaluateFullDerivatives ){
//                floatOuts[ 8 ][ 45 * 45 + i + 9 ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddMicroDeformation[ i ];
//            }
//        }
//
//        for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
//            jacobian[ 45 ][ i + 18 ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddPlasticGradientMicroDeformation[ i ];
//
//            if ( evaluateFullDerivatives ){
//                floatOuts[ 8 ][ 45 * 45 + i + 18 ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddGradientMicroDeformation[ i ];
//            }
//        }
//
//        //The Jacobian terms w.r.t. the plastic multipliers
//        jacobian[ 45 ][ 45 ] = yieldFunctionValues[ 0 ]
//                             + ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddMacroGamma;
//
//        //Set the residuals and Jacobians for the micro-scale terms
//        variableType dMacMicroYielddMicroYield;
//        residual[ 46 ] = currentMicroGamma * yieldFunctionValues[ 1 ]
//                       + tardigradeConstitutiveTools::mac( yieldFunctionValues[ 1 ], dMacMicroYielddMicroYield );
//
//        //The Jacobian terms w.r.t. the plastic deformation
//        for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
//            jacobian[ 46 ][ i ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddPlasticDeformationGradient[ i ];
//
//            if ( evaluateFullDerivatives ){
//                floatOuts[ 8 ][ 46 * 45 + i ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddDeformationGradient[ i ];
//            }
//        }
//
//        for ( unsigned int i = 0; i < currentPlasticMicroDeformation.size(); i++ ){
//            jacobian[ 46 ][ i + 9 ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddPlasticMicroDeformation[ i ];
//
//            if ( evaluateFullDerivatives ){
//                floatOuts[ 8 ][ 46 * 45 + i + 9 ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddMicroDeformation[ i ];
//            }
//        }
//
//        for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
//            jacobian[ 46 ][ i + 18 ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddPlasticGradientMicroDeformation[ i ];
//
//            if ( evaluateFullDerivatives ){
//                floatOuts[ 8 ][ 46 * 45 + i + 18 ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddGradientMicroDeformation[ i ];
//            }
//        }
//
//        //The Jacobian terms w.r.t. the plastic multipliers
//        jacobian[ 46 ][ 46 ] = yieldFunctionValues[ 1 ]
//                             + ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddMicroGamma;
//
//        //Set the residuals and jacobians for the micro gradient terms
//        variableType dMacMicroGradientYielddMicroGradientYield;
//        for ( unsigned int i = 0; i < 3; i++ ){
//            residual[ 47 + i ] = currentMicroGradientGamma[ i ] * yieldFunctionValues[ 2 + i ]
//                               + tardigradeConstitutiveTools::mac( yieldFunctionValues[ 2 + i ], dMacMicroGradientYielddMicroGradientYield );
//
//            //The Jacobian terms w.r.t. the plastic deformation
//            for ( unsigned int j = 0; j < currentPlasticDeformationGradient.size(); j++ ){
//                jacobian[ 47 + i ][ j ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddPlasticDeformationGradient[ i ][ j ];
//
//                if ( evaluateFullDerivatives ){
//                    floatOuts[ 8 ][ ( 47 + i ) * 45 + j ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddDeformationGradient[ i ][ j ];
//                }
//            }
//
//            for ( unsigned int j = 0; j < currentPlasticMicroDeformation.size(); j++ ){
//                jacobian[ 47 + i ][ j + 9 ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddPlasticMicroDeformation[ i ][ j ];
//
//                if ( evaluateFullDerivatives ){
//                    floatOuts[ 8 ][ ( 47 + i ) * 45 + j + 9 ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddMicroDeformation[ i ][ j ];
//                }
//            }
//
//            for ( unsigned int j = 0; j < currentPlasticGradientMicroDeformation.size(); j++ ){
//                jacobian[ 47 + i ][ j + 18 ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddPlasticGradientMicroDeformation[ i ][ j ];
//
//                if ( evaluateFullDerivatives ){
//                    floatOuts[ 8 ][ ( 47 + i ) * 45 + j + 18 ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddGradientMicroDeformation[ i ][ j ];
//                }
//            }
//
//            //The Jacobian terms w.r.t. the plastic multipliers
//            jacobian[ 47 + i ][ 47 + i ] = yieldFunctionValues[ 2 + i ];
//            for ( unsigned int j = 0; j < 3; j++ ){
//                jacobian[ 47 + i ][ 47 + j ] += ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield )
//                                              * dMicroGradientYielddMicroGradientGamma[ i ][ j ];
//            }
//        }

        //Save the stresses
        floatOuts[ 0 ] = currentPK2Stress;
        floatOuts[ 1 ] = currentReferenceMicroStress;
        floatOuts[ 2 ] = currentReferenceHigherOrderStress;
        floatOuts[ 3 ] = { currentMacroStrainISV };
        floatOuts[ 4 ] = { currentMicroStrainISV };
        floatOuts[ 5 ] = currentMicroGradientStrainISV;
        floatOuts[ 6 ] = { currentMacroCohesion };
        floatOuts[ 7 ] = { currentMicroCohesion };
        floatOuts[ 8 ] = currentMicroGradientCohesion;

//        //Try conditioning the residual and jacobian matrix
//        constantType norm;
//        for ( unsigned int i = 0; i < jacobian.size(); i++ ){
//            norm = tardigradeVectorTools::l2norm( jacobian[ i ] );
//
//            residual[ i ] /= norm;
//            jacobian[ i ] /= norm;
//        }
//
//        //Check the condition number of the Jacobian
//        variableVector jflat = tardigradeVectorTools::appendVectors( jacobian );
//        Eigen::Map < const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor> > A(jflat.data(), 55, 55);
//
//        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
//        double cond = svd.singularValues()(0)
//        / svd.singularValues()(svd.singularValues().size()-1);
//
//        std::cout << "\nJACOBAIN CONDITION NUMBER: " << cond << "\n";
//        if ( cond > 1000 ){
//            std::cout << "\nJACOBIAN SINGULAR VALUES:\n" << svd.singularValues() << "\n";
//    
//            std::cout << "\nRESIDUAL\n"; tardigradeVectorTools::print( residual );
//    
//            std::cout << "\nJACOBIAN\n" << A.block<55,18>(0,0) << "\n";
//            std::cout << A.block<55,10>(0,45) << "\n";
//            assert( 1 == 0 );
//        }

        return NULL;
    }

    int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                        const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                        const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                        const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                        std::vector< double > &SDVS,
                        const std::vector< double > &current_ADD_DOF,
                        const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                        const std::vector< double > &previous_ADD_DOF,
                        const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                        std::vector< double > &current_PK2, std::vector< double > &current_SIGMA, std::vector< double > &current_M,
                        std::vector< std::vector< double > > &ADD_TERMS,
                        std::string &output_message
#ifdef DEBUG_MODE
                        , tardigradeSolverTools::homotopyMap &DEBUG
#endif
                        ){
        /*!
         * Evaluate the elasto-plastic constitutive model. Note the format of the header changed to provide a 
         * consistant interface with the material model library.
         *
         * :param const std::vector< double > &time: The current time and the timestep
         *     [ current_t, dt ]
         * :param const std::vector< double > ( &fparams ): The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters,
         *       num_macroHardeningParameters, macroHardeningParameters,
         *       num_microHardeningParameters, microHardeningParameters,
         *       num_microGradientHardeningParameters, microGradientHardeningParameters,
         *       num_macroFlowParameters, macroFlowParameters,
         *       num_microFlowParameters, microFlowParameters,
         *       num_microGradientFlowParameters, microGradientFlowParameters,
         *       num_macroYieldParameters, macroYieldParameters,
         *       num_microYieldParameters, microYieldParameters,
         *       num_microGradientYieldParameters, microGradientYieldParameters,
         *       alphaMacro, alphaMicro, alphaMicroGradient,
         *       relativeTolerance, absoluteTolerance ]
         *
         * :param const double ( &current_grad_u )[ 3 ][ 3 ]: The current displacement gradient
         *     Assumed to be of the form [ [ u_{1,1}, u_{1,2}, u_{1,3} ],
         *                                 [ u_{2,1}, u_{2,2}, u_{2,3} ],
         *                                 [ u_{3,1}, u_{3,2}, u_{3,3} ] ]
         * :param const double ( &current_phi )[ 9 ]: The current micro displacment values.
         *     Assumed to be of the form [ \phi_{11}, \phi_{12}, \phi_{13}, \phi_{21}, \phi_{22}, \phi_{23}, \phi_{31}, \phi_{32}, \phi_{33} ]
         * :param const double ( &current_grad_phi )[ 9 ][ 3 ]: The current micro displacement gradient
         *     Assumed to be of the form [ [ \phi_{11,1}, \phi_{11,2}, \phi_{11,3} ],
         *                                 [ \phi_{12,1}, \phi_{12,2}, \phi_{12,3} ],
         *                                 [ \phi_{13,1}, \phi_{13,2}, \phi_{13,3} ],
         *                                 [ \phi_{21,1}, \phi_{21,2}, \phi_{21,3} ],
         *                                 [ \phi_{22,1}, \phi_{22,2}, \phi_{22,3} ],
         *                                 [ \phi_{23,1}, \phi_{23,2}, \phi_{23,3} ],
         *                                 [ \phi_{31,1}, \phi_{31,2}, \phi_{31,3} ],
         *                                 [ \phi_{32,1}, \phi_{32,2}, \phi_{32,3} ],
         *                                 [ \phi_{33,1}, \phi_{33,2}, \phi_{33,3} ] ]
         * :param const double ( &previous_grad_u )[ 3 ][ 3 ]: The previous displacement gradient.
         * :param const double ( &previous_phi )[ 9 ]: The previous micro displacement.
         * :param const double ( &previous_grad_phi )[ 9 ][ 3 ]: The previous micro displacement gradient.
         * :param std::vector< double > &SDVS: The previously converged values of the state variables
         *     [ previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
         *       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
         *       previousPlasticDeformationGradient - eye, previousPlasticMicroDeformation - eye,
         *       previousPlasticMicroGradient ]
         * :param std::vector< double > &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * :param std::vector< std::vector< double > > &current_ADD_grad_DOF: The current values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * :param std::vector< double > &current_PK2: The current value of the second Piola Kirchhoff stress tensor. The format is
         *     [ S_{11}, S_{12}, S_{13}, S_{21}, S_{22}, S_{23}, S_{31}, S_{32}, S_{33} ]
         * :param std::vector< double > &current_SIGMA: The current value of the reference micro stress. The format is
         *     [ S_{11}, S_{12}, S_{13}, S_{21}, S_{22}, S_{23}, S_{31}, S_{32}, S_{33} ]
         * :param std::vector< double > &current_M: The current value of the reference higher order stress. The format is
         *     [ M_{111}, M_{112}, M_{113}, M_{121}, M_{122}, M_{123}, M_{131}, M_{132}, M_{133},
         *       M_{211}, M_{212}, M_{213}, M_{221}, M_{222}, M_{223}, M_{231}, M_{232}, M_{233},
         *       M_{311}, M_{312}, M_{313}, M_{321}, M_{322}, M_{323}, M_{331}, M_{332}, M_{333} ]
         * :param std::vector< std::vector< double > > &ADD_TERMS: Additional terms ( unused )
         * :param std::string &output_message: The output message string.
         * :param tardigradeSolverTools::homotopyMap DEBUG: The debugging object ( only available if DEBUG_MODE is defined )
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        //Assume 3D
        unsigned int dim = 3;

#ifdef DEBUG_MODE
        tardigradeSolverTools::debugMap tempDEBUG;
        tardigradeSolverTools::iterationMap tempITERATION;
#endif

        //Construct identity matrix
        constantVector eye( dim * dim, 0 );
        tardigradeVectorTools::eye< constantType >( eye );

        //Re-direct the output to a buffer
        std::stringbuf buffer;
        cerr_redirect rd( &buffer );

        /*=============================
        | Extract the incoming values |
        ==============================*/

        //Extract the time
        if ( time.size() != 2 ){
            errorOut result = new errorNode( "evaluate_model",
                                             "The time vector is not of size 2" );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

        constantType Dt = time[ 1 ];

        //Extract the parameters
        parameterVector macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters;
        parameterVector macroFlowParameters, microFlowParameters, microGradientFlowParameters;
        parameterVector macroYieldParameters, microYieldParameters, microGradientYieldParameters;
        parameterVector Amatrix, Bmatrix, Cmatrix, Dmatrix;
        parameterType alphaMacro, alphaMicro, alphaMicroGradient;
        parameterType relativeTolerance, absoluteTolerance;

        errorOut error = extractMaterialParameters( fparams,
                                                    macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters,
                                                    macroFlowParameters, microFlowParameters, microGradientFlowParameters,
                                                    macroYieldParameters, microYieldParameters, microGradientYieldParameters,
                                                    Amatrix, Bmatrix, Cmatrix, Dmatrix, alphaMacro, alphaMicro, alphaMicroGradient,
                                                    relativeTolerance, absoluteTolerance );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the extraction of the material parameters" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

        //Extract the state variables
        variableType previousMacroStrainISV, previousMicroStrainISV;
        variableVector previousMicroGradientStrainISV;

        variableType previousMacroGamma, previousMicroGamma;
        variableVector previousMicroGradientGamma;

        variableVector previousPlasticDeformationGradient, previousPlasticMicroDeformation, previousPlasticGradientMicroDeformation;

        error = extractStateVariables( SDVS,
                                       previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
                                       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
                                       previousPlasticDeformationGradient, previousPlasticMicroDeformation,
                                       previousPlasticGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the extraction of the state variables" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

        /*===============================================
        | Assemble the fundamental deformation measures |
        ================================================*/
//        std::cout << "assembling fundamental deformation measures\n";

        //Compute the fundamental deformation measures from the degrees of freedom
        variableVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

        error = assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                        previousDeformationGradient, previousMicroDeformation,
                                                        previousGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the computation of the previous deformation measures" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

        variableVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

        error = assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                        currentDeformationGradient, currentMicroDeformation,
                                                        currentGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the computation of the current deformation measures" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

#ifdef DEBUG_MODE

        tempDEBUG.emplace( "currentDeformationGradient", currentDeformationGradient );
        tempDEBUG.emplace( "currentMicroDeformation", currentMicroDeformation );
        tempDEBUG.emplace( "currentGradientMicroDeformation", currentGradientMicroDeformation );

#endif

        /*================================
        | Initialize the model variables |
        =================================*/
//        std::cout << "initializing the model variables\n";

        //Assume that the current strain ISVs are the same as the old
        variableType currentMacroStrainISV = previousMacroStrainISV;
        variableType currentMicroStrainISV = previousMicroStrainISV;
        variableVector currentMicroGradientStrainISV = previousMicroGradientStrainISV;

        //Assume that the new plastic multipliers are zero
        variableType currentMacroGamma = 0.;
        variableType currentMicroGamma = 0.;
        variableVector currentMicroGradientGamma( 3, 0. );

        //Set the current plastic deformation to the previously converged values
        variableVector currentPlasticDeformationGradient      = previousPlasticDeformationGradient;
        variableVector currentPlasticMicroDeformation         = previousPlasticMicroDeformation;
        variableVector currentPlasticGradientMicroDeformation = previousPlasticGradientMicroDeformation;

        //Initialize the previous stresses. If the evolution of plastic strain was previously zero these won't matter
        variableVector previousPK2Stress( dim * dim, 0 );
        variableVector previousReferenceMicroStress( dim * dim, 0 );
        variableVector previousReferenceHigherOrderStress( dim * dim * dim, 0 );

        //Initialize the previous elastic plastic flow directions. If the evolution of plastic strain was previously zero these won't matter.
        variableVector previousMacroFlowDirection( dim * dim, 0 );
        variableVector previousMicroFlowDirection( dim * dim, 0 );
        variableVector previousMicroGradientFlowDirection( dim * dim * dim, 0 );

        //Initialize the previous ISV flow directions
        variableType previousdMacroGdMacroCohesion = 0;
        variableType previousdMicroGdMicroCohesion = 0;
        variableMatrix previousdMicroGradientGdMicroGradientCohesion( previousMicroGradientStrainISV.size(),
                                                                      variableVector( previousMicroGradientStrainISV.size(), 0 ) );

        //Set the current macro cohesion values to zero ( if required this will be updated )
        variableType currentdMacroGdMacroC = 0;
        variableType currentdMicroGdMicroC = 0;
        variableMatrix currentdMicroGradientGdMicroGradientC( dim, variableVector( dim, 0 ) );

        //Initialize the previous plastic velocity gradients
        variableVector previousPlasticMacroVelocityGradient( dim * dim, 0 );
        variableVector previousPlasticMicroVelocityGradient( dim * dim, 0 );
        variableVector previousPlasticMicroGradientVelocityGradient( dim * dim * dim, 0 );

        /*==============================================
        | Begin the evolution of the non-linear values |
        ===============================================*/
//        std::cout << "beginning the evolution of the non-linear values\n";

        //Solve for the previous converged elastic deformation measures
        variableVector previousElasticDeformationGradient, previousElasticMicroDeformation, previousElasticGradientMicroDeformation;

        error = computeElasticPartOfDeformation( previousDeformationGradient, previousMicroDeformation,
                                                 previousGradientMicroDeformation, previousPlasticDeformationGradient,
                                                 previousPlasticMicroDeformation, previousPlasticGradientMicroDeformation,
                                                 previousElasticDeformationGradient, previousElasticMicroDeformation,
                                                 previousElasticGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the computation of the previous elastic deformation meaures" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "previousElasticDeformationGradient", previousElasticDeformationGradient );
        tempDEBUG.emplace( "previousElasticMicroDeformation", previousElasticMicroDeformation );
        tempDEBUG.emplace( "previousElasticGradientMicroDeformation", previousElasticGradientMicroDeformation );

        tempDEBUG.emplace( "previousPlasticDeformationGradient", previousPlasticDeformationGradient );
        tempDEBUG.emplace( "previousPlasticMicroDeformation", previousPlasticMicroDeformation );
        tempDEBUG.emplace( "previousPlasticGradientMicroDeformation", previousPlasticGradientMicroDeformation );
#endif

        //Compute the previous cohesion values
        variableType previousMacroCohesion, previousMicroCohesion;
        variableVector previousMicroGradientCohesion;

        error = computeCohesion( previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
                                 macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters,
                                 previousMacroCohesion, previousMicroCohesion, previousMicroGradientCohesion );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the computation of the previous cohesions" );
            result->addNext( error );
            result->print();               //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
        }

#ifdef DEBUG_MODE
        tardigradeSolverTools::floatVector tmp = { previousMacroCohesion };
        tempDEBUG.emplace( "previousMacroCohesion", tmp );
        tmp = { previousMicroCohesion };
        tempDEBUG.emplace( "previousMicroCohesion", tmp );
        tempDEBUG.emplace( "previousMicroGradientCohesion", previousMicroGradientCohesion );
#endif

        //Assume that the new cohesion is the same as the old
        variableType currentMacroCohesion = previousMacroCohesion;
        variableType currentMicroCohesion = previousMicroCohesion;
        variableVector currentMicroGradientCohesion = previousMicroGradientCohesion;

        //Check if the previous increment had plastic yielding and evolve the plastic deformation if so
//        std::cout << "checking for previous plastic yielding\n";
        if ( ( previousMacroGamma > relativeTolerance * fabs( previousMacroGamma ) + absoluteTolerance ) ||
             ( previousMicroGamma > relativeTolerance * fabs( previousMicroGamma ) + absoluteTolerance ) ||
             ( previousMicroGradientGamma[ 0 ] > relativeTolerance * fabs( previousMicroGradientGamma[ 0 ] ) + absoluteTolerance ) ||
             ( previousMicroGradientGamma[ 1 ] > relativeTolerance * fabs( previousMicroGradientGamma[ 1 ] ) + absoluteTolerance ) ||
             ( previousMicroGradientGamma[ 2 ] > relativeTolerance * fabs( previousMicroGradientGamma[ 2 ] ) + absoluteTolerance ) ){

            //Compute the previous stress
            error = tardigradeMicromorphicLinearElasticity::linearElasticityReference( previousElasticDeformationGradient,
                                                                             previousElasticMicroDeformation,
                                                                             previousElasticGradientMicroDeformation,
                                                                             Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                             previousPK2Stress, previousReferenceMicroStress,
                                                                             previousReferenceHigherOrderStress );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model",
                                                 "Error in the computation of the previous stress meaures" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "previousPK2Stress", previousPK2Stress );
            tempDEBUG.emplace( "previousReferenceMicroStress", previousReferenceMicroStress );
            tempDEBUG.emplace( "previousReferenceHigherOrderStress", previousReferenceHigherOrderStress );
#endif

            //Compute the elastic deformation measures for the plastic evolution
            variableVector previousElasticRightCauchyGreen, previousElasticMicroRightCauchyGreen,
                           previousElasticPsi, previousElasticGamma;

            error = computeElasticDeformationMeasures( previousElasticDeformationGradient, previousElasticMicroDeformation,
                                                       previousElasticGradientMicroDeformation,
                                                       previousElasticRightCauchyGreen, previousElasticMicroRightCauchyGreen,
                                                       previousElasticPsi, previousElasticGamma );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model",
                                                 "Error in the computation of the previous derived elastic deformation meaures" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "previousElasticRightCauchyGreen", previousElasticRightCauchyGreen );
            tempDEBUG.emplace( "previousElasticMicroRightCauchyGreen", previousElasticMicroRightCauchyGreen );
            tempDEBUG.emplace( "previousElasticPsi", previousElasticPsi );
            tempDEBUG.emplace( "previousElasticGamma", previousElasticGamma );
#endif

            //Compute the previous plastic flow directions

            error = computeFlowDirections( previousPK2Stress, previousReferenceMicroStress, previousReferenceHigherOrderStress,
                                           previousMacroCohesion, previousMicroCohesion, previousMicroGradientCohesion,
                                           previousElasticRightCauchyGreen, macroFlowParameters, microFlowParameters,
                                           microGradientFlowParameters, previousMacroFlowDirection, previousMicroFlowDirection,
                                           previousMicroGradientFlowDirection, previousdMacroGdMacroCohesion,
                                           previousdMicroGdMicroCohesion, previousdMicroGradientGdMicroGradientCohesion );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model",
                                                 "Error in the computation of the previous plastic flow directions" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "previousMacroFlowDirection", previousMacroFlowDirection );
            tempDEBUG.emplace( "previousMicroFlowDirection", previousMicroFlowDirection );
            tempDEBUG.emplace( "previousMicroGradientFlowDirection", previousMicroGradientFlowDirection );

            tmp = { previousdMacroGdMacroCohesion };
            tempDEBUG.emplace( "previousdMacroGdMacroCohesion", tmp );

            tmp = { previousdMicroGdMicroCohesion };
            tempDEBUG.emplace( "previousdMicroGdMicroCohesion", tmp );

            tempDEBUG.emplace( "previousdMicroGradientGdMicroGradientCohesion",
                               tardigradeVectorTools::appendVectors( previousdMicroGradientGdMicroGradientCohesion ) );
#endif

            //Update the strain ISV values
            error = evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma, currentMicroGradientGamma,
                                                currentdMacroGdMacroC, currentdMicroGdMicroC, currentdMicroGradientGdMicroGradientC,
                                                previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
                                                previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
                                                previousdMacroGdMacroCohesion, previousdMicroGdMicroCohesion,
                                                previousdMicroGradientGdMicroGradientCohesion,
                                                currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                                alphaMacro, alphaMicro, alphaMicroGradient );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model",
                                                 "Error in the evolution of the current strain state variables" );
                result->addNext( error );
                result->print();               //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tmp = { currentMacroStrainISV };
            tempDEBUG.emplace( "currentMacroStrainISV", tmp );

            tmp = { currentMicroStrainISV };
            tempDEBUG.emplace( "currentMicroStrainISV", currentMicroStrainISV );

            tempDEBUG.emplace( "currentMicroGradientStrainISV", currentMicroGradientStrainISV );
#endif

            //Compute the current cohesion values
            error = computeCohesion( currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                     macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters,
                                     currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model",
                                                 "Error in the evolution of the updated cohesion" );
                result->addNext( error );
                result->print();               //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tmp = { currentMacroCohesion };
            tempDEBUG.emplace( "currentMacroCohesion", tmp );

            tmp = { currentMicroCohesion };
            tempDEBUG.emplace( "currentMicroCohesion", currentMicroCohesion );

            tempDEBUG.emplace( "currentMicroGradientCohesion", currentMicroGradientCohesion );
#endif

            //Compute the previous plastic velocity gradients
            error = computePlasticVelocityGradients( previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
                                                     previousElasticRightCauchyGreen, previousElasticMicroRightCauchyGreen,
                                                     previousElasticPsi, previousElasticGamma, previousMacroFlowDirection,
                                                     previousMicroFlowDirection, previousMicroGradientFlowDirection,
                                                     previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                                     previousPlasticMicroGradientVelocityGradient );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model",
                                                 "Error in the computation of the previous plastic velocity gradient" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }
#ifdef DEBUG_MODE
            tempDEBUG.emplace( "previousMacroFlowDirection", previousMacroFlowDirection );
            tempDEBUG.emplace( "previousMicroFlowDirection", previousMicroFlowDirection );
            tempDEBUG.emplace( "previousMicroGradientFlowDirection", previousMicroGradientFlowDirection );
#endif

            //Update the current plastic deformation measures
            variableVector currentPlasticMacroVelocityGradient( dim * dim, 0 );
            variableVector currentPlasticMicroVelocityGradient( dim * dim, 0 );
            variableVector currentPlasticMicroGradientVelocityGradient( dim * dim * dim, 0 );

            error = evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient, currentPlasticMicroVelocityGradient,
                                              currentPlasticMicroGradientVelocityGradient, previousPlasticDeformationGradient,
                                              previousPlasticMicroDeformation, previousPlasticGradientMicroDeformation,
                                              previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                              previousPlasticMicroGradientVelocityGradient,
                                              currentPlasticDeformationGradient, currentPlasticMicroDeformation,
                                              currentPlasticGradientMicroDeformation,
                                              alphaMacro, alphaMicro, alphaMicroGradient );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model",
                                                 "Error in the computation of the current plastic deformation" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }
        }

        //Update the elastic deformation
//        std::cout << "updating the elastic deformation\n";
        variableVector currentElasticDeformationGradient, currentElasticMicroDeformation, currentElasticGradientMicroDeformation;

        error = computeElasticPartOfDeformation( currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation,
                                                 currentPlasticDeformationGradient, currentPlasticMicroDeformation,
                                                 currentPlasticGradientMicroDeformation, currentElasticDeformationGradient,
                                                 currentElasticMicroDeformation, currentElasticGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the computation of the current elastic deformation" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "currentElasticDeformationGradient", currentElasticDeformationGradient );
        tempDEBUG.emplace( "currentElasticMicroDeformation", currentElasticMicroDeformation );
        tempDEBUG.emplace( "currentElasticGradientMicroDeformation", currentElasticGradientMicroDeformation );

        tempDEBUG.emplace( "currentPlasticDeformationGradient", currentPlasticDeformationGradient );
        tempDEBUG.emplace( "currentPlasticMicroDeformation", currentPlasticMicroDeformation );
        tempDEBUG.emplace( "currentPlasticGradientMicroDeformation", currentPlasticGradientMicroDeformation );
#endif

//        std::cout << "\n\nDEFORMATION MEASURES\n";
//        std::cout << "  "; tardigradeVectorTools::print( currentDeformationGradient );
//        std::cout << "  "; tardigradeVectorTools::print( currentMicroDeformation );
//        std::cout << "  "; tardigradeVectorTools::print( currentGradientMicroDeformation );
//
//        std::cout << "TRIAL ELASTIC DEFORMATION MEASURES\n";
//        std::cout << "  "; tardigradeVectorTools::print( currentElasticDeformationGradient );
//        std::cout << "  "; tardigradeVectorTools::print( currentElasticMicroDeformation );
//        std::cout << "  "; tardigradeVectorTools::print( currentElasticGradientMicroDeformation );
//
//        std::cout << "TRIAL PLASTIC DEFORMATION MEASURES\n";
//        std::cout << "  "; tardigradeVectorTools::print( currentPlasticDeformationGradient );
//        std::cout << "  "; tardigradeVectorTools::print( currentPlasticMicroDeformation );
//        std::cout << "  "; tardigradeVectorTools::print( currentPlasticGradientMicroDeformation );
//
//        std::cout << "PREVIOUS PLASTIC DEFORMATION MEASURES\n";
//        std::cout << "  "; tardigradeVectorTools::print( previousPlasticDeformationGradient );
//        std::cout << "  "; tardigradeVectorTools::print( previousPlasticMicroDeformation );
//        std::cout << "  "; tardigradeVectorTools::print( previousPlasticGradientMicroDeformation );

        //Compute the right Cauchy-Green deformation gradient
        variableVector currentElasticRightCauchyGreen;

        error = tardigradeConstitutiveTools::computeRightCauchyGreen( currentElasticDeformationGradient, currentElasticRightCauchyGreen );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the computation of the current elastic right Cauchy-Green deformation" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        //Compute the new stress values
        variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

        error = tardigradeMicromorphicLinearElasticity::linearElasticityReference( currentElasticDeformationGradient,
                                                                         currentElasticMicroDeformation,
                                                                         currentElasticGradientMicroDeformation,
                                                                         Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                         currentPK2Stress, currentReferenceMicroStress,
                                                                         currentReferenceHigherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the computation of the current stress measures" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

//        std::cout << "TRIAL STRESSES\n";
//        std::cout << "  "; tardigradeVectorTools::print( currentPK2Stress );
//        std::cout << "  "; tardigradeVectorTools::print( currentReferenceMicroStress );
//        std::cout << "  "; tardigradeVectorTools::print( currentReferenceHigherOrderStress );

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "intermediatePK2Stress", currentPK2Stress );
        tempDEBUG.emplace( "intermediateReferenceMicroStress", currentReferenceMicroStress );
        tempDEBUG.emplace( "intermediateReferenceHigherOrderStress", currentReferenceHigherOrderStress );
#endif

//        std::cout << "current stress measures\n";
//        std::cout << "currentPK2Stress:\n"; tardigradeVectorTools::print( currentPK2Stress );
//        std::cout << "currentReferenceMicroStress:\n"; tardigradeVectorTools::print( currentReferenceMicroStress );
//        std::cout << "currentReferenceHigherOrderStress:\n"; tardigradeVectorTools::print( currentReferenceHigherOrderStress );

        //Evaluate the yield functions
//        std::cout << "evaluating the yield functions\n";
        variableVector currentYieldFunctionValues;

        error = evaluateYieldFunctions( currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress,
                                        currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion,
                                        currentElasticRightCauchyGreen,
                                        macroYieldParameters, microYieldParameters, microGradientYieldParameters,
                                        currentYieldFunctionValues
#ifdef DEBUG_MODE
                                        , tempDEBUG
#endif
                                       );

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "currentYieldFunctionValues", currentYieldFunctionValues );
        tempITERATION.emplace( "pre_iteration_values", tempDEBUG );
        DEBUG.emplace( "pre_iteration_values", tempITERATION );
        tempDEBUG.clear();
        tempITERATION.clear();
#endif

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the computation of the current yield function values" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

//        std::cout << "TRIAL YIELD VALUES\n";
//        std::cout << "  "; tardigradeVectorTools::print( currentYieldFunctionValues );

//        std::cout << "initial yield function values:\n"; tardigradeVectorTools::print( currentYieldFunctionValues );

        /*============================
        | Begin the non-linear solve |
        ============================*/
//        std::cout << "beginning the nonlinear solve\n";

        //Check if any of the surfaces are yielding and begin the non-linear solver if they are
        tardigradeSolverTools::floatVector solutionVector;
        tardigradeSolverTools::intVector activePlasticity( currentYieldFunctionValues.size(), 0 );
        bool convergenceFlag = true;

        for ( unsigned int i = 0; i < currentYieldFunctionValues.size(); i++ ){
            if ( currentYieldFunctionValues[ i ] > absoluteTolerance ){

                convergenceFlag = false;
                
                activePlasticity[ i ] = 1;
            }
        }
//        std::cout << "convergenceFlag: " << convergenceFlag << "\n";
//        std::cout << "activePlasticity: "; tardigradeVectorTools::print( activePlasticity );

        if ( !convergenceFlag ){
            tardigradeSolverTools::floatMatrix floatArgs =
                {
                    { Dt },
                    currentDeformationGradient,
                    currentMicroDeformation,
                    currentGradientMicroDeformation,
                    { previousMacroGamma },
                    { previousMicroGamma },
                    previousMicroGradientGamma,
                    previousPlasticDeformationGradient,
                    previousPlasticMicroDeformation,
                    previousPlasticGradientMicroDeformation,
                    { previousMacroStrainISV },
                    { previousMicroStrainISV },
                    previousMicroGradientStrainISV,
                    { previousdMacroGdMacroCohesion },
                    { previousdMicroGdMicroCohesion },
                    tardigradeVectorTools::appendVectors( previousdMicroGradientGdMicroGradientCohesion ),
                    previousPlasticMacroVelocityGradient,
                    previousPlasticMicroVelocityGradient,
                    previousPlasticMicroGradientVelocityGradient,
                    macroFlowParameters,
                    microFlowParameters,
                    microGradientFlowParameters,
                    macroHardeningParameters,
                    microHardeningParameters,
                    microGradientHardeningParameters,
                    macroYieldParameters,
                    microYieldParameters,
                    microGradientYieldParameters,
                    Amatrix,
                    Bmatrix,
                    Cmatrix,
                    Dmatrix,
                    { alphaMacro },
                    { alphaMicro },
                    { alphaMicroGradient }
                };
        
            tardigradeSolverTools::floatMatrix floatOuts =
                {
                    currentPK2Stress,
                    currentReferenceMicroStress,
                    currentReferenceHigherOrderStress,
                    { currentMacroStrainISV },
                    { currentMicroStrainISV },
                    currentMicroGradientStrainISV,
                    {}, //The plastic deformation gradient
                    {}, //The plastic micro deformation
                    {}  //The gradient of plastic micro deformation
                };

            tardigradeSolverTools::intMatrix intOuts = { { } };

            tardigradeSolverTools::stdFncNLFJ func
                = static_cast< tardigradeSolverTools::NonLinearFunctionWithJacobian >( computePlasticMultiplierResidual );

            tardigradeSolverTools::floatVector x0 =
                {
                    currentMacroGamma,
                    currentMicroGamma,
                    currentMicroGradientGamma[ 0 ],
                    currentMicroGradientGamma[ 1 ],
                    currentMicroGradientGamma[ 2 ]
                };

            tardigradeSolverTools::floatVector solutionVector;

            bool convergeFlag, fatalErrorFlag;
           
            tardigradeSolverTools::intMatrix intArgs = { { 0 } };

            tardigradeSolverTools::solverType linearSolver;
            tardigradeSolverTools::floatMatrix J;
            tardigradeSolverTools::intVector boundVariableIndices = {  0,  1,  2,  3,  4 };
            tardigradeSolverTools::intVector boundSigns           = {  0,  0,  0,  0,  0 };
            tardigradeSolverTools::floatVector boundValues        = {  0,  0,  0,  0,  0 };

            error = tardigradeSolverTools::homotopySolver( func, x0, solutionVector, convergeFlag, fatalErrorFlag,
                                                 floatOuts, intOuts, floatArgs, intArgs, linearSolver, J,
                                                 boundVariableIndices, boundSigns, boundValues, true,
#ifdef DEBUG_MODE
                                                 DEBUG,
#endif
                                                 20, relativeTolerance, absoluteTolerance,
                                                 1e-4, 5, 1.0, .01 );

            if ( ( error ) && fatalErrorFlag ){ //Fatal error
                errorOut result = new errorNode( "evaluate_model",
                                                 "Error in nonlinear solve" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }
            else if ( ( error ) && !convergeFlag ){ //Solution didn't converge
                errorOut result = new errorNode( "evaluate_model",
                                                 "Solution failed to converge. Requesting timestep cutback" );
                result->addNext( error );
                result->print();          //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing

                //Additional messages
                output_message += "\ncurrentDeformationGradient:\n";
                for ( unsigned int kk = 0; kk < currentDeformationGradient.size(); kk++ ){
                    output_message += std::to_string( currentDeformationGradient[ kk ] ) + ", ";
                }
                output_message += "\n";
                output_message += "\ncurrentMicroDeformation:\n";
                for ( unsigned int kk = 0; kk < currentMicroDeformation.size(); kk++ ){
                    output_message += std::to_string( currentMicroDeformation[ kk ] ) + ", ";
                }
                output_message += "\n";
                output_message += "\ncurrentMicroGradientDeformation:\n";
                for ( unsigned int kk = 0; kk < currentGradientMicroDeformation.size(); kk++ ){
                    output_message += std::to_string( currentGradientMicroDeformation[ kk ] ) + ", ";
                }
                output_message += "\n";
                output_message += "\nSDVS:\n";
                for ( unsigned int kk = 0; kk < SDVS.size(); kk++ ){
                    output_message += std::to_string( SDVS[ kk ] ) + ", ";
                }
                output_message += "\n";
                return 1;
            }

            //Extract the deformation measures
            currentPlasticDeformationGradient      = floatOuts[ 6 ];
            currentPlasticMicroDeformation         = floatOuts[ 7 ];
            currentPlasticGradientMicroDeformation = floatOuts[ 8 ];
            currentMacroGamma                      = solutionVector[ 0 ];
            currentMicroGamma                      = solutionVector[ 1 ];
            currentMicroGradientGamma              = variableVector( solutionVector.begin() + 2, solutionVector.begin() + 5 );

            if ( currentMacroGamma < -absoluteTolerance ){
                errorOut result = new errorNode( "evaluate_model", "The macro-scale plastic multiplier is negative" );
                result->print();
                output_message = buffer.str();
                return 1;
            }

            if ( currentMicroGamma < -absoluteTolerance ){
                errorOut result = new errorNode( "evaluate_model", "The micro-scale plastic multiplier is negative" );
                result->print();
                output_message = buffer.str();
                return 1;
            }

            for ( unsigned int i = 0; i < currentMicroGradientGamma.size(); i++ ){
                if ( currentMicroGradientGamma[ i ] < -absoluteTolerance ){
                    errorOut result = new errorNode( "evaluate_model", "A micro-gradient plastic multiplier is negative" );
                    result->print();
                    output_message = buffer.str();
                    return 1;
                }
            }

            //Extract the stresses and state variables
            currentPK2Stress                       = floatOuts[ 0 ];
            currentReferenceMicroStress            = floatOuts[ 1 ];
            currentReferenceHigherOrderStress      = floatOuts[ 2 ];
            currentMacroStrainISV                  = floatOuts[ 3 ][ 0 ];
            currentMicroStrainISV                  = floatOuts[ 4 ][ 0 ];
            currentMicroGradientStrainISV          = floatOuts[ 5 ];

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "convergedPlasticDeformationGradient", currentPlasticDeformationGradient );
            tempDEBUG.emplace( "convergedPlasticMicroDeformation", currentPlasticMicroDeformation );
            tempDEBUG.emplace( "convergedPlasticGradientMicroDeformation", currentPlasticGradientMicroDeformation );

            tardigradeSolverTools::floatVector temp = { currentMacroStrainISV };
            tempDEBUG.emplace( "convergedMacroStrainISV", temp );

            temp = { currentMicroStrainISV };
            tempDEBUG.emplace( "convergedMicroStrainISV", temp );

            tempDEBUG.emplace( "convergedMicroGradientStrainISV", currentMicroGradientStrainISV );

            temp = { currentMacroGamma };
            tempDEBUG.emplace( "convergedMacroGamma", temp );

            temp = { currentMicroGamma };
            tempDEBUG.emplace( "convergedMicroGamma", temp );

            tempDEBUG.emplace( "convergedMicroGradientGamma", currentMicroGradientGamma );

            tempDEBUG.emplace( "intermediatePK2Stress", currentPK2Stress );
            tempDEBUG.emplace( "intermediateReferenceMicroStress", currentReferenceMicroStress );
            tempDEBUG.emplace( "intermediateReferenceHigherOrderStress", currentReferenceHigherOrderStress );

            tempITERATION.emplace( "converged_values", tempDEBUG );
            DEBUG.emplace( "converged_values", tempITERATION );

            tempITERATION.clear();
            tempDEBUG.clear();
#endif
        }

        /*============================
        | Assemble the output values |
        ============================*/
        
        //Assemble the SDV vector
        std::vector< double > currentStrainISVS = { currentMacroStrainISV, currentMicroStrainISV };
        currentStrainISVS = tardigradeVectorTools::appendVectors( { currentStrainISVS, currentMicroGradientStrainISV } );

        std::vector< double > currentGammas = { currentMacroGamma, currentMicroGamma };
        currentGammas = tardigradeVectorTools::appendVectors( { currentGammas, currentMicroGradientGamma } );

        SDVS = tardigradeVectorTools::appendVectors( { currentStrainISVS, currentGammas, currentPlasticDeformationGradient - eye,
                                             currentPlasticMicroDeformation - eye, currentPlasticGradientMicroDeformation } );

        //Output the stresses mapped to the true reference configuration
        error = tardigradeMicromorphicTools::pullBackCauchyStress( currentPK2Stress, currentPlasticDeformationGradient, current_PK2 );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in pullback operation on the PK2 stress" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        error = tardigradeMicromorphicTools::pullBackMicroStress( currentReferenceMicroStress, currentPlasticDeformationGradient, current_SIGMA );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in pullback operation on the reference symmetric micro-stress" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        error = tardigradeMicromorphicTools::pullBackHigherOrderStress( currentReferenceHigherOrderStress, 
                                                              currentPlasticDeformationGradient, 
                                                              currentPlasticMicroDeformation, current_M );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in pullback operation on the reference higher order stress" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

//        std::cout << "SDVS:\n"; tardigradeVectorTools::print( SDVS );
        //Model evaluation successful. Return.
        return 0;
    }

    int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                        const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                        const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                        const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                        std::vector< double > &SDVS,
                        const std::vector< double > &current_ADD_DOF,
                        const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                        const std::vector< double > &previous_ADD_DOF,
                        const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                        std::vector< double > &current_PK2, std::vector< double > &current_SIGMA, std::vector< double > &current_M,
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
                        ){
        /*!
         * Evaluate the elasto-plastic constitutive model. Note the format of the header changed to provide a 
         * consistant interface with the material model library.
         *
         * :param const std::vector< double > &time: The current time and the timestep
         *     [ current_t, dt ]
         * :param const std::vector< double > ( &fparams ): The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters,
         *       num_macroHardeningParameters, macroHardeningParameters,
         *       num_microHardeningParameters, microHardeningParameters,
         *       num_microGradientHardeningParameters, microGradientHardeningParameters,
         *       num_macroFlowParameters, macroFlowParameters,
         *       num_microFlowParameters, microFlowParameters,
         *       num_microGradientFlowParameters, microGradientFlowParameters,
         *       num_macroYieldParameters, macroYieldParameters,
         *       num_microYieldParameters, microYieldParameters,
         *       num_microGradientYieldParameters, microGradientYieldParameters,
         *       alphaMacro, alphaMicro, alphaMicroGradient,
         *       relativeTolerance, absoluteTolerance ]
         *
         * :param const double ( &current_grad_u )[ 3 ][ 3 ]: The current displacement gradient
         *     Assumed to be of the form [ [ u_{1,1}, u_{1,2}, u_{1,3} ],
         *                                 [ u_{2,1}, u_{2,2}, u_{2,3} ],
         *                                 [ u_{3,1}, u_{3,2}, u_{3,3} ] ]
         * :param const double ( &current_phi )[ 9 ]: The current micro displacment values.
         *     Assumed to be of the form [ \phi_{11}, \phi_{12}, \phi_{13}, \phi_{21}, \phi_{22}, \phi_{23}, \phi_{31}, \phi_{32}, \phi_{33} ]
         * :param const double ( &current_grad_phi )[ 9 ][ 3 ]: The current micro displacement gradient
         *     Assumed to be of the form [ [ \phi_{11,1}, \phi_{11,2}, \phi_{11,3} ],
         *                                 [ \phi_{12,1}, \phi_{12,2}, \phi_{12,3} ],
         *                                 [ \phi_{13,1}, \phi_{13,2}, \phi_{13,3} ],
         *                                 [ \phi_{21,1}, \phi_{21,2}, \phi_{21,3} ],
         *                                 [ \phi_{22,1}, \phi_{22,2}, \phi_{22,3} ],
         *                                 [ \phi_{23,1}, \phi_{23,2}, \phi_{23,3} ],
         *                                 [ \phi_{31,1}, \phi_{31,2}, \phi_{31,3} ],
         *                                 [ \phi_{32,1}, \phi_{32,2}, \phi_{32,3} ],
         *                                 [ \phi_{33,1}, \phi_{33,2}, \phi_{33,3} ] ]
         * :param const double ( &previous_grad_u )[ 3 ][ 3 ]: The previous displacement gradient.
         * :param const double ( &previous_phi )[ 9 ]: The previous micro displacement.
         * :param const double ( &previous_grad_phi )[ 9 ][ 3 ]: The previous micro displacement gradient.
         * :param std::vector< double > &SDVS: The previously converged values of the state variables
         *     [ previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
         *       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
         *       previousPlasticDeformationGradient - eye, previousPlasticMicroDeformation - eye,
         *       previousPlasticMicroGradient ]
         * :param std::vector< double > &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * :param std::vector< std::vector< double > > &current_ADD_grad_DOF: The current values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * :param std::vector< double > &current_PK2: The current value of the second Piola Kirchhoff stress tensor. The format is
         *     [ S_{11}, S_{12}, S_{13}, S_{21}, S_{22}, S_{23}, S_{31}, S_{32}, S_{33} ]
         * :param std::vector< double > &current_SIGMA: The current value of the reference micro stress. The format is
         *     [ S_{11}, S_{12}, S_{13}, S_{21}, S_{22}, S_{23}, S_{31}, S_{32}, S_{33} ]
         * :param std::vector< double > &current_M: The current value of the reference higher order stress. The format is
         *     [ M_{111}, M_{112}, M_{113}, M_{121}, M_{122}, M_{123}, M_{131}, M_{132}, M_{133},
         *       M_{211}, M_{212}, M_{213}, M_{221}, M_{222}, M_{223}, M_{231}, M_{232}, M_{233},
         *       M_{311}, M_{312}, M_{313}, M_{321}, M_{322}, M_{323}, M_{331}, M_{332}, M_{333} ]
         * :param std::vector< std::vector< double > > &DPK2Dgrad_u: The Jacobian of the PK2 stress w.r.t. the 
         *     gradient of macro displacement.
         * :param std::vector< std::vector< double > > &DPK2Dphi: The Jacobian of the PK2 stress w.r.t. the
         *     micro displacement.
         * :param std::vector< std::vector< double > > &DPK2Dgrad_phi: The Jacobian of the PK2 stress w.r.t.
         *     the gradient of the micro displacement.
         * :param std::vector< std::vector< double > > &DSIGMAdgrad_u: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the macro displacement.
         * :param std::vector< std::vector< double > > &DSIGMAdphi: The Jacobian of the reference symmetric micro
         *     stress w.r.t. the micro displacement.
         * :param std::vector< std::vector< double > > &DSIGMAdgrad_phi: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the micro displacement.
         * :param std::vector< std::vector< double > > &DMDgrad_u: The Jacobian of the reference higher order
         *     stress w.r.t. the gradient of the macro displacement.
         * :param std::vector< std::vector< double > > &DMDphi: The Jacobian of the reference higher order stress
         *     w.r.t. the micro displacement.
         * :param std::vector< std::vector< double > > &DMDgrad_phi: The Jacobian of the reference higher order stress
         *     w.r.t. the gradient of the micro displacement.
         * :param std::vector< std::vector< double > > &ADD_TERMS: Additional terms ( unused )
         * :param std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS: The jacobians of the additional
         *     terms w.r.t. the deformation. This is currently being used to support the gradient enhanced damage work
         *     by returning the Jacobians of the plastic deformation gradients w.r.t. the deformation measures. The
         *     ordering is: DFpDgrad_u, DFpDphi, DFpDgrad_phi, DchipDgrad_u, DchipDphi, DchipDgrad_phi, Dgrad_chipDgrad_u, Dgrad_chipDchi, Dgrad_chipDgrad_chi
         * :param std::string &output_message: The output message string.
         * :param tardigradeSolverTools::homotopyMap DEBUG: The debugging map ( only available if DEBUG_MODE is defined )
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        //Assume 3D
        unsigned int dim = 3;

#ifdef DEBUG_MODE
        tardigradeSolverTools::debugMap tempDEBUG;
        tardigradeSolverTools::iterationMap tempITERATION;
#endif

        //Construct identity matrix
        constantVector eye( dim * dim, 0 );
        tardigradeVectorTools::eye< constantType >( eye );

        //Re-direct the output to a buffer
        std::stringbuf buffer;
        cerr_redirect rd( &buffer );

        /*=============================
        | Extract the incoming values |
        ==============================*/

        //Extract the time
        if ( time.size() != 2 ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "The time vector is not of size 2" );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

        constantType Dt = time[ 1 ];

        //Extract the parameters
        parameterVector macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters;
        parameterVector macroFlowParameters, microFlowParameters, microGradientFlowParameters;
        parameterVector macroYieldParameters, microYieldParameters, microGradientYieldParameters;
        parameterVector Amatrix, Bmatrix, Cmatrix, Dmatrix;
        parameterType alphaMacro, alphaMicro, alphaMicroGradient;
        parameterType relativeTolerance, absoluteTolerance;

        errorOut error = extractMaterialParameters( fparams,
                                                    macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters,
                                                    macroFlowParameters, microFlowParameters, microGradientFlowParameters,
                                                    macroYieldParameters, microYieldParameters, microGradientYieldParameters,
                                                    Amatrix, Bmatrix, Cmatrix, Dmatrix, alphaMacro, alphaMicro, alphaMicroGradient,
                                                    relativeTolerance, absoluteTolerance );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "Error in the extraction of the material parameters" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

        //Extract the state variables
        variableType previousMacroStrainISV, previousMicroStrainISV;
        variableVector previousMicroGradientStrainISV;

        variableType previousMacroGamma, previousMicroGamma;
        variableVector previousMicroGradientGamma;

        variableVector previousPlasticDeformationGradient, previousPlasticMicroDeformation, previousPlasticGradientMicroDeformation;

        error = extractStateVariables( SDVS,
                                       previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
                                       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
                                       previousPlasticDeformationGradient, previousPlasticMicroDeformation,
                                       previousPlasticGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "Error in the extraction of the state variables" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

        /*===============================================
        | Assemble the fundamental deformation measures |
        ================================================*/

        //Compute the fundamental deformation measures from the degrees of freedom
        variableVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

        error = assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                        previousDeformationGradient, previousMicroDeformation,
                                                        previousGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "Error in the computation of the previous deformation measures" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

        variableVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;
        variableMatrix dDeformationGradientdGradientMacroDisplacement, dMicroDeformationdMicroDisplacement,
                       dGradientMicroDeformationdGradientMicroDisplacement;

        error = assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                        currentDeformationGradient, currentMicroDeformation,
                                                        currentGradientMicroDeformation,
                                                        dDeformationGradientdGradientMacroDisplacement,
                                                        dMicroDeformationdMicroDisplacement,
                                                        dGradientMicroDeformationdGradientMicroDisplacement );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "Error in the computation of the current deformation measures" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "currentDeformationGradient_", currentDeformationGradient );
        tempDEBUG.emplace( "currentMicroDeformation_", currentMicroDeformation );
        tempDEBUG.emplace( "currentGradientMicroDeformation_", currentGradientMicroDeformation );

        tempDEBUG.emplace( "totaldDeformationGradientdGradientMacroDisplacement",
                           tardigradeVectorTools::appendVectors( dDeformationGradientdGradientMacroDisplacement ) );
        tempDEBUG.emplace( "totaldMicroDeformationdMicroDisplacement",
                           tardigradeVectorTools::appendVectors( dMicroDeformationdMicroDisplacement ) );
        tempDEBUG.emplace( "totaldGradientMicroDeformationdGradientMicroDisplacement",
                           tardigradeVectorTools::appendVectors( dGradientMicroDeformationdGradientMicroDisplacement ) );
#endif

        /*================================
        | Initialize the model variables |
        =================================*/

        //Assume that the current strain ISVs are the same as the old
        variableType currentMacroStrainISV = previousMacroStrainISV;
        variableType currentMicroStrainISV = previousMicroStrainISV;
        variableVector currentMicroGradientStrainISV = previousMicroGradientStrainISV;

        //Assume that the new plastic multipliers are zero
        variableType currentMacroGamma = 0.;
        variableType currentMicroGamma = 0.;
        variableVector currentMicroGradientGamma( 3, 0. );

        //Set the current plastic deformation to the previously converged values
        variableVector currentPlasticDeformationGradient      = previousPlasticDeformationGradient;
        variableVector currentPlasticMicroDeformation         = previousPlasticMicroDeformation;
        variableVector currentPlasticGradientMicroDeformation = previousPlasticGradientMicroDeformation;

        //Initialize the jacobians of the plastic deformation measures
        variableMatrix dPlasticDeformationGradientdDeformationGradient( dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dPlasticDeformationGradientdMicroDeformation( dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dPlasticDeformationGradientdGradientMicroDeformation( dim * dim, variableVector( dim * dim * dim, 0 ) );

        variableMatrix dPlasticMicroDeformationdDeformationGradient( dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dPlasticMicroDeformationdMicroDeformation( dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dPlasticMicroDeformationdGradientMicroDeformation( dim * dim, variableVector( dim * dim * dim, 0 ) );

        variableMatrix dPlasticGradientMicroDeformationdDeformationGradient( dim * dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dPlasticGradientMicroDeformationdMicroDeformation( dim * dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dPlasticGradientMicroDeformationdGradientMicroDeformation( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        //Initialize the previous stresses. If the evolution of plastic strain was previously zero these won't matter
        variableVector previousPK2Stress( dim * dim, 0 );
        variableVector previousReferenceMicroStress( dim * dim, 0 );
        variableVector previousReferenceHigherOrderStress( dim * dim * dim, 0 );

        //Initialize the previous elastic plastic flow directions. If the evolution of plastic strain was previously zero these won't matter.
        variableVector previousMacroFlowDirection( dim * dim, 0 );
        variableVector previousMicroFlowDirection( dim * dim, 0 );
        variableVector previousMicroGradientFlowDirection( dim * dim * dim, 0 );

        //Initialize the previous ISV flow directions
        variableType previousdMacroGdMacroCohesion = 0;
        variableType previousdMicroGdMicroCohesion = 0;
        variableMatrix previousdMicroGradientGdMicroGradientCohesion( previousMicroGradientStrainISV.size(),
                                                                      variableVector( previousMicroGradientStrainISV.size(), 0 ) );

        //Set the current macro cohesion values to zero ( if required this will be updated )
        variableType currentdMacroGdMacroC = 0;
        variableType currentdMicroGdMicroC = 0;
        variableMatrix currentdMicroGradientGdMicroGradientC( dim, variableVector( dim, 0 ) );

        //Initialize the previous plastic velocity gradients
        variableVector previousPlasticMacroVelocityGradient( dim * dim, 0 );
        variableVector previousPlasticMicroVelocityGradient( dim * dim, 0 );
        variableVector previousPlasticMicroGradientVelocityGradient( dim * dim * dim, 0 );

        //Initialize the intermediate Jacobians
        variableMatrix dPK2StressdDeformationGradient
            = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dPK2StressdMicroDeformation
            = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dPK2StressdGradientMicroDeformation
            = variableMatrix( dim * dim, variableVector( dim * dim * dim, 0 ) );

        variableMatrix dReferenceMicroStressdDeformationGradient
            = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dReferenceMicroStressdMicroDeformation
            = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dReferenceMicroStressdGradientMicroDeformation
            = variableMatrix( dim * dim, variableVector( dim * dim * dim, 0 ) ); 

        variableMatrix dReferenceHigherOrderStressdDeformationGradient
            = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dReferenceHigherOrderStressdMicroDeformation
            = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dReferenceHigherOrderStressdGradientMicroDeformation
            = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        /*==============================================
        | Begin the evolution of the non-linear values |
        ===============================================*/

        //Solve for the previous converged elastic deformation measures
        variableVector previousElasticDeformationGradient, previousElasticMicroDeformation, previousElasticGradientMicroDeformation;

        error = computeElasticPartOfDeformation( previousDeformationGradient, previousMicroDeformation,
                                                 previousGradientMicroDeformation, previousPlasticDeformationGradient,
                                                 previousPlasticMicroDeformation, previousPlasticGradientMicroDeformation,
                                                 previousElasticDeformationGradient, previousElasticMicroDeformation,
                                                 previousElasticGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "Error in the computation of the previous elastic deformation meaures" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
            return 2;
        }

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "previousElasticDeformationGradient", previousElasticDeformationGradient );
        tempDEBUG.emplace( "previousElasticMicroDeformation", previousElasticMicroDeformation );
        tempDEBUG.emplace( "previousElasticGradientMicroDeformation", previousElasticGradientMicroDeformation );

        tempDEBUG.emplace( "previousPlasticDeformationGradient", previousPlasticDeformationGradient );
        tempDEBUG.emplace( "previousPlasticMicroDeformation", previousPlasticMicroDeformation );
        tempDEBUG.emplace( "previousPlasticGradientMicroDeformation", previousPlasticGradientMicroDeformation );
#endif

        //Compute the previous cohesion values
        variableType previousMacroCohesion, previousMicroCohesion;
        variableVector previousMicroGradientCohesion;

        error = computeCohesion( previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
                                 macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters,
                                 previousMacroCohesion, previousMicroCohesion, previousMicroGradientCohesion );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "Error in the computation of the previous cohesions" );
            result->addNext( error );
            result->print();               //Print the error message
            output_message = buffer.str(); //Save the output to enable message printing
        }

#ifdef DEBUG_MODE
        tardigradeSolverTools::floatVector tmp = { previousMacroCohesion };
        tempDEBUG.emplace( "previousMacroCohesion", tmp );
        tmp = { previousMicroCohesion };
        tempDEBUG.emplace( "previousMicroCohesion", tmp );
        tempDEBUG.emplace( "previousMicroGradientCohesion", previousMicroGradientCohesion );
#endif

        //Assume that the new cohesion is the same as the old
        variableType currentMacroCohesion = previousMacroCohesion;
        variableType currentMicroCohesion = previousMicroCohesion;
        variableVector currentMicroGradientCohesion = previousMicroGradientCohesion;

        //Check if the previous increment had plastic yielding and evolve the plastic deformation if so
        if ( ( previousMacroGamma > relativeTolerance * fabs( previousMacroGamma ) + absoluteTolerance ) ||
             ( previousMicroGamma > relativeTolerance * fabs( previousMicroGamma ) + absoluteTolerance ) ||
             ( previousMicroGradientGamma[ 0 ] > relativeTolerance * fabs( previousMicroGradientGamma[ 0 ] ) + absoluteTolerance ) ||
             ( previousMicroGradientGamma[ 1 ] > relativeTolerance * fabs( previousMicroGradientGamma[ 1 ] ) + absoluteTolerance ) ||
             ( previousMicroGradientGamma[ 2 ] > relativeTolerance * fabs( previousMicroGradientGamma[ 2 ] ) + absoluteTolerance ) ){

            //Compute the previous stress
            error = tardigradeMicromorphicLinearElasticity::linearElasticityReference( previousElasticDeformationGradient,
                                                                             previousElasticMicroDeformation,
                                                                             previousElasticGradientMicroDeformation,
                                                                             Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                             previousPK2Stress, previousReferenceMicroStress,
                                                                             previousReferenceHigherOrderStress );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model (jacobian)",
                                                 "Error in the computation of the previous stress meaures" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "previousPK2Stress", previousPK2Stress );
            tempDEBUG.emplace( "previousReferenceMicroStress", previousReferenceMicroStress );
            tempDEBUG.emplace( "previousReferenceHigherOrderStress", previousReferenceHigherOrderStress );
#endif

            //Compute the elastic deformation measures for the plastic evolution
            variableVector previousElasticRightCauchyGreen, previousElasticMicroRightCauchyGreen,
                           previousElasticPsi, previousElasticGamma;

            error = computeElasticDeformationMeasures( previousElasticDeformationGradient, previousElasticMicroDeformation,
                                                       previousElasticGradientMicroDeformation,
                                                       previousElasticRightCauchyGreen, previousElasticMicroRightCauchyGreen,
                                                       previousElasticPsi, previousElasticGamma );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model (jacobian)",
                                                 "Error in the computation of the previous derived elastic deformation meaures" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "previousElasticRightCauchyGreen", previousElasticRightCauchyGreen );
            tempDEBUG.emplace( "previousElasticMicroRightCauchyGreen", previousElasticMicroRightCauchyGreen );
            tempDEBUG.emplace( "previousElasticPsi", previousElasticPsi );
            tempDEBUG.emplace( "previousElasticGamma", previousElasticGamma );
#endif

            //Compute the previous plastic flow directions

            error = computeFlowDirections( previousPK2Stress, previousReferenceMicroStress, previousReferenceHigherOrderStress,
                                           previousMacroCohesion, previousMicroCohesion, previousMicroGradientCohesion,
                                           previousElasticRightCauchyGreen, macroFlowParameters, microFlowParameters,
                                           microGradientFlowParameters, previousMacroFlowDirection, previousMicroFlowDirection,
                                           previousMicroGradientFlowDirection, previousdMacroGdMacroCohesion,
                                           previousdMicroGdMicroCohesion, previousdMicroGradientGdMicroGradientCohesion );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model (jacobian)",
                                                 "Error in the computation of the previous plastic flow directions" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "previousMacroFlowDirection", previousMacroFlowDirection );
            tempDEBUG.emplace( "previousMicroFlowDirection", previousMicroFlowDirection );
            tempDEBUG.emplace( "previousMicroGradientFlowDirection", previousMicroGradientFlowDirection );

            tmp = { previousdMacroGdMacroCohesion };
            tempDEBUG.emplace( "previousdMacroGdMacroCohesion", tmp );

            tmp = { previousdMicroGdMicroCohesion };
            tempDEBUG.emplace( "previousdMicroGdMicroCohesion", tmp );

            tempDEBUG.emplace( "previousdMicroGradientGdMicroGradientCohesion",
                               tardigradeVectorTools::appendVectors( previousdMicroGradientGdMicroGradientCohesion ) );
#endif

            //Update the strain ISV values
            error = evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma, currentMicroGradientGamma,
                                                currentdMacroGdMacroC, currentdMicroGdMicroC, currentdMicroGradientGdMicroGradientC,
                                                previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
                                                previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
                                                previousdMacroGdMacroCohesion, previousdMicroGdMicroCohesion,
                                                previousdMicroGradientGdMicroGradientCohesion,
                                                currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                                alphaMacro, alphaMicro, alphaMicroGradient );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model (jacobian)",
                                                 "Error in the evolution of the current strain state variables" );
                result->addNext( error );
                result->print();               //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tmp = { currentMacroStrainISV };
            tempDEBUG.emplace( "currentMacroStrainISV", tmp );

            tmp = { currentMicroStrainISV };
            tempDEBUG.emplace( "currentMicroStrainISV", tmp );

            tempDEBUG.emplace( "currentMicroGradientStrainISV", currentMicroGradientStrainISV );
#endif

            //Compute the current cohesion values
            error = computeCohesion( currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                     macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters,
                                     currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model (jacobian)",
                                                 "Error in the evolution of the updated cohesion" );
                result->addNext( error );
                result->print();               //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tmp = { currentMacroCohesion };
            tempDEBUG.emplace( "currentMacroCohesion", tmp );

            tmp = { currentMicroCohesion };
            tempDEBUG.emplace( "currentMicroCohesion", tmp );

            tempDEBUG.emplace( "currentMicroGradientCohesion", currentMicroGradientCohesion );
#endif

            //Compute the previous plastic velocity gradients
            error = computePlasticVelocityGradients( previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
                                                     previousElasticRightCauchyGreen, previousElasticMicroRightCauchyGreen,
                                                     previousElasticPsi, previousElasticGamma, previousMacroFlowDirection,
                                                     previousMicroFlowDirection, previousMicroGradientFlowDirection,
                                                     previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                                     previousPlasticMicroGradientVelocityGradient );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model (jacobian)",
                                                 "Error in the computation of the previous plastic velocity gradient" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "previousPlasticMacroVelocityGradient", previousPlasticMacroVelocityGradient );
            tempDEBUG.emplace( "previousPlasticMicroVelocityGradient", previousPlasticMicroVelocityGradient );
            tempDEBUG.emplace( "previousPlasticMicroGradientVelocityGradient", previousPlasticMicroGradientVelocityGradient );
#endif

            //Update the current plastic deformation measures
            variableVector currentPlasticMacroVelocityGradient( dim * dim, 0 );
            variableVector currentPlasticMicroVelocityGradient( dim * dim, 0 );
            variableVector currentPlasticMicroGradientVelocityGradient( dim * dim * dim, 0 );

            error = evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient, currentPlasticMicroVelocityGradient,
                                              currentPlasticMicroGradientVelocityGradient, previousPlasticDeformationGradient,
                                              previousPlasticMicroDeformation, previousPlasticGradientMicroDeformation,
                                              previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                              previousPlasticMicroGradientVelocityGradient,
                                              currentPlasticDeformationGradient, currentPlasticMicroDeformation,
                                              currentPlasticGradientMicroDeformation,
                                              alphaMacro, alphaMicro, alphaMicroGradient );

            if ( error ){
                errorOut result = new errorNode( "evaluate_model (jacobian)",
                                                 "Error in the computation of the current plastic deformation" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }
        }

        //Update the elastic deformation
        variableVector currentElasticDeformationGradient, currentElasticMicroDeformation, currentElasticGradientMicroDeformation;

        variableMatrix dElasticDeformationGradientdDeformationGradient, dElasticDeformationGradientdPlasticDeformationGradient,
                       dElasticMicroDeformationdMicroDeformation, dElasticMicroDeformationdPlasticMicroDeformation,
                       dElasticGradientMicroDeformationdGradientMicroDeformation,
                       dElasticGradientMicroDeformationdPlasticGradientMicroDeformation,
                       dElasticGradientMicroDeformationdPlasticDeformationGradient, dElasticGradientMicroDeformationdMicroDeformation,
                       dElasticGradientMicroDeformationdPlasticMicroDeformation;

        error = computeElasticPartOfDeformation( currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation,
                                                 currentPlasticDeformationGradient, currentPlasticMicroDeformation,
                                                 currentPlasticGradientMicroDeformation, currentElasticDeformationGradient,
                                                 currentElasticMicroDeformation, currentElasticGradientMicroDeformation,
                                                 dElasticDeformationGradientdDeformationGradient,
                                                 dElasticDeformationGradientdPlasticDeformationGradient,
                                                 dElasticMicroDeformationdMicroDeformation,
                                                 dElasticMicroDeformationdPlasticMicroDeformation,
                                                 dElasticGradientMicroDeformationdGradientMicroDeformation,
                                                 dElasticGradientMicroDeformationdPlasticGradientMicroDeformation,
                                                 dElasticGradientMicroDeformationdPlasticDeformationGradient,
                                                 dElasticGradientMicroDeformationdMicroDeformation,
                                                 dElasticGradientMicroDeformationdPlasticMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "Error in the computation of the current elastic deformation" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "currentElasticDeformationGradient", currentElasticDeformationGradient );
        tempDEBUG.emplace( "currentElasticMicroDeformation", currentElasticMicroDeformation );
        tempDEBUG.emplace( "currentElasticGradientMicroDeformation", currentElasticGradientMicroDeformation );

        tempDEBUG.emplace( "currentPlasticDeformationGradient", currentPlasticDeformationGradient );
        tempDEBUG.emplace( "currentPlasticMicroDeformation", currentPlasticMicroDeformation );
        tempDEBUG.emplace( "currentPlasticGradientMicroDeformation", currentPlasticGradientMicroDeformation );
#endif

        //Compute the right Cauchy-Green deformation gradient
        variableVector currentElasticRightCauchyGreen;

        error = tardigradeConstitutiveTools::computeRightCauchyGreen( currentElasticDeformationGradient, currentElasticRightCauchyGreen );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "Error in the computation of the current elastic right Cauchy-Green deformation" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        //Compute the new stress values
        variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

        variableMatrix dPK2StressdElasticDeformationGradient, dPK2StressdElasticMicroDeformation,
                       dPK2StressdElasticGradientMicroDeformation,
                       dReferenceMicroStressdElasticDeformationGradient, dReferenceMicroStressdElasticMicroDeformation,
                       dReferenceMicroStressdElasticGradientMicroDeformation,
                       dReferenceHigherOrderStressdElasticDeformationGradient,
                       dReferenceHigherOrderStressdElasticGradientMicroDeformation;

        error = tardigradeMicromorphicLinearElasticity::linearElasticityReference( currentElasticDeformationGradient,
                                                                         currentElasticMicroDeformation,
                                                                         currentElasticGradientMicroDeformation,
                                                                         Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                         currentPK2Stress, currentReferenceMicroStress,
                                                                         currentReferenceHigherOrderStress,
                                                                         dPK2StressdElasticDeformationGradient,
                                                                         dPK2StressdElasticMicroDeformation,
                                                                         dPK2StressdElasticGradientMicroDeformation,
                                                                         dReferenceMicroStressdElasticDeformationGradient,
                                                                         dReferenceMicroStressdElasticMicroDeformation,
                                                                         dReferenceMicroStressdElasticGradientMicroDeformation,
                                                                         dReferenceHigherOrderStressdElasticDeformationGradient,
                                                                         dReferenceHigherOrderStressdElasticGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model (jacobian)",
                                             "Error in the computation of the current stress measures" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "intermediatePK2Stress", currentPK2Stress );
        tempDEBUG.emplace( "intermediateReferenceMicroStress", currentReferenceMicroStress );
        tempDEBUG.emplace( "intermediateReferenceHigherOrderStress", currentReferenceHigherOrderStress );
#endif

        //Evaluate the yield functions
        variableVector currentYieldFunctionValues;

        error = evaluateYieldFunctions( currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress,
                                        currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion,
                                        currentElasticRightCauchyGreen,
                                        macroYieldParameters, microYieldParameters, microGradientYieldParameters,
                                        currentYieldFunctionValues
#ifdef DEBUG_MODE
                                        , tempDEBUG
#endif
                                       );

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "currentYieldFunctionValues", currentYieldFunctionValues );
        tempITERATION.emplace( "pre_iteration_values", tempDEBUG );
        DEBUG.emplace( "pre_iteration_values", tempITERATION );
        tempITERATION.clear();
        tempDEBUG.clear();
#endif

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in the computation of the current yield function values" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        /*============================
        | Begin the non-linear solve |
        ============================*/

        //Check if any of the surfaces are yielding and begin the non-linear solver if they are
        tardigradeSolverTools::floatVector solutionVector;
        tardigradeSolverTools::intVector activePlasticity( currentYieldFunctionValues.size(), 0 );
        bool convergenceFlag = true;

        for ( unsigned int i = 0; i < currentYieldFunctionValues[ i ]; i++ ){
            if ( currentYieldFunctionValues[ i ] > absoluteTolerance ){

                convergenceFlag = false;
                
                activePlasticity[ i ] = 1;
            }
        }

        if ( !convergenceFlag ){
            tardigradeSolverTools::floatMatrix floatArgs =
                {
                    { Dt },
                    currentDeformationGradient,
                    currentMicroDeformation,
                    currentGradientMicroDeformation,
                    { previousMacroGamma },
                    { previousMicroGamma },
                    previousMicroGradientGamma,
                    previousPlasticDeformationGradient,
                    previousPlasticMicroDeformation,
                    previousPlasticGradientMicroDeformation,
                    { previousMacroStrainISV },
                    { previousMicroStrainISV },
                    previousMicroGradientStrainISV,
                    { previousdMacroGdMacroCohesion },
                    { previousdMicroGdMicroCohesion },
                    tardigradeVectorTools::appendVectors( previousdMicroGradientGdMicroGradientCohesion ),
                    previousPlasticMacroVelocityGradient,
                    previousPlasticMicroVelocityGradient,
                    previousPlasticMicroGradientVelocityGradient,
                    macroFlowParameters,
                    microFlowParameters,
                    microGradientFlowParameters,
                    macroHardeningParameters,
                    microHardeningParameters,
                    microGradientHardeningParameters,
                    macroYieldParameters,
                    microYieldParameters,
                    microGradientYieldParameters,
                    Amatrix,
                    Bmatrix,
                    Cmatrix,
                    Dmatrix,
                    { alphaMacro },
                    { alphaMicro },
                    { alphaMicroGradient }
                };

            tardigradeSolverTools::floatVector dPlasticDeformationdSolutionVector, dStressdSolutionVector,
                                     dPlasticDeformationdDeformation, dStressdDeformation,
                                     dResidualdDeformation;

            tardigradeSolverTools::floatMatrix floatOuts =
                {
                    currentPK2Stress,
                    currentReferenceMicroStress,
                    currentReferenceHigherOrderStress,
                    { currentMacroStrainISV },
                    { currentMicroStrainISV },
                    currentMicroGradientStrainISV,
                    {}, //The plastic deformation gradient
                    {}, //The plastic micro deformation
                    {}, //The gradient of plastic micro deformation
                    dPlasticDeformationdSolutionVector,
                    dStressdSolutionVector,
                    dPlasticDeformationdDeformation,
                    dStressdDeformation,
                    dResidualdDeformation
                };

            tardigradeSolverTools::intMatrix intOuts = { { } };

            tardigradeSolverTools::stdFncNLFJ func
                = static_cast< tardigradeSolverTools::NonLinearFunctionWithJacobian >( computePlasticMultiplierResidual );

            tardigradeSolverTools::floatVector x0 =
                {
                    currentMacroGamma,
                    currentMicroGamma,
                    currentMicroGradientGamma[ 0 ],
                    currentMicroGradientGamma[ 1 ],
                    currentMicroGradientGamma[ 2 ]
                };

            tardigradeSolverTools::floatVector solutionVector;

            bool convergeFlag, fatalErrorFlag;
           
            tardigradeSolverTools::intMatrix intArgs = { { 1 } };

            tardigradeSolverTools::solverType linearSolver;
            tardigradeSolverTools::floatMatrix J;
            tardigradeSolverTools::intVector boundVariableIndices = {  0,  1,  2,  3,  4 };
            tardigradeSolverTools::intVector boundSigns           = {  0,  0,  0,  0,  0 };
            tardigradeSolverTools::floatVector boundValues        = {  0,  0,  0,  0,  0 };

            error = tardigradeSolverTools::homotopySolver( func, x0, solutionVector, convergeFlag, fatalErrorFlag,
                                                 floatOuts, intOuts, floatArgs, intArgs, linearSolver, J,
                                                 boundVariableIndices, boundSigns, boundValues, true,
#ifdef DEBUG_MODE
                                                 DEBUG,
#endif
                                                 20, relativeTolerance, absoluteTolerance,
                                                 1e-4, 5, 1.0, .01 );

            if ( ( error ) && fatalErrorFlag ){ //Fatal error
                errorOut result = new errorNode( "evaluate_model",
                                                 "Error in nonlinear solve" );
                result->addNext( error );
                result->print();           //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing
                return 2;
            }
            else if ( ( error ) && !convergeFlag ){ //Solution didn't converge
                errorOut result = new errorNode( "evaluate_model",
                                                 "Solution failed to converge. Requesting timestep cutback" );
                result->addNext( error );
                result->print();          //Print the error message
                output_message = buffer.str(); //Save the output to enable message passing

                //Additional messages
                output_message += "\nx0\n{ ";
                for ( unsigned int kk = 0; kk < x0.size(); kk++ ){
                    output_message += std::to_string( x0[ kk ] ) + ", ";
                }
                output_message += "} \n";

                output_message += "\nfloatArgs\n";
                for ( unsigned int kk = 0; kk < floatArgs.size(); kk++ ){
                    output_message += "{ ";
                    for ( unsigned int ll = 0; ll < floatArgs[ kk ].size(); ll++ ){
                        output_message += std::to_string( floatArgs[ kk ][ ll ] ) + ", ";
                    }
                    output_message += "},\n";
                }

                output_message += "\nintArgs\n";
                for ( unsigned int kk = 0; kk < intArgs.size(); kk++ ){
                    output_message += "{ ";
                    for ( unsigned int ll = 0; ll < intArgs[ kk ].size(); ll++ ){
                        output_message += std::to_string( intArgs[ kk ][ ll ] ) + ", ";
                    }
                    output_message += "}\n";
                }
                
                return 1;
            }

            //Extract the deformation measures
            currentPlasticDeformationGradient      = floatOuts[ 6 ];
            currentPlasticMicroDeformation         = floatOuts[ 7 ];
            currentPlasticGradientMicroDeformation = floatOuts[ 8 ];
            currentMacroGamma                      = solutionVector[ 0 ];
            currentMicroGamma                      = solutionVector[ 1 ];
            currentMicroGradientGamma              = variableVector( solutionVector.begin() + 2, solutionVector.begin() + 5 );

            if ( currentMacroGamma < -absoluteTolerance ){
                errorOut result = new errorNode( "evaluate_model", "The macro-scale plastic multiplier is negative" );
                result->print();
                output_message = buffer.str();
                return 1;
            }

            if ( currentMicroGamma < -absoluteTolerance ){
                errorOut result = new errorNode( "evaluate_model", "The micro-scale plastic multiplier is negative" );
                result->print();
                output_message = buffer.str();
                return 1;
            }

            for ( unsigned int i = 0; i < currentMicroGradientGamma.size(); i++ ){
                if ( currentMicroGradientGamma[ i ] < -absoluteTolerance ){
                    errorOut result = new errorNode( "evaluate_model", "A micro-gradient plastic multiplier is negative" );
                    result->print();
                    output_message = buffer.str();
                    return 1;
                }
            }

            //Extract the stresses and state variables
            currentPK2Stress                       = floatOuts[  0 ];
            currentReferenceMicroStress            = floatOuts[  1 ];
            currentReferenceHigherOrderStress      = floatOuts[  2 ];
            currentMacroStrainISV                  = floatOuts[  3 ][ 0 ];
            currentMicroStrainISV                  = floatOuts[  4 ][ 0 ];
            currentMicroGradientStrainISV          = floatOuts[  5 ];

            //Extract the stresses and Jacobians
            dPlasticDeformationdSolutionVector     = floatOuts[  9 ];
            dStressdSolutionVector                 = floatOuts[ 10 ];
            dPlasticDeformationdDeformation        = floatOuts[ 11 ];
            dStressdDeformation                    = floatOuts[ 12 ];
            dResidualdDeformation                  = floatOuts[ 13 ];

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "convergedPlasticDeformationGradient", currentPlasticDeformationGradient );
            tempDEBUG.emplace( "convergedPlasticMicroDeformation", currentPlasticMicroDeformation );
            tempDEBUG.emplace( "convergedPlasticGradientMicroDeformation", currentPlasticGradientMicroDeformation );

            tardigradeSolverTools::floatVector temp = { currentMacroStrainISV };
            tempDEBUG.emplace( "convergedMacroStrainISV", temp );

            temp = { currentMicroStrainISV };
            tempDEBUG.emplace( "convergedMicroStrainISV", temp );

            tempDEBUG.emplace( "convergedMicroGradientStrainISV", currentMicroGradientStrainISV );

            temp = { currentMacroGamma };
            tempDEBUG.emplace( "convergedMacroGamma", temp );

            temp = { currentMicroGamma };
            tempDEBUG.emplace( "convergedMicroGamma", temp );

            tempDEBUG.emplace( "convergedMicroGradientGamma", currentMicroGradientGamma );

            tempDEBUG.emplace( "intermediatePK2Stress", currentPK2Stress );
            tempDEBUG.emplace( "intermediateReferenceMicroStress", currentReferenceMicroStress );
            tempDEBUG.emplace( "intermediateReferenceHigherOrderStress", currentReferenceHigherOrderStress );
#endif

            //Form the Jacobian of the stresses w.r.t. the fundamental deformation measures

            //Compute decomposition of the final Jacobian
            tardigradeSolverTools::floatVector dResidualdSolutionVector = tardigradeVectorTools::appendVectors( J );
            Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                dRdX( dResidualdSolutionVector.data(), 5, 5 );

            linearSolver = tardigradeSolverTools::solverType( dRdX );

            Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                dRdD( dResidualdDeformation.data(), 5, 45 );

            Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                dSdX( dStressdSolutionVector.data(), 45, 5 );

            Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                dSdD( dStressdDeformation.data(), 45, 45 );

            Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                dPDdX( dPlasticDeformationdSolutionVector.data(), 45, 5 );

            Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                dPDdD( dPlasticDeformationdDeformation.data(), 45, 45 );

            //Solve for the total derivative of the residual w.r.t. the deformation 
            tardigradeSolverTools::floatVector DSolutionVectorDDeformation( 5 * 45, 0 );
            Eigen::Map< Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                DXDD( DSolutionVectorDDeformation.data(), 5, 45 );

            //Perform the linear solve
            DXDD = -linearSolver.solve( dRdD );

            tardigradeSolverTools::floatVector DStressDDeformation( 45 * 45, 0 );
            Eigen::Map< Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                DSDD( DStressDDeformation.data(), 45, 45 );

            DSDD = dSdD + dSdX * DXDD;

            tardigradeSolverTools::floatVector DPlasticDeformationDDeformation( 45 * 45, 0 );
            Eigen::Map< Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                DPDDD( DPlasticDeformationDDeformation.data(), 45, 45 );

            DPDDD = dPDdD + dPDdX * DXDD;

            //Extract the total derivatives of the plastic deformation measures and stresses

            for ( unsigned int i = 0; i < 9; i++ ){
                for ( unsigned int j = 0; j < 9; j++ ){

                    //The plastic deformation measures
                    dPlasticDeformationGradientdDeformationGradient[ i ][ j ] = DPlasticDeformationDDeformation[ 45 * i + j ];
                    dPlasticDeformationGradientdMicroDeformation[ i ][ j ]    = DPlasticDeformationDDeformation[ 45 * i + j + 9 ];

                    dPlasticMicroDeformationdDeformationGradient[ i ][ j ]    = DPlasticDeformationDDeformation[ 45 * ( i + 9 ) + j ];
                    dPlasticMicroDeformationdMicroDeformation[ i ][ j ]       = DPlasticDeformationDDeformation[ 45 * ( i + 9 ) + j + 9 ];

                    //The stress measures
                    dPK2StressdDeformationGradient[ i ][ j ]            = DStressDDeformation[ 45 * i + j ];
                    dPK2StressdMicroDeformation[ i ][ j ]               = DStressDDeformation[ 45 * i + j + 9 ];

                    dReferenceMicroStressdDeformationGradient[ i ][ j ] = DStressDDeformation[ 45 * ( i + 9 ) + j ];
                    dReferenceMicroStressdMicroDeformation[ i ][ j ]    = DStressDDeformation[ 45 * ( i + 9 ) + j + 9 ];
                }

                for ( unsigned int j = 0; j < 27; j++ ){

                    //The plastic deformation measures
                    dPlasticDeformationGradientdGradientMicroDeformation[ i ][ j ]
                        = DPlasticDeformationDDeformation[ 45 * i + j + 18 ];
                    dPlasticMicroDeformationdGradientMicroDeformation[ i ][ j ]
                        = DPlasticDeformationDDeformation[ 45 * ( i + 9 ) + j + 18 ];

                    //The stress measures
                    dPK2StressdGradientMicroDeformation[ i ][ j ]            = DStressDDeformation[ 45 * i + j + 18 ];
                    dReferenceMicroStressdGradientMicroDeformation[ i ][ j ] = DStressDDeformation[ 45 * ( i + 9 ) + j + 18 ];
                }
            }

            for ( unsigned int i = 0; i < 27; i++ ){
                for ( unsigned int j = 0; j < 9; j++ ){

                    //The plastic deformation measures
                    dPlasticGradientMicroDeformationdDeformationGradient[ i ][ j ]
                        = DPlasticDeformationDDeformation[ 45 * ( i + 18 ) + j ];
                    dPlasticGradientMicroDeformationdMicroDeformation[ i ][ j ]
                        = DPlasticDeformationDDeformation[ 45 * ( i + 18 ) + j + 9 ];

                    //The stress measures
                    dReferenceHigherOrderStressdDeformationGradient[ i ][ j ]
                        = DStressDDeformation[ 45 * ( i + 18 ) + j ];
                    dReferenceHigherOrderStressdMicroDeformation[ i ][ j ]
                        = DStressDDeformation[ 45 * ( i + 18 ) + j + 9 ];
                }

                for ( unsigned int j = 0; j < 27; j++ ){

                    //The plastic deformation measures
                    dPlasticGradientMicroDeformationdGradientMicroDeformation[ i ][ j ]
                        = DPlasticDeformationDDeformation[ 45 * ( i + 18 ) + j + 18 ];

                    //The stress measures
                    dReferenceHigherOrderStressdGradientMicroDeformation[ i ][ j ]
                        = DStressDDeformation[ 45 * ( i + 18 ) + j + 18 ];
                }
            }
        }
        else{

            dPK2StressdDeformationGradient = tardigradeVectorTools::dot( dPK2StressdElasticDeformationGradient,
                                                               dElasticDeformationGradientdDeformationGradient );
            dPK2StressdMicroDeformation = tardigradeVectorTools::dot( dPK2StressdElasticMicroDeformation,
                                                            dElasticMicroDeformationdMicroDeformation )
                                        + tardigradeVectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                                            dElasticGradientMicroDeformationdMicroDeformation );
            dPK2StressdGradientMicroDeformation = tardigradeVectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            dReferenceMicroStressdDeformationGradient = tardigradeVectorTools::dot( dReferenceMicroStressdElasticDeformationGradient,
                                                                          dElasticDeformationGradientdDeformationGradient );
            dReferenceMicroStressdMicroDeformation = tardigradeVectorTools::dot( dReferenceMicroStressdElasticMicroDeformation,
                                                                       dElasticMicroDeformationdMicroDeformation )
                                                   + tardigradeVectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                                                       dElasticGradientMicroDeformationdMicroDeformation );
            dReferenceMicroStressdGradientMicroDeformation
                = tardigradeVectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            dReferenceHigherOrderStressdDeformationGradient = tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticDeformationGradient,
                                                                                dElasticDeformationGradientdDeformationGradient );
            dReferenceHigherOrderStressdMicroDeformation = tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                                                             dElasticGradientMicroDeformationdMicroDeformation );
            dReferenceHigherOrderStressdGradientMicroDeformation
                = tardigradeVectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdGradientMicroDeformation );


        }

#ifdef DEBUG_MODE
        //Save the total derivatives of the plastic deformation measures w.r.t. the fundamental deformation measures
        tempDEBUG.emplace( "totaldPlasticDeformationGradientdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdDeformationGradient ) );
        tempDEBUG.emplace( "totaldPlasticDeformationGradientdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdMicroDeformation ) );
        tempDEBUG.emplace( "totaldPlasticDeformationGradientdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdGradientMicroDeformation ) );

        tempDEBUG.emplace( "totaldPlasticMicroDeformationdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dPlasticMicroDeformationdDeformationGradient ) );
        tempDEBUG.emplace( "totaldPlasticMicroDeformationdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticMicroDeformationdMicroDeformation ) );
        tempDEBUG.emplace( "totaldPlasticMicroDeformationdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticMicroDeformationdGradientMicroDeformation ) );

        tempDEBUG.emplace( "totaldPlasticGradientMicroDeformationdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dPlasticGradientMicroDeformationdDeformationGradient ) );
        tempDEBUG.emplace( "totaldPlasticGradientMicroDeformationdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticGradientMicroDeformationdMicroDeformation ) );
        tempDEBUG.emplace( "totaldPlasticGradientMicroDeformationdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPlasticGradientMicroDeformationdGradientMicroDeformation ) );

        //Save the total derivatives of the intermediate stresses w.r.t. the fundamental deformation measures
        tempDEBUG.emplace( "totaldPK2StressdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dPK2StressdDeformationGradient ) );
        tempDEBUG.emplace( "totaldPK2StressdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPK2StressdMicroDeformation ) );
        tempDEBUG.emplace( "totaldPK2StressdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dPK2StressdGradientMicroDeformation ) );

        tempDEBUG.emplace( "totaldReferenceMicroStressdDeformationGradient",
                            tardigradeVectorTools::appendVectors( dReferenceMicroStressdDeformationGradient ) );
        tempDEBUG.emplace( "totaldReferenceMicroStressdMicroDeformation",
                            tardigradeVectorTools::appendVectors( dReferenceMicroStressdMicroDeformation ) );
        tempDEBUG.emplace( "totaldReferenceMicroStressdGradientMicroDeformation",
                            tardigradeVectorTools::appendVectors( dReferenceMicroStressdGradientMicroDeformation ) );

        tempDEBUG.emplace( "totaldReferenceHigherOrderStressdDeformationGradient",
                           tardigradeVectorTools::appendVectors( dReferenceHigherOrderStressdDeformationGradient ) );
        tempDEBUG.emplace( "totaldReferenceHigherOrderStressdMicroDeformation",
                           tardigradeVectorTools::appendVectors( dReferenceHigherOrderStressdMicroDeformation ) );
        tempDEBUG.emplace( "totaldReferenceHigherOrderStressdGradientMicroDeformation",
                           tardigradeVectorTools::appendVectors( dReferenceHigherOrderStressdGradientMicroDeformation ) );

        tempITERATION.emplace( "converged_values", tempDEBUG );
        DEBUG.emplace( "converged_values", tempITERATION );
        tempITERATION.clear();
        tempDEBUG.clear();
#endif

        /*============================
        | Assemble the output values |
        ============================*/
        
        //Assemble the SDV vector
        std::vector< double > currentStrainISVS = { currentMacroStrainISV, currentMicroStrainISV };
        currentStrainISVS = tardigradeVectorTools::appendVectors( { currentStrainISVS, currentMicroGradientStrainISV } );

        std::vector< double > currentGammas = { currentMacroGamma, currentMicroGamma };
        currentGammas = tardigradeVectorTools::appendVectors( { currentGammas, currentMicroGradientGamma } );

        SDVS = tardigradeVectorTools::appendVectors( { currentStrainISVS, currentGammas, currentPlasticDeformationGradient - eye,
                                             currentPlasticMicroDeformation - eye, currentPlasticGradientMicroDeformation } );

        //Output the stresses mapped to the true reference configuration
        variableMatrix dPK2dIntermediatePK2, dPK2dPlasticDeformationGradient;
        error = tardigradeMicromorphicTools::pullBackCauchyStress( currentPK2Stress, currentPlasticDeformationGradient, current_PK2,
                                                         dPK2dIntermediatePK2, dPK2dPlasticDeformationGradient );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in pullback operation on the PK2 stress" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        //Assemble the Jacobians
        variableMatrix dReferencePK2StressdDeformationGradient
            = tardigradeVectorTools::dot( dPK2dIntermediatePK2, dPK2StressdDeformationGradient )
            + tardigradeVectorTools::dot( dPK2dPlasticDeformationGradient, dPlasticDeformationGradientdDeformationGradient );

        variableMatrix dReferencePK2StressdMicroDeformation
            = tardigradeVectorTools::dot( dPK2dIntermediatePK2, dPK2StressdMicroDeformation )
            + tardigradeVectorTools::dot( dPK2dPlasticDeformationGradient, dPlasticDeformationGradientdMicroDeformation );

        variableMatrix dReferencePK2StressdGradientMicroDeformation
            = tardigradeVectorTools::dot( dPK2dIntermediatePK2, dPK2StressdGradientMicroDeformation )
            + tardigradeVectorTools::dot( dPK2dPlasticDeformationGradient, dPlasticDeformationGradientdGradientMicroDeformation );

        variableMatrix dSIGMAdIntermediateSIGMA, dSIGMAdPlasticDeformationGradient;
        error = tardigradeMicromorphicTools::pullBackMicroStress( currentReferenceMicroStress, currentPlasticDeformationGradient, current_SIGMA,
                                                         dSIGMAdIntermediateSIGMA, dSIGMAdPlasticDeformationGradient );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in pullback operation on the reference symmetric micro-stress" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        //Assemble the Jacobians
        variableMatrix dSIGMAStressdDeformationGradient
            = tardigradeVectorTools::dot( dSIGMAdIntermediateSIGMA, dReferenceMicroStressdDeformationGradient )
            + tardigradeVectorTools::dot( dSIGMAdPlasticDeformationGradient, dPlasticDeformationGradientdDeformationGradient );

        variableMatrix dSIGMAStressdMicroDeformation
            = tardigradeVectorTools::dot( dSIGMAdIntermediateSIGMA, dReferenceMicroStressdMicroDeformation )
            + tardigradeVectorTools::dot( dSIGMAdPlasticDeformationGradient, dPlasticDeformationGradientdMicroDeformation );

        variableMatrix dSIGMAStressdGradientMicroDeformation
            = tardigradeVectorTools::dot( dSIGMAdIntermediateSIGMA, dReferenceMicroStressdGradientMicroDeformation )
            + tardigradeVectorTools::dot( dSIGMAdPlasticDeformationGradient, dPlasticDeformationGradientdGradientMicroDeformation );

        variableMatrix dMdIntermediateM, dMdPlasticDeformationGradient, dMdPlasticMicroDeformation;
        error = tardigradeMicromorphicTools::pullBackHigherOrderStress( currentReferenceHigherOrderStress, 
                                                              currentPlasticDeformationGradient, 
                                                              currentPlasticMicroDeformation, current_M,
                                                              dMdIntermediateM, dMdPlasticDeformationGradient,
                                                              dMdPlasticMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in pullback operation on the reference higher order stress" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        //Assemble the Jacobians
        variableMatrix dMStressdDeformationGradient
            = tardigradeVectorTools::dot( dMdIntermediateM, dReferenceHigherOrderStressdDeformationGradient )
            + tardigradeVectorTools::dot( dMdPlasticDeformationGradient, dPlasticDeformationGradientdDeformationGradient )
            + tardigradeVectorTools::dot( dMdPlasticMicroDeformation, dPlasticMicroDeformationdDeformationGradient );

        variableMatrix dMStressdMicroDeformation
            = tardigradeVectorTools::dot( dMdIntermediateM, dReferenceHigherOrderStressdMicroDeformation )
            + tardigradeVectorTools::dot( dMdPlasticDeformationGradient, dPlasticDeformationGradientdMicroDeformation )
            + tardigradeVectorTools::dot( dMdPlasticMicroDeformation, dPlasticMicroDeformationdMicroDeformation );

        variableMatrix dMStressdGradientMicroDeformation
            = tardigradeVectorTools::dot( dMdIntermediateM, dReferenceHigherOrderStressdGradientMicroDeformation )
            + tardigradeVectorTools::dot( dMdPlasticDeformationGradient, dPlasticDeformationGradientdGradientMicroDeformation )
            + tardigradeVectorTools::dot( dMdPlasticMicroDeformation, dPlasticMicroDeformationdGradientMicroDeformation );

        //Assemble the Jacobians w.r.t. the degrees of freedom
        DPK2Dgrad_u   = tardigradeVectorTools::dot( dReferencePK2StressdDeformationGradient,
                                          dDeformationGradientdGradientMacroDisplacement );
        DPK2Dphi      = tardigradeVectorTools::dot( dReferencePK2StressdMicroDeformation,
                                          dMicroDeformationdMicroDisplacement );
        DPK2Dgrad_phi = tardigradeVectorTools::dot( dReferencePK2StressdGradientMicroDeformation,
                                          dGradientMicroDeformationdGradientMicroDisplacement );

        DSIGMADgrad_u   = tardigradeVectorTools::dot( dSIGMAStressdDeformationGradient, dDeformationGradientdGradientMacroDisplacement );
        DSIGMADphi      = tardigradeVectorTools::dot( dSIGMAStressdMicroDeformation, dMicroDeformationdMicroDisplacement );
        DSIGMADgrad_phi = tardigradeVectorTools::dot( dSIGMAStressdGradientMicroDeformation, dGradientMicroDeformationdGradientMicroDisplacement );

        DMDgrad_u   = tardigradeVectorTools::dot( dMStressdDeformationGradient, dDeformationGradientdGradientMacroDisplacement );
        DMDphi      = tardigradeVectorTools::dot( dMStressdMicroDeformation, dMicroDeformationdMicroDisplacement );
        DMDgrad_phi = tardigradeVectorTools::dot( dMStressdGradientMicroDeformation, dGradientMicroDeformationdGradientMicroDisplacement );

        //Assemble the additional Jacobian terms
        ADD_JACOBIANS.resize(9);
        // Macro plastic deformation gradient terms
        ADD_JACOBIANS[0] = tardigradeVectorTools::dot( dPlasticDeformationGradientdDeformationGradient, dDeformationGradientdGradientMacroDisplacement );
        ADD_JACOBIANS[1] = tardigradeVectorTools::dot( dPlasticDeformationGradientdMicroDeformation, dMicroDeformationdMicroDisplacement );
        ADD_JACOBIANS[2] = tardigradeVectorTools::dot( dPlasticDeformationGradientdGradientMicroDeformation, dGradientMicroDeformationdGradientMicroDisplacement );

        // Micro plastic displacement terms
        ADD_JACOBIANS[3] = tardigradeVectorTools::dot( dPlasticMicroDeformationdDeformationGradient, dDeformationGradientdGradientMacroDisplacement );
        ADD_JACOBIANS[4] = tardigradeVectorTools::dot( dPlasticMicroDeformationdMicroDeformation, dMicroDeformationdMicroDisplacement );
        ADD_JACOBIANS[5] = tardigradeVectorTools::dot( dPlasticMicroDeformationdGradientMicroDeformation, dGradientMicroDeformationdGradientMicroDisplacement );

        // Micro plastic displacement terms
        ADD_JACOBIANS[6] = tardigradeVectorTools::dot( dPlasticGradientMicroDeformationdDeformationGradient, dDeformationGradientdGradientMacroDisplacement );
        ADD_JACOBIANS[7] = tardigradeVectorTools::dot( dPlasticGradientMicroDeformationdMicroDeformation, dMicroDeformationdMicroDisplacement );
        ADD_JACOBIANS[8] = tardigradeVectorTools::dot( dPlasticGradientMicroDeformationdGradientMicroDeformation, dGradientMicroDeformationdGradientMicroDisplacement );

        //Model evaluation successful. Return.
        return 0;
    }


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
                                        constantType &relativeTolerance, constantType &absoluteTolerance ){
        /*!
         * Extract the parameters from the parameter vector
         *
         * :param const std::vector< double > &fparams: The incoming parameter vector
         * :param parameterVector &macroHardeningParameters: The parameters used in the hardening of the macro Strain ISV
         * :param parameterVector &microHardeningParameters: The parameters used in the hardening of the micro Strain ISV
         * :param parameterVector &microGradientHardeningParameters: The parameters used in the hardening of the micro Gradient Strain ISV
         * :param parameterVector &macroFlowParameters: The parameters used in the macro flow direction computation.
         * :param parameterVector &microFlowParameters: The parameters used in the micro flow direction computation
         * :param parameterVector &microGradientFlowParameters: The parameters used in the micro Gradient flow direction computation.
         * :param parameterVector &macroYieldParameters: The parameters used in the macro yielding computation.
         * :param parameterVector &microYieldParameters: The parameters used in the micro yielding computation
         * :param parameterVector &microGradientYieldParameters: The parameters used in the micro Gradient yielding computation.
         * :param parameterVector &Amatrix: The A stiffness matrix.
         * :param parameterVector &Bmatrix: The B stiffness matrix.
         * :param parameterVector &Cmatrix: The C stiffness matrix.
         * :param parameterVector &Dmatrix: The D stiffness matrix.
         * :param parameterVector &alphaMacro: The integration parameter for the macro plasticity.
         * :param parameterVector &alphaMicro: The integration parameter for the micro plasticity.
         * :param parameterVector &alphaMicroGradient: The integration parameter for the micro gradient plasticity.
         * :param constantType &relativeTolerance: The relative tolerance for the solver.
         * :param constantType &absoluteTolerance: The absolute tolerance for the solver.
         */

        if ( fparams.size() == 0 ){
            return new errorNode( "extractMaterialParameters",
                                  "The material parameters vector has a length of 0" );
        }

        unsigned int start = 0;
        unsigned int span;

        std::vector< parameterVector > outputs( 13 );

        //Extract the material parameters
        for ( unsigned int i = 0; i < outputs.size(); i++ ){
            span = ( unsigned int )std::floor( fparams[ start ]  + 0.5 ); //Extract the span of the parameter set

            if ( fparams.size() < start + 1 + span ){
                std::string outstr = "fparams is not long enough to contain all of the required parameters:\n";
                outstr +=            "    filling variable " + std::to_string( i ) + "\n";
                outstr +=            "    size =          "  + std::to_string( fparams.size() ) + "\n";
                outstr +=            "    required size = "  + std::to_string( start + 1 + span );

                return new errorNode( "extractMaterialParameters",
                                      outstr.c_str() );
            }

            outputs[ i ] = parameterVector( fparams.begin() + start + 1, fparams.begin() + start + 1 + span );

            start = start + 1 + span;
        }

        //Set the output values
        macroHardeningParameters         = outputs[  0 ];
        microHardeningParameters         = outputs[  1 ];
        microGradientHardeningParameters = outputs[  2 ];
        macroFlowParameters              = outputs[  3 ];
        microFlowParameters              = outputs[  4 ];
        microGradientFlowParameters      = outputs[  5 ];
        macroYieldParameters             = outputs[  6 ];
        microYieldParameters             = outputs[  7 ];
        microGradientYieldParameters     = outputs[  8 ];

        //Form the stiffness tensors
        errorOut error;
        if ( outputs[ 9 ].size() == 2 ){
            error = tardigradeMicromorphicLinearElasticity::formIsotropicA( outputs[ 9 ][ 0 ], outputs[ 9 ][ 1 ], Amatrix );
        }
        else{
            std::string outstr = "Unrecognized number of parameters ( " + std::to_string( outputs[ 9 ].size() ) + " ) for the A stiffness tensor";
            return new errorNode( "extractMaterialParameters",
                                  outstr.c_str() );
        }
        
        if ( error ){
            errorOut result = new errorNode( "extractMaterialParameters", "Error in computation of the A stiffness tensor" );
            result->addNext( error );
            return result;
        }

        if ( outputs[ 10 ].size() == 5 ){
            error = tardigradeMicromorphicLinearElasticity::formIsotropicB( outputs[ 10 ][ 0 ], outputs[ 10 ][ 1 ], outputs[ 10 ][ 2 ],
                                                                 outputs[ 10 ][ 3 ], outputs[ 10 ][ 4 ], Bmatrix );
        }
        else{
            std::string outstr = "Unrecognized number of parameters ( " + std::to_string( outputs[ 10 ].size() ) + " ) for the B stiffness tensor";
            return new errorNode( "extractMaterialParameters",
                                  outstr.c_str() );
        }

        if ( error ){
            errorOut result = new errorNode( "extractMaterialParameters", "Error in computation of the B stiffness tensor" );
            result->addNext( error );
            return result;
        }

        if ( outputs[ 11 ].size() == 11 ){
            error = tardigradeMicromorphicLinearElasticity::formIsotropicC( outputs[ 11 ], Cmatrix );
        }
        else{
            std::string outstr = "Unrecognized number of parameters ( " + std::to_string( outputs[ 11 ].size() ) + " ) for the C stiffness tensor";
            return new errorNode( "extractMaterialParameters",
                                  outstr.c_str() );
        }

        if ( error ){
            errorOut result = new errorNode( "extractMaterialParameters", "Error in computation of the C stiffness tensor" );
            result->addNext( error );
            return result;
        }

        if ( outputs[ 12 ].size() == 2 ){
            error = tardigradeMicromorphicLinearElasticity::formIsotropicD( outputs[ 12 ][ 0 ], outputs[ 12 ][ 1 ], Dmatrix );
        }
        else{
            std::string outstr = "Unrecognized number of parameters ( " + std::to_string( outputs[ 12 ].size() ) + " ) for the D stiffness tensor";
            return new errorNode( "extractMaterialParameters",
                                  outstr.c_str() );
        }

        if ( error ){
            errorOut result = new errorNode( "extractMaterialParameters", "Error in computation of the D stiffness tensor" );
            result->addNext( error );
            return result;
        }

        //Extract the integration and tolerance parameters
        if ( fparams.size() < start + 5 ){
            return new errorNode( "extractMaterialParameters",
                                  "fparams does not store the integration parameters." );
        }

        alphaMacro         = fparams[ start + 0 ];
        alphaMicro         = fparams[ start + 1 ];
        alphaMicroGradient = fparams[ start + 2 ];
        relativeTolerance  = fparams[ start + 3 ];
        absoluteTolerance  = fparams[ start + 4 ];

        return NULL;
    }

    errorOut extractStateVariables( std::vector< double > &SDVS,
                                    variableType &previousMacroStrainISV, variableType &previousMicroStrainISV,
                                    variableVector &previousMicroGradientStrainISV,
                                    variableType &previousMacroGamma, variableType &previousMicroGamma,
                                    variableVector &previousMicroGradientGamma,
                                    variableVector &previousPlasticDeformationGradient,
                                    variableVector &previousPlasticMicroDeformation,
                                    variableVector &previousPlasticGradientMicroDeformation ){
        /*!
         * Extract the state variables from the state variable vector.
         *
         * :param std::vector< double > &SDVS: The state variable vector.
         * :param variableType &previousMacroStrainISV: The previous value of the macro Strain ISV
         * :param variableType &previousMicroStrainISV: The previous value of the micro Strain ISV
         * :param variableVector &previousMicroGradientStrainISV: The previous value of the micro gradient Strain ISV
         * :param variableType &previousMacroGamma: The previous value of the macro gamma
         * :param variableType &previousMicroGamma: The previous value of the micro gamma
         * :param variableVector &previousMicroGradientGamma: The previous value of the micro gradient gammas
         * :param variableVector &previousPlasticDeformationGradient: The previous value of the plastic deformation 
         *     gradient.
         * :param variableVector &previousPlasticMicroDeformation: The previous value of the plastic micro deformation.
         * :param variableVector &previousPlasticGradientMicroDeformation: The previous value of the plastic gradient of 
         *     the micro deformation.
         */

        //Assume 3D
        unsigned int dim = 3;
        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        if ( SDVS.size() != 55 ){
            std::string outstr = "The SDVS vector must have 55 elements ( has " + std::to_string( SDVS.size() ) + " )";
            return new errorNode( "extractStateVariables", outstr.c_str() );
        }

        previousMacroStrainISV = SDVS[ 0 ];
        previousMicroStrainISV = SDVS[ 1 ];
        previousMicroGradientStrainISV = { SDVS[ 2 ], SDVS[ 3 ], SDVS[ 4 ] };

        previousMacroGamma = SDVS[ 5 ];
        previousMicroGamma = SDVS[ 6 ];
        previousMicroGradientGamma = { SDVS[ 7 ], SDVS[ 8 ], SDVS[ 9 ] };

        previousPlasticDeformationGradient = eye + variableVector( SDVS.begin() + 10, SDVS.begin() + 19 );
        previousPlasticMicroDeformation    = eye + variableVector( SDVS.begin() + 19, SDVS.begin() + 28 );
        previousPlasticGradientMicroDeformation = variableVector( SDVS.begin() + 28, SDVS.begin() + 55 );

        return NULL;
    }

    errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation ){
        /*!
         * Assemble the fundamental deformation meaures from the degrees of freedom.
         *
         * :param const double ( &grad_u )[ 3 ][ 3 ]: The macro displacement gradient w.r.t. the reference configuration.
         * :param const double ( &phi )[ 9 ]: The micro displacement.
         * :param const double ( &grad_phi )[ 9 ][ 3 ]: The gradient of the micro displacement w.r.t. the reference configuration.
         * :param variableVector &deformationGradient: The deformation gradient
         * :param variableVector &microDeformation: The micro deformation
         * :param variableVector &gradientMicroDeformation: The gradient of the micro deformation.
         */


        //Extract the degrees of freedom
        variableMatrix displacementGradient = { { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] },
                                                { grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] },
                                                { grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] } };

        variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                             phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                             phi[ 6 ], phi[ 7 ], phi[ 8 ] };

        variableMatrix gradientMicroDisplacement = { { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] },
                                                     { grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] },
                                                     { grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] },
                                                     { grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] },
                                                     { grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] },
                                                     { grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] },
                                                     { grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] },
                                                     { grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] },
                                                     { grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] } };

        errorOut error = tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures",
                                             "Error in assembly of the deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures",
                                             "Error in assembly of the micro deformation" );
            result->addNext( error );
            return result;
        }

        error = tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures",
                                             "Error in assembly of the gradient of the micro deformation" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation, variableMatrix &dFdGradU,
                                                     variableMatrix &dChidPhi, variableMatrix &dGradChidGradPhi ){
        /*!
         * Assemble the fundamental deformation meaures from the degrees of freedom.
         *
         * :param const double ( &grad_u )[ 3 ][ 3 ]: The macro displacement gradient w.r.t. the reference configuration.
         * :param const double ( &phi )[ 9 ]: The micro displacement.
         * :param const double ( &grad_phi )[ 9 ][ 3 ]: The gradient of the micro displacement w.r.t. the reference configuration.
         * :param variableVector &deformationGradient: The deformation gradient
         * :param variableVector &microDeformation: The micro deformation
         * :param variableVector &gradientMicroDeformation: The gradient of the micro deformation.
         * :param variableMatrix &dFdGradU: The Jacobian of the deformation gradient w.r.t. the gradient of the displacement
         * :param variableMatrix &dChidPhi: The Jacobian of the micro deformation w.r.t. the micro displacement
         * :param variableMatrix &dGradChidGradPhi: The Jacobian of the gradient of the micro deformation w.r.t.
         *      the gradient of the micro displacement
         */


        //Extract the degrees of freedom
        variableMatrix displacementGradient = { { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] },
                                                { grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] },
                                                { grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] } };

        variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                             phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                             phi[ 6 ], phi[ 7 ], phi[ 8 ] };

        variableMatrix gradientMicroDisplacement = { { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] },
                                                     { grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] },
                                                     { grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] },
                                                     { grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] },
                                                     { grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] },
                                                     { grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] },
                                                     { grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] },
                                                     { grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] },
                                                     { grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] } };

        errorOut error = tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient, dFdGradU );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures (jacobian)",
                                             "Error in assembly of the deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation, dChidPhi );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures (jacobian)",
                                             "Error in assembly of the micro deformation" );
            result->addNext( error );
            return result;
        }

        error = tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation,
                                                                     dGradChidGradPhi );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures (jacobian)",
                                             "Error in assembly of the gradient of the micro deformation" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut evaluateYieldFunctions( const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                     const variableVector &referenceHigherOrderStress, const variableType &macroCohesion,
                                     const variableType &microCohesion, const variableVector &microGradientCohesion,
                                     const variableVector &elasticRightCauchyGreen,
                                     const parameterVector &macroYieldParameters, const parameterVector &microYieldParameters,
                                     const parameterVector &microGradientYieldParameters, variableVector &yieldFunctionValues
#ifdef DEBUG_MODE
                                     , tardigradeSolverTools::debugMap &DEBUG
#endif
                                    ){
        /*!
         * Evaluate all of the yield functions.
         *
         * :param const variableVector &PK2Stress: The second Piola Kirchhoff stress
         * :param const variableVector &referenceMicroStress: The micro stress in the reference configuration.
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableType &macroCohesion: The macro cohesion value.
         * :param const variableType &microCohesion: The micro cohesion value.
         * :param const variableVector &microGradientCohesion: The micro gradient cohesion value.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation tensor.
         * :param const parameterVector &macroYieldParameters: The macro yield parameters.
         * :param const parameterVector &microYieldParameters: The micro yield parameters.
         * :param const parameterVector &microGradientYieldParameters: The micro gradient yield parameters.
         * :param variableVector &yieldFunctionValues: The current values of the yield functions.
         * :param std::map< std::string, tardigradeSolverTools::floatVector > &DEBUG: The debug map. Only output when in debug mode.
         */

        if ( macroYieldParameters.size() != 2 ){
            std::cout << "macroYieldParameters: "; tardigradeVectorTools::print( macroYieldParameters );
            return new errorNode( "evaluateYieldFunctions",
                                  "The macro yield functions must have a length of 2" );
        }

        if ( microYieldParameters.size() != 2 ){
            return new errorNode( "evaluateYieldFunctions",
                                  "The micro yield functions must have a length of 2" );
        }

        if ( microGradientYieldParameters.size() != 2 ){
            return new errorNode( "evaluateYieldFunctions",
                                  "The micro Gradient yield functions must have a length of 2" );
        }

        yieldFunctionValues = variableVector( 5, 0 );
        errorOut error = computeSecondOrderDruckerPragerYieldEquation( PK2Stress, macroCohesion, elasticRightCauchyGreen,
                                                                       macroYieldParameters[ 0 ], macroYieldParameters[ 1 ],
                                                                       yieldFunctionValues[ 0 ] );

        if ( error ){
            errorOut result = new errorNode( "evaluateYieldFunctions",
                                             "Error in the computation of the macro yield equation" );
            result->addNext( error );
            return result;
        }

#ifdef DEBUG_MODE
        tardigradeSolverTools::floatVector tmp = { yieldFunctionValues[ 0 ] };
        DEBUG.emplace( "macroYieldFunction", tmp );
#endif

        error = computeSecondOrderDruckerPragerYieldEquation( referenceMicroStress, microCohesion,
                                                              elasticRightCauchyGreen,
                                                              microYieldParameters[ 0 ], microYieldParameters[ 1 ],
                                                              yieldFunctionValues[ 1 ] );

        if ( error ){
            errorOut result = new errorNode( "evaluateYieldFunctions",
                                             "Error in the computation of the micro yield equation" );
            result->addNext( error );
            return result;
        }

#ifdef DEBUG_MODE
        tmp = { yieldFunctionValues[ 1 ] };
        DEBUG.emplace( "microYieldFunction", tmp );
#endif

        variableVector yftmp;
        error = computeHigherOrderDruckerPragerYieldEquation( referenceHigherOrderStress, microGradientCohesion,
                                                              elasticRightCauchyGreen,
                                                              microGradientYieldParameters[ 0 ],
                                                              microGradientYieldParameters[ 1 ],
                                                              yftmp );

        if ( error ){
            errorOut result = new errorNode( "evaluateYieldFunctions",
                                             "Error in the computation of the micro yield equation" );
            result->addNext( error );
            return result;
        }

#ifdef DEBUG_MODE
        DEBUG.emplace( "microGradientYieldFunction", yftmp );
#endif

        yieldFunctionValues[ 2 ] = yftmp[ 0 ];
        yieldFunctionValues[ 3 ] = yftmp[ 1 ];
        yieldFunctionValues[ 4 ] = yftmp[ 2 ];

        return NULL;
    }

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
                                   ){
        /*!
         * Evaluate all of the yield functions.
         *
         * :param const variableVector &PK2Stress: The second Piola Kirchhoff stress
         * :param const variableVector &referenceMicroStress: The micro stress in the reference configuration.
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableType &macroCohesion: The macro cohesion value.
         * :param const variableType &microCohesion: The micro cohesion value.
         * :param const variableVector &microGradientCohesion: The micro gradient cohesion value.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation tensor.
         * :param const parameterVector &macroYieldParameters: The macro yield parameters.
         * :param const parameterVector &microYieldParameters: The micro yield parameters.
         * :param const parameterVector &microGradientYieldParameters: The micro gradient yield parameters.
         * :param variableVector &yieldFunctionValues: The current values of the yield functions.
         * :param variableVector &dMacroFdPK2: The Jacobian of the macro yield function w.r.t. the PK2 stress
         * :param variableType &dMacroFdMacroC: The Jacobian of the macro yield function w.r.t. the macro cohesion
         * :param variableVector &dMacroFdElasticRCG: The Jacobian of the macro yield function w.r.t. the elastic 
         *     right Cauchy-Green deformation tensor.
         * :param variableVector &dMicroFdSigma: The Jacobian of the micro yield function w.r.t. the reference
         *     symmetric micro stress.
         * :param variableType &dMicroFdMicroC: The Jacobian of the micro yield function w.r.t. the micro cohesion
         * :param variableVector &dMicroFdElasticRCG: The Jacobian of the micro yield function w.r.t. the elastic 
         *     right Cauchy-Green deformation tensor.
         * :param variableMatrix &dMicroGradientFdM: The Jacobian of the micro gradient yield function w.r.t. the reference
         *     higher order stress
         * :param variableMatrix &dMicroGradientFdMicroGradientC: The Jacobian of the micro gradient yield function w.r.t. 
         *     the micro gradient cohesion
         * :param variableMatrix &dMicroGradientFdElasticRCG: The Jacobian of the micro gradient yield function w.r.t. the elastic 
         *     right Cauchy-Green deformation tensor.
         * :param std::map< std::string, tardigradeSolverTools::floatVector > &DEBUG: The debug map. Only output when in debug mode.
         */

        if ( macroYieldParameters.size() != 2 ){
            std::cout << "macroYieldParameters: "; tardigradeVectorTools::print( macroYieldParameters );
            return new errorNode( "evaluateYieldFunctions (jacobian)",
                                  "The macro yield functions must have a length of 2" );
        }

        if ( microYieldParameters.size() != 2 ){
            return new errorNode( "evaluateYieldFunctions (jacobian)",
                                  "The micro yield functions must have a length of 2" );
        }

        if ( microGradientYieldParameters.size() != 2 ){
            return new errorNode( "evaluateYieldFunctions (jacobian)",
                                  "The micro Gradient yield functions must have a length of 2" );
        }

        yieldFunctionValues = variableVector( 5, 0 );
        errorOut error = computeSecondOrderDruckerPragerYieldEquation( PK2Stress, macroCohesion, elasticRightCauchyGreen,
                                                                       macroYieldParameters[ 0 ], macroYieldParameters[ 1 ],
                                                                       yieldFunctionValues[ 0 ], dMacroFdPK2, dMacroFdMacroC,
                                                                       dMacroFdElasticRCG );

        if ( error ){
            errorOut result = new errorNode( "evaluateYieldFunctions (jacobian)",
                                             "Error in the computation of the macro yield equation" );
            result->addNext( error );
            return result;
        }

#ifdef DEBUG_MODE
        tardigradeSolverTools::floatVector tmp = { yieldFunctionValues[ 0 ] };
            DEBUG.emplace( "macroYieldFunction", tmp );
#endif

        error = computeSecondOrderDruckerPragerYieldEquation( referenceMicroStress, microCohesion,
                                                              elasticRightCauchyGreen,
                                                              microYieldParameters[ 0 ], microYieldParameters[ 1 ],
                                                              yieldFunctionValues[ 1 ], dMicroFdSigma, dMicroFdMicroC,
                                                              dMicroFdElasticRCG );

        if ( error ){
            errorOut result = new errorNode( "evaluateYieldFunctions",
                                             "Error in the computation of the micro yield equation" );
            result->addNext( error );
            return result;
        }

#ifdef DEBUG_MODE
        tmp = { yieldFunctionValues[ 1 ] };
        DEBUG.emplace( "microYieldFunction", tmp );
#endif

        variableVector yftmp;
        error = computeHigherOrderDruckerPragerYieldEquation( referenceHigherOrderStress, microGradientCohesion,
                                                              elasticRightCauchyGreen,
                                                              microGradientYieldParameters[ 0 ],
                                                              microGradientYieldParameters[ 1 ],
                                                              yftmp, dMicroGradientFdM, dMicroGradientFdMicroGradientC,
                                                              dMicroGradientFdElasticRCG );

        if ( error ){
            errorOut result = new errorNode( "evaluateYieldFunctions",
                                             "Error in the computation of the micro yield equation" );
            result->addNext( error );
            return result;
        }

#ifdef DEBUG_MODE
        DEBUG.emplace( "microGradientYieldFunction", yftmp );
#endif

        yieldFunctionValues[ 2 ] = yftmp[ 0 ];
        yieldFunctionValues[ 3 ] = yftmp[ 1 ];
        yieldFunctionValues[ 4 ] = yftmp[ 2 ];

        //Handle the special case when the stresses are very small which is undefined
        if ( tardigradeVectorTools::dot( PK2Stress, PK2Stress ) < 1e-9 ){
            dMacroFdPK2        = variableVector( PK2Stress.size(), 0 );
            dMacroFdElasticRCG = variableVector( elasticRightCauchyGreen.size(), 0 );
        }

        if ( tardigradeVectorTools::dot( referenceMicroStress, referenceMicroStress ) < 1e-9 ){
            dMicroFdSigma      = variableVector( referenceMicroStress.size(), 0 );
            dMicroFdElasticRCG = variableVector( elasticRightCauchyGreen.size(), 0 );
        }

        if ( tardigradeVectorTools::dot( referenceHigherOrderStress, referenceHigherOrderStress ) < 1e-9 ){
            dMicroGradientFdM          = variableMatrix( yftmp.size(), variableVector( referenceHigherOrderStress.size(), 0 ) );
            dMicroGradientFdElasticRCG = variableMatrix( yftmp.size(), variableVector( elasticRightCauchyGreen.size(), 0 ) );
        }

        return NULL;
    }

    errorOut computeCohesion( const variableType &macroStrainISV, const variableType &microStrainISV,
                              const variableVector &microGradientStrainISV,
                              const parameterVector &macroHardeningParameters, const parameterVector &microHardeningParameters,
                              const parameterVector &microGradientHardeningParameters,
                              variableType &macroCohesion, variableType &microCohesion,
                              variableVector &microGradientCohesion ){
        /*!
         * Compute the cohesion value from the strain-like ISVs
         *
         * :param const variableType &macroStrainISV: The macro strain-like ISV
         * :param const variableType &microStrainISV: The micro strain-like ISV
         * :param const variableVector &microGradientStrainISV: The micro gradient strain-like ISV
         * :param const parameterVector &macroHardeningParameters: The hardening parameters for the macro cohesion.
         * :param const parameterVector &microHardeningParameters: The hardening parameters for the micro cohesion.
         * :param const parameterVector &microGradientHardeningParameters: The hardening parameters for the micro
         *     gradient cohesion.
         * :param variableType &macroCohesion: The macro cohesion
         * :param variableType &microCohesion: The micro cohesion
         * :param variableVector &microGradientCohesion: The micro gradient cohesion
         */

        if ( macroHardeningParameters.size() != 2 ){
            return new errorNode( "computeCohesion",
                                  "The macro hardening parameters must have a size of 2" );
        }

        if ( microHardeningParameters.size() != 2 ){
            return new errorNode( "computeCohesion",
                                  "The micro hardening parameters must have a size of 2" );
        }

        if ( microGradientHardeningParameters.size() != 2 ){
            return new errorNode( "computeCohesion",
                                  "The micro gradient hardening parameters must have a size of 2" );
        }

        if ( std::isnan( macroStrainISV ) ){
            return new errorNode( __func__, "The macro-strain ISV is nan" );
        }

        if ( std::isnan( microStrainISV ) ){
            return new errorNode( __func__, "The micro-strain ISV is nan" );
        }

        for ( auto mGISV = microGradientStrainISV.begin( ); mGISV != microGradientStrainISV.end( ); mGISV++ ){
            if ( std::isnan( *mGISV ) ){
                return new errorNode( __func__, "The " + std::to_string( mGISV - microGradientStrainISV.begin( ) ) + "th index of the micro gradient strain ISV is nan" );
            }
        }

        macroCohesion = macroHardeningParameters[ 0 ] + macroHardeningParameters[ 1 ] * macroStrainISV;

        microCohesion = microHardeningParameters[ 0 ] + microHardeningParameters[ 1 ] * microStrainISV;

        microGradientCohesion = microGradientHardeningParameters[ 0 ]
                              + microGradientHardeningParameters[ 1 ] * microGradientStrainISV;

        return NULL;
    }

    errorOut computeCohesion( const variableType &macroStrainISV, const variableType &microStrainISV,
                              const variableVector &microGradientStrainISV,
                              const parameterVector &macroHardeningParameters, const parameterVector &microHardeningParameters,
                              const parameterVector &microGradientHardeningParameters,
                              variableType &macroCohesion, variableType &microCohesion,
                              variableVector &microGradientCohesion,
                              variableType &dMacroCdMacroStrainISV, variableType &dMicroCdMicroStrainISV,
                              variableMatrix &dMicroGradientCdMicroGradientStrainISV ){
        /*!
         * Compute the cohesion value from the strain-like ISVs
         *
         * :param const variableType &macroStrainISV: The macro strain-like ISV
         * :param const variableType &microStrainISV: The micro strain-like ISV
         * :param const variableVector &microGradientStrainISV: The micro gradient strain-like ISV
         * :param const parameterVector &macroHardeningParameters: The hardening parameters for the macro cohesion.
         * :param const parameterVector &microHardeningParameters: The hardening parameters for the micro cohesion.
         * :param const parameterVector &microGradientHardeningParameters: The hardening parameters for the micro
         *     gradient cohesion.
         * :param variableType &macroCohesion: The macro cohesion
         * :param variableType &microCohesion: The micro cohesion
         * :param variableVector &microGradientCohesion: The micro gradient cohesion
         * :param variableType &dMacroCdMacroStrainISV: The Jacobian of the macro cohesion w.r.t. the macro strain ISV
         * :param variableType &dMairoCdMicroStrainISV: The Jacobian of the micro cohesion w.r.t. the micro strain ISV
         * :param variableMatrix &dMicroGradientCdMicroGradientStrainISV: The Jacobian of the micro gradient cohesion
         *      w.r.t. the micro gradient strain ISV
         */

        errorOut error = computeCohesion( macroStrainISV, microStrainISV, microGradientStrainISV,
                                          macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters,
                                          macroCohesion, microCohesion, microGradientCohesion );

        if ( error ){
            errorOut result = new errorNode( "computeCohesion (jacobian)",
                                             "Error in computation of the cohesion" );
            result->addNext( error );
            return result;
        }

        //Assemble the Jacobians
        constantMatrix eye = tardigradeVectorTools::eye< constantType >( microGradientStrainISV.size() );

        dMacroCdMacroStrainISV = macroHardeningParameters[ 1 ];
        dMicroCdMicroStrainISV = microHardeningParameters[ 1 ];
        dMicroGradientCdMicroGradientStrainISV = microGradientHardeningParameters[ 1 ] * eye;

        return NULL;
    }

    errorOut computePlasticMultiplierResidual( const tardigradeSolverTools::floatVector &x, const tardigradeSolverTools::floatMatrix &floatArgs,
                                               const tardigradeSolverTools::intMatrix &intArgs, tardigradeSolverTools::floatVector &residual,
                                               tardigradeSolverTools::floatMatrix &jacobian, tardigradeSolverTools::floatMatrix &floatOuts,
                                               tardigradeSolverTools::intMatrix &intOuts
#ifdef DEBUG_MODE
                                               , tardigradeSolverTools::debugMap &DEBUG
#endif
                                             ){
        /*!
         * Compute the residual on the plastic multiplier residual
         * 
         * :param tardigradeSolverTools::floatVector &x: The unknown vector. Organized as
         *     [ macroGamma, microGamma, microGradientGamma ]
         * :param const tardigradeSolverTools::floatMatrix &floatArgs: The floating point arguments which do not vary
         *     during the solve.
         * :param const tardigradeSolverTools::intMatrix &intArgs: The integer arguments which do not vary during 
         *     the solve.
         * :param tardigradeSolverTools::floatVector &residual: The value of the residual. This will be the 
         *     the values in the x vector - the estimated amount of plastic deformation
         * :param tardigradeSolverTools::floatMatrix &jacobian: The jacobian matrix
         * :param tardigradeSolverTools::floatMatrix &floatOuts: The floating point values that do change during the solve.
         * :param tardigradeSolverTools::intMatrix &intOuts: The integer values that do change during the solve.
         * :param std::map< std::string, tardigradeSolverTools::floatVector > &DEBUG: The debug map. Only available if
         *     DEBUG_MODE is defined.
         *
         * Ordering of floatArgs
         * floatArgs[  0 ] = Dt, The change in time
         * floatArgs[  1 ] = currentDeformationGradient, The current value of the deformation gradient
         * floatArgs[  2 ] = currentMicroDeformation, The current value of the micro-deformation
         * floatArgs[  3 ] = currentGradientMicroDeformation, The current value of the gradient of the
         *     micro deformation in the reference configuration.
         *     multipliers
         * floatArgs[  4 ] = previousMacroGamma, The previous value of the macro plastic multiplier
         * floatArgs[  5 ] = previousMicroGamma, The previous value of the micro plastic multiplier
         * floatArgs[  6 ] = previousMicroGradientGamma, The previous values of the micro plastic
         *     multipliers
         * floatArgs[  7 ] = previousPlasticDeformationGradient, The previous value of the plastic
         *     deformation gradient.
         * floatArgs[  8 ] = previousPlasticMicroDeformation, The previous value of the plastic micro
         *     deforamtion.
         * floatArgs[  9 ] = previousPlasticGradientMicroDeformation, The previous value of the plastic
         *     intermediate configuration gradient of the micro deformation.
         * floatArgs[ 10 ] = previousMacroStrainISV, The previous value of the macro strain-like ISV.
         * floatArgs[ 11 ] = previousMicroStrainISV, The previous value of the micro strain-like ISV.
         * floatArgs[ 12 ] = previousMicroGradientStrainISV, The previous values of the micro gradient
         *     strain-like ISVs.
         * floatArgs[ 13 ] = previousdMacroGdMicroCohesion, The previous value of the Jacobian of the
         *     macro flow direction w.r.t. the macro cohesion value.
         * floatArgs[ 14 ] = previousdMicroGdMicroCohesion, The previous value of the Jacobian of the
         *     micro flow direction w.r.t. the micro cohesion value.
         * floatArgs[ 15 ] = previousdMicroGradientGdMicroGradientCohesion, The previous value of the
         *     Jacobian of the micro gradient flow direction w.r.t. the micro gradient cohesion value.
         * floatArgs[ 16 ] = previousPlasticMacroVelocityGradient, The previous plastic macro velocity
         *     gradient.
         * floatArgs[ 17 ] = previousPlasticMicroVelocityGradient, The previous plastic micro velocity
         *     gradient.
         * floatArgs[ 18 ] = previousPlasticMicroGradientVelocityGradient, The previous plastic micro
         *     gradient velocity gradient.
         * floatArgs[ 19 ] = macroFlowParameters, The macro flow parameters.
         * floatArgs[ 20 ] = microFlowParameters, The micro flow parameters.
         * floatArgs[ 21 ] = microGradientFlowParameters, The micro gradient flow parameters.
         * floatArgs[ 22 ] = macroHardeningParameters, The macro hardening parameters.
         * floatArgs[ 23 ] = microHardeningParameters, The micro hardening parameters.
         * floatArgs[ 24 ] = microGradientHardeningParameters, The micro gradient hardening parameters.
         * floatArgs[ 25 ] = macroYieldParameters, The yield parameters for the macro yield surface
         * floatArgs[ 26 ] = microYieldParameters, The yield parameters for the micro yield surface
         * floatArgs[ 27 ] = microGradientYieldParameters, The yield parameters for the micro gradient yield surface
         * floatArgs[ 28 ] = Amatrix, The A stiffness tensor.
         * floatArgs[ 29 ] = Bmatrix, The B stiffness tensor.
         * floatArgs[ 30 ] = Cmatrix, The C stiffness tensor.
         * floatArgs[ 31 ] = Dmatrix, The D stiffness tensor.
         * floatArgs[ 32 ] = alphaMacro, The macro integration parameter.
         * floatArgs[ 33 ] = alphaMicro, The micro integration parameter.
         * floatArgs[ 34 ] = alphaMicroGradient, The micro gradient integration parameter.
         *
         * Ordering of intArgs
         * intArgs[ 0 ] = evaluateFullDerivatives, Flag which indicates if the full derivatives should be
         *     computed.
         * 
         * Ordering of floatOuts
         * floatOuts[  0 ] = currentPK2Stress, The current value of the second Piola-Kirchoff stress
         * floatOuts[  1 ] = currentReferenceMicroStress, The current value of the reference micro
         *     stress.
         * floatOuts[  2 ] = currentReferenceHigherOrderStress, The current value of the reference
         *     higher order stress.
         * floatOuts[  3 ] = currentMacroStrainISV, The current value of the macro strain-like ISV
         * floatOuts[  4 ] = currentMacroStrainISV, The current value of the micro strain-like ISV
         * floatOuts[  5 ] = currentMacroStrainISV, The current value of the micro gradient strain-like ISV
         * floatOuts[  6 ] = currentPlasticDeformationGradient, The current value of the plastic deformation gradient.
         * floatOuts[  7 ] = currentPlasticMicroDeformation, The current value of the plastic micro deforamtion
         * floatOuts[  8 ] = currentPlasticGradientMicroDeformation, The current value of the plastic gradient of the
         *     micro deformation.
         * floatOuts[  9 ] = dPlasticDeformationdx, The partial derivative of the plastic deformation
         *     w.r.t. the solution vector.
         * floatOuts[ 10 ] = dStressdx, The partial derivative of the stresses w.r.t. the 
         *     solution vector.
         * floatOuts[ 11 ] = dPlasticDeformationdDeformation, The partial derivative of the 
         *     plastic deformation w.r.t. the fundamental deformation measures.
         * floatOuts[ 12 ] = dStressdDeformation, The partial derivative of the stresses w.r.t.
         *     the fundamental deformation measures.
         * floatOuts[ 13 ] = dResidualdDeformation, The partial derivative of the residual w.r.t.
         *     the fundamental deformation measures.
         *
         * Ordering of intOuts
         * intOuts[ 0 ] = flags from the sub-Newton-Raphson process. There are two values:
         *     plasticDeformationConvergenceFlag, The convergence flag from the solution of the plastic deformation.
         *     plasticDeformationFatalErrorFlag, The fatal error flag from the solution of the plastic deformation.
         */

        if ( x.size() != 5 ){
            return new errorNode( "computePlasticMultiplierResidual",
                                  "The x vector must have a length of 5" );
        }

        if ( floatArgs.size() != 35 ){
            return new errorNode( "computePlasticMultiplierResidual",
                                  "The floating point argument matrix floatArgs must have a length of 35" );
        }

        if ( intArgs.size() != 1 ){
            return new errorNode( "computePlasticMultiplierResidual",
                                  "The integer argument matrix intArgs must have a length of 1" );
        }

        if ( intArgs[ 0 ].size() != 1 ){
            return new errorNode( "computePlasticMultiplierResidual",
                                  "The deformation measure evaluation flag must have a size of 1" );
        }

        bool evaluateFullDerivatives = false;
        if ( intArgs[ 0 ][ 0 ] > 0 ){
            evaluateFullDerivatives = true;
        }

        if ( evaluateFullDerivatives ){
            if ( floatOuts.size() != 14 ){
                return new errorNode( "computePlasticMultiplierResidual",
                                      "The floating point output matrix floatOuts must have a length of 14" );
            }
        }
        else{ 
            if ( floatOuts.size() != 9 ){
                return new errorNode( "computePlasticMultiplierResidual",
                                      "The floating point output matrix floatOuts must have a length of 9" );
            }
        }

        if ( intOuts.size() != 1 ){
            return new errorNode( "computePlasticMultiplierResidual",
                                  "The integer output matrix intOuts must have a length of 1" );
        }

        /*=============================
        | Extract the incoming values |
        =============================*/

        const variableVector currentMacroGamma( x.begin(), x.begin() + 1 );
        const variableVector currentMicroGamma( x.begin() + 1, x.begin() + 2 );
        const variableVector currentMicroGradientGamma( x.begin() + 2, x.begin() + 5 );

//        unsigned int ii = 0;
//        const constantType    *Dt                                            = &floatArgs[ ii++ ][ 0 ];
//        const variableVector  *currentDeformationGradient                    = &floatArgs[ ii++ ];
//        const variableVector  *currentMicroDeformation                       = &floatArgs[ ii++ ];
//        const variableVector  *currentGradientMicroDeformation               = &floatArgs[ ii++ ];
//        const variableType    *previousMacroGamma                            = &floatArgs[ ii++ ][ 0 ];
//        const variableType    *previousMicroGamma                            = &floatArgs[ ii++ ][ 0 ];
//        const variableVector  *previousMicroGradientGamma                    = &floatArgs[ ii++ ];
        const variableVector  *previousPlasticDeformationGradient            = &floatArgs[ 7 ];
        const variableVector  *previousPlasticMicroDeformation               = &floatArgs[ 8 ];
        const variableVector  *previousPlasticGradientMicroDeformation       = &floatArgs[ 9 ];
//        const variableType    *previousMacroStrainISV                        = &floatArgs[ ii++ ][ 0 ];
//        const variableType    *previousMicroStrainISV                        = &floatArgs[ ii++ ][ 0 ];
//        const variableVector  *previousMicroGradientStrainISV                = &floatArgs[ ii++ ];
//        const variableType    *previousdMacroGdMacroCohesion                 = &floatArgs[ ii++ ][ 0 ];
//        const variableType    *previousdMicroGdMicroCohesion                 = &floatArgs[ ii++ ][ 0 ];
//        const variableMatrix   previousdMicroGradientGdMicroGradientCohesion = tardigradeVectorTools::inflate( floatArgs[ ii++ ], 3, 3 );
//        const variableVector  *previousPlasticMacroVelocityGradient          = &floatArgs[ ii++ ];
//        const variableVector  *previousPlasticMicroVelocityGradient          = &floatArgs[ ii++ ];
//        const variableVector  *previousPlasticMicroGradientVelocityGradient  = &floatArgs[ ii++ ];
//        const parameterVector *macroFlowParameters                           = &floatArgs[ ii++ ];
//        const parameterVector *microFlowParameters                           = &floatArgs[ ii++ ];
//        const parameterVector *microGradientFlowParameters                   = &floatArgs[ ii++ ];
//        const parameterVector *macroHardeningParameters                      = &floatArgs[ ii++ ];
//        const parameterVector *microHardeningParameters                      = &floatArgs[ ii++ ];
//        const parameterVector *microGradientHardeningParameters              = &floatArgs[ ii++ ];
        const parameterVector *macroYieldParameters                          = &floatArgs[ 25 ];
        const parameterVector *microYieldParameters                          = &floatArgs[ 26 ];
        const parameterVector *microGradientYieldParameters                  = &floatArgs[ 27 ];
//        const parameterVector *Amatrix                                       = &floatArgs[ ii++ ];
//        const parameterVector *Bmatrix                                       = &floatArgs[ ii++ ];
//        const parameterVector *Cmatrix                                       = &floatArgs[ ii++ ];
//        const parameterVector *Dmatrix                                       = &floatArgs[ ii++ ];
//        const parameterType   *alphaMacro                                    = &floatArgs[ ii++ ][ 0 ];
//        const parameterType   *alphaMicro                                    = &floatArgs[ ii++ ][ 0 ];
//        const parameterType   *alphaMicroGradient                            = &floatArgs[ ii++ ][ 0 ];

        //Construct the inputs for the solve for the plastic deformation measure

        tardigradeSolverTools::floatMatrix floatArgsPlasticDeformation =
            {
                floatArgs[  0 ], //Dt
                floatArgs[  1 ], //Current deformation gradient
                floatArgs[  2 ], //Current micro deformation
                floatArgs[  3 ], //Current gradient micro deformation
                { currentMacroGamma },
                { currentMicroGamma },
                currentMicroGradientGamma,
                floatArgs[  4 ], //previous macro gamma
                floatArgs[  5 ], //previous micro gamma
                floatArgs[  6 ], //previous micro gradient gamma
                floatArgs[  7 ], //Previous plastic deformation gradient
                floatArgs[  8 ], //Previous plastic micro deformation
                floatArgs[  9 ], //Previous plastic gradient micro deformation
                floatArgs[ 10 ], //Previous macro strain ISV
                floatArgs[ 11 ], //Previous micro strain ISV
                floatArgs[ 12 ], //Previous micro gradient strain ISV
                floatArgs[ 13 ], //Previous dMacroGdMacroCohesion
                floatArgs[ 14 ], //Previous dMicroGdMicroCohesion
                floatArgs[ 15 ], //Previous dMicroGradientGdMicroGradientCohesion
                floatArgs[ 16 ], //PreviousPlasticMacroVelocityGradient
                floatArgs[ 17 ], //PreviousPlasticMicroVelocityGradient
                floatArgs[ 18 ], //PreviousPlasticMicroGradientVelocityGradient
                floatArgs[ 19 ], //macro flow parameters
                floatArgs[ 20 ], //micro flow parameters
                floatArgs[ 21 ], //micro gradient flow parameters
                floatArgs[ 22 ], //macro hardening parameters
                floatArgs[ 23 ], //micro hardening parameters
                floatArgs[ 24 ], //micro gradient hardening parameters
                floatArgs[ 28 ], //A matrix
                floatArgs[ 29 ], //B matrix
                floatArgs[ 30 ], //C matrix
                floatArgs[ 31 ], //D matrix
                floatArgs[ 32 ], //alpha macro
                floatArgs[ 33 ], //alpha micro
                floatArgs[ 34 ], //alpha micro gradient
            };

        /*==================================
        | Compute the plastic deformations |
        ==================================*/

//        std::cout << "    x: "; tardigradeVectorTools::print( x );

        //Assemble the values required for the non-linear solve

        tardigradeSolverTools::intMatrix intArgsPlasticDeformation = { { 1 } };

        tardigradeSolverTools::intMatrix intOutsPlasticDeformation = { };

        tardigradeSolverTools::floatMatrix floatOutsPlasticDeformation =
            {
                {}, {}, {},
                {}, {}, {},
                {}, {}, {},
                {}, {},
                {}, {}, {}, {},
                {}, {}
            };

        //Wrap the plastic deformation measure function
        tardigradeSolverTools::stdFncNLFJ func
                = static_cast< tardigradeSolverTools::NonLinearFunctionWithJacobian >( computePlasticDeformationResidual );

        tardigradeSolverTools::floatVector plasticDeformationX0
            = tardigradeVectorTools::appendVectors( { *previousPlasticDeformationGradient,
                                            *previousPlasticMicroDeformation,
                                            *previousPlasticGradientMicroDeformation } );

        tardigradeSolverTools::floatVector currentPlasticDeformation( 45, 0 );

        bool convergeFlag, fatalErrorFlag;

        tardigradeSolverTools::solverType plasticDeformationLinearSolver;
        tardigradeSolverTools::floatMatrix plasticDeformationJacobian;

#ifdef DEBUG_MODE
        tardigradeSolverTools::iterationMap plasticDeformationDEBUG;
#endif

        //Solve for the plastic deformation measures.
        errorOut error = tardigradeSolverTools::newtonRaphson( func, plasticDeformationX0, currentPlasticDeformation,
                                                     convergeFlag, fatalErrorFlag,
                                                     floatOutsPlasticDeformation, intOutsPlasticDeformation,
                                                     floatArgsPlasticDeformation, intArgsPlasticDeformation,
                                                     plasticDeformationLinearSolver, plasticDeformationJacobian,
#ifdef DEBUG_MODE
                                                     plasticDeformationDEBUG,
#endif
                                                     20, 1e-9, 1e-9, 1e-4, 5, false );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMultiplierResidual", "Error in solution of plastic deformation" );
            result->addNext( error );
            intOuts[ 0 ] = { ( int )convergeFlag, ( int )fatalErrorFlag };
            if ( ( !fatalErrorFlag ) && ( !convergeFlag ) ){
                //Alert the Newton-Raphson solver that there is a sub-level convergence problem.
                residual.resize( 101 );
                jacobian.resize( 212 );
            }
            return result;
        }

        //Solve for the Jacobian of the plastic deformation w.r.t. the gammas
        
        tardigradeSolverTools::floatVector plasticDeformationJacobianVector = tardigradeVectorTools::appendVectors( plasticDeformationJacobian );
        Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > > 
            dPDResidualdPD( plasticDeformationJacobianVector.data(), 45, 45 );

        Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
            dPDResidualdGammas( floatOutsPlasticDeformation[ 15 ].data(), 45, 5 );

        tardigradeSolverTools::floatVector dCurrentPlasticDeformationdGammas( 45 * 5, 0 );
        Eigen::Map< Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
            dPDdG( dCurrentPlasticDeformationdGammas.data(), 45, 5 );

        tardigradeSolverTools::solverType linearSolver( dPDResidualdPD );

        if ( linearSolver.rank() < 45 ){
            return new errorNode( "computePlasticMultiplierResidual", "The plastic deformation Jacobian is not full rank" );
        }

        dPDdG = -linearSolver.solve( dPDResidualdGammas );

        tardigradeSolverTools::floatVector dCurrentPlasticDeformationdDeformation;
        if ( evaluateFullDerivatives ){
            
            //Map the derivative of the residual w.r.t. the deformation to an Eigen Matrix
            Eigen::Map< Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                dPDResidualdF( floatOutsPlasticDeformation[ 13 ].data(), 45, 45 );

            //Map the derivative of the plastic deformation w.r.t. the total deformation to an Eigen Matrix
            dCurrentPlasticDeformationdDeformation.resize( 45 * 45 );
            Eigen::Map< Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
                dPDdF( dCurrentPlasticDeformationdDeformation.data(), 45, 45 );

            //Solve for the gradient of the plastic deformation w.r.t. the total deformation
            dPDdF = -linearSolver.solve( dPDResidualdF );

#ifdef DEBUG_MODE
            DEBUG.emplace( "dCurrentPlasticDeformationdDeformation", dCurrentPlasticDeformationdDeformation );
#endif
        }

#ifdef DEBUG_MODE
        DEBUG.emplace( "currentPlasticDeformation", currentPlasticDeformation );
        DEBUG.emplace( "dCurrentPlasticDeformationdGammas", dCurrentPlasticDeformationdGammas );
#endif

        //Solve for the Jacobian of the stresses w.r.t. the gammas

        tardigradeSolverTools::floatMatrix dStressdGammas
            = tardigradeVectorTools::inflate( tardigradeVectorTools::matrixMultiply( floatOutsPlasticDeformation[ 11 ],
                                                                 dCurrentPlasticDeformationdGammas, 45, 45, 45, 5 ), 45, 5 );

        if ( evaluateFullDerivatives ){
            floatOuts[  9 ] = dCurrentPlasticDeformationdGammas;
            floatOuts[ 10 ] = tardigradeVectorTools::appendVectors( dStressdGammas );
            floatOuts[ 11 ] = dCurrentPlasticDeformationdDeformation;
            floatOuts[ 12 ] = floatOutsPlasticDeformation[ 12 ]
                            + tardigradeVectorTools::matrixMultiply( floatOutsPlasticDeformation[ 11 ],
                                                           dCurrentPlasticDeformationdDeformation,
                                                           45, 45, 45, 45 );
        }

#ifdef DEBUG_MODE
        DEBUG.emplace( "stresses", tardigradeVectorTools::appendVectors( { floatOutsPlasticDeformation[ 0 ],
                                                                 floatOutsPlasticDeformation[ 1 ],
                                                                 floatOutsPlasticDeformation[ 2 ]  } ) );
        DEBUG.emplace( "dStressdGammas", tardigradeVectorTools::appendVectors( dStressdGammas ) );
#endif

        //Construct the Jacobian of the elastic RCG w.r.t. the gammas
        variableMatrix dPlasticDeformationGradientdGammas =
            {
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 0, dCurrentPlasticDeformationdGammas.begin() + 5 * 1 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 1, dCurrentPlasticDeformationdGammas.begin() + 5 * 2 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 2, dCurrentPlasticDeformationdGammas.begin() + 5 * 3 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 3, dCurrentPlasticDeformationdGammas.begin() + 5 * 4 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 4, dCurrentPlasticDeformationdGammas.begin() + 5 * 5 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 5, dCurrentPlasticDeformationdGammas.begin() + 5 * 6 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 6, dCurrentPlasticDeformationdGammas.begin() + 5 * 7 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 7, dCurrentPlasticDeformationdGammas.begin() + 5 * 8 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 8, dCurrentPlasticDeformationdGammas.begin() + 5 * 9 ),
            };

        variableMatrix dElasticRightCauchyGreendGammas
            = tardigradeVectorTools::dot( tardigradeVectorTools::inflate( floatOutsPlasticDeformation[ 10 ], 9, 9 ), dPlasticDeformationGradientdGammas );

        tardigradeSolverTools::floatVector dElasticRightCauchyGreendDeformation;
        if ( evaluateFullDerivatives ){

            dElasticRightCauchyGreendDeformation = tardigradeSolverTools::floatVector( 9 * 45, 0 );
            for ( unsigned int i = 0; i < 9; i++ ){
                for ( unsigned int j = 0; j < 9; j++ ){
                    dElasticRightCauchyGreendDeformation[ 45 * i + j ] = floatOutsPlasticDeformation[ 14 ][ 9 * i + j ];
                }
            }
           
            tardigradeSolverTools::floatVector temp( dCurrentPlasticDeformationdDeformation.begin(),
                                           dCurrentPlasticDeformationdDeformation.begin() + 9 * 45 ); 

            dElasticRightCauchyGreendDeformation +=
                tardigradeVectorTools::matrixMultiply( floatOutsPlasticDeformation[ 10 ], temp, 9, 9, 9, 45 );

#ifdef DEBUG_MODE
            DEBUG.emplace( "dElasticRightCauchyGreendDeformation", dElasticRightCauchyGreendDeformation );
#endif
        }

#ifdef DEBUG_MODE
        DEBUG.emplace( "currentElasticRightCauchyGreen", floatOutsPlasticDeformation[ 9 ] );
        DEBUG.emplace( "dElasticRightCauchyGreendGammas", tardigradeVectorTools::appendVectors( dElasticRightCauchyGreendGammas ) );
#endif

#ifdef DEBUG_MODE
        tardigradeSolverTools::debugMap yieldFunctionDEBUG;
#endif

        /*=============================
        | Compute the cohesion values |
        =============================*/

        //TODO: Currently the cohesion values are computed in the plastic deformation gamma solve. They don't need to be.

        variableMatrix dCohesiondGammas = tardigradeVectorTools::inflate( floatOutsPlasticDeformation[ 16 ], 5, 5 );

        /*===================================
        | Compute the yield function values |
        ===================================*/

        //Compute the yield function values
        variableVector yieldFunctionValues;
        variableType   dMacroFdMacroC, dMicroFdMicroC;
        variableVector dMacroFdPK2, dMacroFdElasticRCG, dMicroFdSigma, dMicroFdElasticRCG;
        variableMatrix dMicroGradientFdM, dMGFdMGC, dMicroGradientFdElasticRCG;

        error = evaluateYieldFunctions( floatOutsPlasticDeformation[ 0 ], floatOutsPlasticDeformation[ 1 ], //The stresses
                                        floatOutsPlasticDeformation[ 2 ],
                                        floatOutsPlasticDeformation[ 6 ][ 0 ], floatOutsPlasticDeformation[ 7 ][ 0 ], //The cohesions
                                        floatOutsPlasticDeformation[ 8 ],
                                        floatOutsPlasticDeformation[ 9 ], //The elastic right cauchy green deformation tensor
                                       *macroYieldParameters, *microYieldParameters, *microGradientYieldParameters,
                                        yieldFunctionValues,
                                        dMacroFdPK2, dMacroFdMacroC, dMacroFdElasticRCG,
                                        dMicroFdSigma, dMicroFdMicroC, dMicroFdElasticRCG,
                                        dMicroGradientFdM, dMGFdMGC, dMicroGradientFdElasticRCG
#ifdef DEBUG_MODE
                                        , yieldFunctionDEBUG
#endif
                                      );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMultiplierResidual", "Error in the computation of the yield functions" );
            result->addNext( error );
            fatalErrorFlag = true;
            return result;
        }

        //Construct the Jacobians of the yield functions
        variableVector zero9( 9, 0 );
        variableVector zero27( 27, 0 );
        variableMatrix dYieldFunctionValuesdStresses =
            {
                    tardigradeVectorTools::appendVectors( { dMacroFdPK2,         zero9,                 zero27 } ),
                    tardigradeVectorTools::appendVectors( {       zero9, dMicroFdSigma,                 zero27 } ),
                    tardigradeVectorTools::appendVectors( {       zero9,         zero9, dMicroGradientFdM[ 0 ] } ),
                    tardigradeVectorTools::appendVectors( {       zero9,         zero9, dMicroGradientFdM[ 1 ] } ),
                    tardigradeVectorTools::appendVectors( {       zero9,         zero9, dMicroGradientFdM[ 2 ] } )
            };

        variableMatrix dYieldFunctionValuesdCohesion =
            {
                { dMacroFdMacroC,             0.,                 0.,                 0.,                 0. },
                {             0., dMicroFdMicroC,                 0.,                 0.,                 0. },
                {             0.,             0., dMGFdMGC[ 0 ][ 0 ], dMGFdMGC[ 0 ][ 1 ], dMGFdMGC[ 0 ][ 2 ] },
                {             0.,             0., dMGFdMGC[ 1 ][ 0 ], dMGFdMGC[ 1 ][ 1 ], dMGFdMGC[ 1 ][ 2 ] },
                {             0.,             0., dMGFdMGC[ 2 ][ 0 ], dMGFdMGC[ 2 ][ 1 ], dMGFdMGC[ 2 ][ 2 ] }
            };

        variableMatrix dYieldFunctionValuesdElasticRightCauchyGreen =
            {
                dMacroFdElasticRCG,
                dMicroFdElasticRCG,
                dMicroGradientFdElasticRCG[ 0 ],
                dMicroGradientFdElasticRCG[ 1 ],
                dMicroGradientFdElasticRCG[ 2 ]
            };

        variableMatrix dYieldFunctionValuesdGammas =
            tardigradeVectorTools::dot( dYieldFunctionValuesdStresses, dStressdGammas )
          + tardigradeVectorTools::dot( dYieldFunctionValuesdCohesion, dCohesiondGammas )
          + tardigradeVectorTools::dot( dYieldFunctionValuesdElasticRightCauchyGreen, dElasticRightCauchyGreendGammas );

        variableVector dYieldFunctionValuesdDeformation;
        if ( evaluateFullDerivatives ){

            dYieldFunctionValuesdDeformation
                = tardigradeVectorTools::matrixMultiply( tardigradeVectorTools::appendVectors( dYieldFunctionValuesdStresses ),
                                                                           floatOuts[ 12 ], 5, 45, 45, 45 )
                + tardigradeVectorTools::matrixMultiply( tardigradeVectorTools::appendVectors( dYieldFunctionValuesdElasticRightCauchyGreen ),
                                                                           dElasticRightCauchyGreendDeformation, 5, 9, 9, 45 );
        }

#ifdef DEBUG_MODE
        DEBUG.emplace( "yieldFunctionValues", yieldFunctionValues );
        DEBUG.emplace( "dYieldFunctionValuesdGammas", tardigradeVectorTools::appendVectors( dYieldFunctionValuesdGammas ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dYieldFunctionValuesdDeformation", dYieldFunctionValuesdDeformation );
        }
#endif

        //Construct the residual and the Jacobian
        residual = tardigradeSolverTools::floatVector( 5, 0 );
        jacobian = tardigradeSolverTools::floatMatrix( 5, tardigradeSolverTools::floatVector( 5, 0 ) );

        variableType dMacYieldFunctionValuedYieldFunctionValue_positive,
                     dMacYieldFunctionValuedYieldFunctionValue_negative;
        tardigradeSolverTools::floatVector rowEye;

        if ( evaluateFullDerivatives ){
            floatOuts[ 13 ].resize( 5 * 45 );
        }

        for ( unsigned int i = 0; i < 5; i++ ){
            residual[ i ] = tardigradeConstitutiveTools::mac( yieldFunctionValues[ i ], dMacYieldFunctionValuedYieldFunctionValue_positive )
                          - x[ i ] * tardigradeConstitutiveTools::mac( -yieldFunctionValues[ i ] );

            if ( tardigradeVectorTools::fuzzyEquals( dMacYieldFunctionValuedYieldFunctionValue_positive,  0. ) ){
                dMacYieldFunctionValuedYieldFunctionValue_negative = 1;
            }
            else{
                dMacYieldFunctionValuedYieldFunctionValue_negative = 0;
            }

            rowEye = tardigradeSolverTools::floatVector( 5, 0 );
            rowEye[ i ] = 1.;

            jacobian[ i ] = ( dMacYieldFunctionValuedYieldFunctionValue_positive + x[ i ] * dMacYieldFunctionValuedYieldFunctionValue_negative ) * dYieldFunctionValuesdGammas[ i ]
                          - tardigradeConstitutiveTools::mac( -yieldFunctionValues[ i ] ) * rowEye;

            if ( evaluateFullDerivatives ){
                for ( unsigned int j = 0; j < 45; j++ ){
                    floatOuts[ 13 ][ 45 * i + j ] = ( dMacYieldFunctionValuedYieldFunctionValue_positive + x[ i ] * dMacYieldFunctionValuedYieldFunctionValue_negative )
                                                  * dYieldFunctionValuesdDeformation[ 45 * i + j ];
                }
            }
        }

        //Update the output vector
        floatOuts[ 0 ] = floatOutsPlasticDeformation[ 0 ];
        floatOuts[ 1 ] = floatOutsPlasticDeformation[ 1 ];
        floatOuts[ 2 ] = floatOutsPlasticDeformation[ 2 ];
        floatOuts[ 3 ] = floatOutsPlasticDeformation[ 3 ];
        floatOuts[ 4 ] = floatOutsPlasticDeformation[ 4 ];
        floatOuts[ 5 ] = floatOutsPlasticDeformation[ 5 ];
        floatOuts[ 6 ] = tardigradeSolverTools::floatVector( currentPlasticDeformation.begin() +  0, currentPlasticDeformation.begin() + 9 );
        floatOuts[ 7 ] = tardigradeSolverTools::floatVector( currentPlasticDeformation.begin() +  9, currentPlasticDeformation.begin() + 18 );
        floatOuts[ 8 ] = tardigradeSolverTools::floatVector( currentPlasticDeformation.begin() + 18, currentPlasticDeformation.begin() + 45 );

        return NULL;
    }

    errorOut computePlasticMultiplierLagrangian( const tardigradeSolverTools::floatVector &x, const tardigradeSolverTools::floatMatrix &floatArgs,
                                                 const tardigradeSolverTools::intMatrix &intArgs, tardigradeSolverTools::floatType &lagrangian,
                                                 tardigradeSolverTools::floatVector &jacobian, tardigradeSolverTools::floatMatrix &floatOuts,
                                                 tardigradeSolverTools::intMatrix &intOuts
#ifdef DEBUG_MODE
                                                 , tardigradeSolverTools::debugMap &DEBUG
#endif
                                             ){
        /*!
         * Compute the lagrangian of the plastic multiplier
         * 
         * :param tardigradeSolverTools::floatVector &x: The unknown vector. Organized as
         *     [ macroGamma, microGamma, microGradientGamma, lambda1, lambda2, lambda3, lambda4, lambda5 ]
         * :param const tardigradeSolverTools::floatMatrix &floatArgs: The floating point arguments which do not vary
         *     during the solve.
         * :param const tardigradeSolverTools::intMatrix &intArgs: The integer arguments which do not vary during 
         *     the solve.
         * :param tardigradeSolverTools::floatType &lagrangian: The value of the lagrangian. 
         * :param tardigradeSolverTools::floatVector &jacobian: The jacobian matrix
         * :param tardigradeSolverTools::floatMatrix &floatOuts: The floating point values that do change during the solve.
         * :param tardigradeSolverTools::intMatrix &intOuts: The integer values that do change during the solve.
         * :param std::map< std::string, tardigradeSolverTools::floatVector > &DEBUG: The debug map. Only available if
         *     DEBUG_MODE is defined.
         *
         * Ordering of floatArgs
         * floatArgs[  0 ] = Dt, The change in time
         * floatArgs[  1 ] = currentDeformationGradient, The current value of the deformation gradient
         * floatArgs[  2 ] = currentMicroDeformation, The current value of the micro-deformation
         * floatArgs[  3 ] = currentGradientMicroDeformation, The current value of the gradient of the
         *     micro deformation in the reference configuration.
         *     multipliers
         * floatArgs[  4 ] = previousMacroGamma, The previous value of the macro plastic multiplier
         * floatArgs[  5 ] = previousMicroGamma, The previous value of the micro plastic multiplier
         * floatArgs[  6 ] = previousMicroGradientGamma, The previous values of the micro plastic
         *     multipliers
         * floatArgs[  7 ] = previousPlasticDeformationGradient, The previous value of the plastic
         *     deformation gradient.
         * floatArgs[  8 ] = previousPlasticMicroDeformation, The previous value of the plastic micro
         *     deforamtion.
         * floatArgs[  9 ] = previousPlasticGradientMicroDeformation, The previous value of the plastic
         *     intermediate configuration gradient of the micro deformation.
         * floatArgs[ 10 ] = previousMacroStrainISV, The previous value of the macro strain-like ISV.
         * floatArgs[ 11 ] = previousMicroStrainISV, The previous value of the micro strain-like ISV.
         * floatArgs[ 12 ] = previousMicroGradientStrainISV, The previous values of the micro gradient
         *     strain-like ISVs.
         * floatArgs[ 13 ] = previousdMacroGdMicroCohesion, The previous value of the Jacobian of the
         *     macro flow direction w.r.t. the macro cohesion value.
         * floatArgs[ 14 ] = previousdMicroGdMicroCohesion, The previous value of the Jacobian of the
         *     micro flow direction w.r.t. the micro cohesion value.
         * floatArgs[ 15 ] = previousdMicroGradientGdMicroGradientCohesion, The previous value of the
         *     Jacobian of the micro gradient flow direction w.r.t. the micro gradient cohesion value.
         * floatArgs[ 16 ] = previousPlasticMacroVelocityGradient, The previous plastic macro velocity
         *     gradient.
         * floatArgs[ 17 ] = previousPlasticMicroVelocityGradient, The previous plastic micro velocity
         *     gradient.
         * floatArgs[ 18 ] = previousPlasticMicroGradientVelocityGradient, The previous plastic micro
         *     gradient velocity gradient.
         * floatArgs[ 19 ] = macroFlowParameters, The macro flow parameters.
         * floatArgs[ 20 ] = microFlowParameters, The micro flow parameters.
         * floatArgs[ 21 ] = microGradientFlowParameters, The micro gradient flow parameters.
         * floatArgs[ 22 ] = macroHardeningParameters, The macro hardening parameters.
         * floatArgs[ 23 ] = microHardeningParameters, The micro hardening parameters.
         * floatArgs[ 24 ] = microGradientHardeningParameters, The micro gradient hardening parameters.
         * floatArgs[ 25 ] = macroYieldParameters, The yield parameters for the macro yield surface
         * floatArgs[ 26 ] = microYieldParameters, The yield parameters for the micro yield surface
         * floatArgs[ 27 ] = microGradientYieldParameters, The yield parameters for the micro gradient yield surface
         * floatArgs[ 28 ] = Amatrix, The A stiffness tensor.
         * floatArgs[ 29 ] = Bmatrix, The B stiffness tensor.
         * floatArgs[ 30 ] = Cmatrix, The C stiffness tensor.
         * floatArgs[ 31 ] = Dmatrix, The D stiffness tensor.
         * floatArgs[ 32 ] = alphaMacro, The macro integration parameter.
         * floatArgs[ 33 ] = alphaMicro, The micro integration parameter.
         * floatArgs[ 34 ] = alphaMicroGradient, The micro gradient integration parameter.
         *
         * Ordering of intArgs
         * intArgs[ 0 ] = evaluateFullDerivatives, Flag which indicates if the full derivatives should be
         *     computed.
         * 
         * Ordering of floatOuts
         * floatOuts[  0 ] = currentPK2Stress, The current value of the second Piola-Kirchoff stress
         * floatOuts[  1 ] = currentReferenceMicroStress, The current value of the reference micro
         *     stress.
         * floatOuts[  2 ] = currentReferenceHigherOrderStress, The current value of the reference
         *     higher order stress.
         * floatOuts[  3 ] = currentMacroStrainISV, The current value of the macro strain-like ISV
         * floatOuts[  4 ] = currentMacroStrainISV, The current value of the micro strain-like ISV
         * floatOuts[  5 ] = currentMacroStrainISV, The current value of the micro gradient strain-like ISV
         * floatOuts[  6 ] = currentPlasticDeformationGradient, The current value of the plastic deformation gradient.
         * floatOuts[  7 ] = currentPlasticMicroDeformation, The current value of the plastic micro deforamtion
         * floatOuts[  8 ] = currentPlasticGradientMicroDeformation, The current value of the plastic gradient of the
         *     micro deformation.
         * floatOuts[  9 ] = dPlasticDeformationdx, The partial derivative of the plastic deformation
         *     w.r.t. the solution vector.
         * floatOuts[ 10 ] = dStressdx, The partial derivative of the stresses w.r.t. the 
         *     solution vector.
         * floatOuts[ 11 ] = dPlasticDeformationdDeformation, The partial derivative of the 
         *     plastic deformation w.r.t. the fundamental deformation measures.
         * floatOuts[ 12 ] = dStressdDeformation, The partial derivative of the stresses w.r.t.
         *     the fundamental deformation measures.
         * floatOuts[ 13 ] = dResidualdDeformation, The partial derivative of the residual w.r.t.
         *     the fundamental deformation measures.
         *
         * Ordering of intOuts
         * intOuts[ 0 ] = flags from the sub-Newton-Raphson process. There are two values:
         *     plasticDeformationConvergenceFlag, The convergence flag from the solution of the plastic deformation.
         *     plasticDeformationFatalErrorFlag, The fatal error flag from the solution of the plastic deformation.
         */

        if ( x.size() != 10 ){
            return new errorNode( "computePlasticMultiplierLagrangian",
                                  "The x vector must have a length of 10" );
        }

        if ( floatArgs.size() != 35 ){
            return new errorNode( "computePlasticMultiplierLagrangian",
                                  "The floating point argument matrix floatArgs must have a length of 35" );
        }

        if ( intArgs.size() != 0 ){
            return new errorNode( "computePlasticMultiplierLagrangian",
                                  "The integer argument matrix intArgs must have a length of 0" );
        }

        if ( floatOuts.size() != 0 ){
            return new errorNode( "computePlasticMultiplierLagrangian",
                                  "The floating point output matrix floatOuts must have a length of 0" );
        }

        if ( intOuts.size() != 1 ){
            return new errorNode( "computePlasticMultiplierLagrangian",
                                  "The integer output matrix intOuts must have a length of 1" );
        }

        /*=============================
        | Extract the incoming values |
        =============================*/

        const variableVector currentMacroGamma( x.begin(), x.begin() + 1 );
        const variableVector currentMicroGamma( x.begin() + 1, x.begin() + 2 );
        const variableVector currentMicroGradientGamma( x.begin() + 2, x.begin() + 5 );
        const variableVector macroLambda( x.begin() + 5, x.begin() + 6 );
        const variableVector microLambda( x.begin() + 6, x.begin() + 7 );
        const variableVector microGradientLambda( x.begin() + 7, x.begin() + 10 );

//        unsigned int ii = 0;
//        const constantType    *Dt                                            = &floatArgs[ ii++ ][ 0 ];
//        const variableVector  *currentDeformationGradient                    = &floatArgs[ ii++ ];
//        const variableVector  *currentMicroDeformation                       = &floatArgs[ ii++ ];
//        const variableVector  *currentGradientMicroDeformation               = &floatArgs[ ii++ ];
//        const variableType    *previousMacroGamma                            = &floatArgs[ ii++ ][ 0 ];
//        const variableType    *previousMicroGamma                            = &floatArgs[ ii++ ][ 0 ];
//        const variableVector  *previousMicroGradientGamma                    = &floatArgs[ ii++ ];
        const variableVector  *previousPlasticDeformationGradient            = &floatArgs[ 7 ];
        const variableVector  *previousPlasticMicroDeformation               = &floatArgs[ 8 ];
        const variableVector  *previousPlasticGradientMicroDeformation       = &floatArgs[ 9 ];
//        const variableType    *previousMacroStrainISV                        = &floatArgs[ ii++ ][ 0 ];
//        const variableType    *previousMicroStrainISV                        = &floatArgs[ ii++ ][ 0 ];
//        const variableVector  *previousMicroGradientStrainISV                = &floatArgs[ ii++ ];
//        const variableType    *previousdMacroGdMacroCohesion                 = &floatArgs[ ii++ ][ 0 ];
//        const variableType    *previousdMicroGdMicroCohesion                 = &floatArgs[ ii++ ][ 0 ];
//        const variableMatrix   previousdMicroGradientGdMicroGradientCohesion = tardigradeVectorTools::inflate( floatArgs[ ii++ ], 3, 3 );
//        const variableVector  *previousPlasticMacroVelocityGradient          = &floatArgs[ ii++ ];
//        const variableVector  *previousPlasticMicroVelocityGradient          = &floatArgs[ ii++ ];
//        const variableVector  *previousPlasticMicroGradientVelocityGradient  = &floatArgs[ ii++ ];
//        const parameterVector *macroFlowParameters                           = &floatArgs[ ii++ ];
//        const parameterVector *microFlowParameters                           = &floatArgs[ ii++ ];
//        const parameterVector *microGradientFlowParameters                   = &floatArgs[ ii++ ];
//        const parameterVector *macroHardeningParameters                      = &floatArgs[ ii++ ];
//        const parameterVector *microHardeningParameters                      = &floatArgs[ ii++ ];
//        const parameterVector *microGradientHardeningParameters              = &floatArgs[ ii++ ];
        const parameterVector *macroYieldParameters                          = &floatArgs[ 25 ];
        const parameterVector *microYieldParameters                          = &floatArgs[ 26 ];
        const parameterVector *microGradientYieldParameters                  = &floatArgs[ 27 ];
//        const parameterVector *Amatrix                                       = &floatArgs[ ii++ ];
//        const parameterVector *Bmatrix                                       = &floatArgs[ ii++ ];
//        const parameterVector *Cmatrix                                       = &floatArgs[ ii++ ];
//        const parameterVector *Dmatrix                                       = &floatArgs[ ii++ ];
//        const parameterType   *alphaMacro                                    = &floatArgs[ ii++ ][ 0 ];
//        const parameterType   *alphaMicro                                    = &floatArgs[ ii++ ][ 0 ];
//        const parameterType   *alphaMicroGradient                            = &floatArgs[ ii++ ][ 0 ];

        //Construct the inputs for the solve for the plastic deformation measure

        tardigradeSolverTools::floatMatrix floatArgsPlasticDeformation =
            {
                floatArgs[  0 ], //Dt
                floatArgs[  1 ], //Current deformation gradient
                floatArgs[  2 ], //Current micro deformation
                floatArgs[  3 ], //Current gradient micro deformation
                { currentMacroGamma },
                { currentMicroGamma },
                currentMicroGradientGamma,
                floatArgs[  4 ], //previous macro gamma
                floatArgs[  5 ], //previous micro gamma
                floatArgs[  6 ], //previous micro gradient gamma
                floatArgs[  7 ], //Previous plastic deformation gradient
                floatArgs[  8 ], //Previous plastic micro deformation
                floatArgs[  9 ], //Previous plastic gradient micro deformation
                floatArgs[ 10 ], //Previous macro strain ISV
                floatArgs[ 11 ], //Previous micro strain ISV
                floatArgs[ 12 ], //Previous micro gradient strain ISV
                floatArgs[ 13 ], //Previous dMacroGdMacroCohesion
                floatArgs[ 14 ], //Previous dMicroGdMicroCohesion
                floatArgs[ 15 ], //Previous dMicroGradientGdMicroGradientCohesion
                floatArgs[ 16 ], //PreviousPlasticMacroVelocityGradient
                floatArgs[ 17 ], //PreviousPlasticMicroVelocityGradient
                floatArgs[ 18 ], //PreviousPlasticMicroGradientVelocityGradient
                floatArgs[ 19 ], //macro flow parameters
                floatArgs[ 20 ], //micro flow parameters
                floatArgs[ 21 ], //micro gradient flow parameters
                floatArgs[ 22 ], //macro hardening parameters
                floatArgs[ 23 ], //micro hardening parameters
                floatArgs[ 24 ], //micro gradient hardening parameters
                floatArgs[ 28 ], //A matrix
                floatArgs[ 29 ], //B matrix
                floatArgs[ 30 ], //C matrix
                floatArgs[ 31 ], //D matrix
                floatArgs[ 32 ], //alpha macro
                floatArgs[ 33 ], //alpha micro
                floatArgs[ 34 ], //alpha micro gradient
            };

        /*==================================
        | Compute the plastic deformations |
        ==================================*/

        //Assemble the values required for the non-linear solve

        tardigradeSolverTools::intMatrix intArgsPlasticDeformation = { { 1 } };

        tardigradeSolverTools::intMatrix intOutsPlasticDeformation = { };

        tardigradeSolverTools::floatMatrix floatOutsPlasticDeformation =
            {
                {}, {}, {},
                {}, {}, {},
                {}, {}, {},
                {}, {},
                {}, {}, {}, {},
                {}, {}
            };

        //Wrap the plastic deformation measure function
        tardigradeSolverTools::stdFncNLFJ func
                = static_cast< tardigradeSolverTools::NonLinearFunctionWithJacobian >( computePlasticDeformationResidual );

        tardigradeSolverTools::floatVector plasticDeformationX0
            = tardigradeVectorTools::appendVectors( { *previousPlasticDeformationGradient,
                                            *previousPlasticMicroDeformation,
                                            *previousPlasticGradientMicroDeformation } );

        tardigradeSolverTools::floatVector currentPlasticDeformation( 45, 0 );

        bool convergeFlag, fatalErrorFlag;

        tardigradeSolverTools::solverType plasticDeformationLinearSolver;
        tardigradeSolverTools::floatMatrix plasticDeformationJacobian;

#ifdef DEBUG_MODE
        tardigradeSolverTools::iterationMap plasticDeformationDEBUG;
#endif

        //Solve for the plastic deformation measures.
        errorOut error = tardigradeSolverTools::newtonRaphson( func, plasticDeformationX0, currentPlasticDeformation,
                                                     convergeFlag, fatalErrorFlag,
                                                     floatOutsPlasticDeformation, intOutsPlasticDeformation,
                                                     floatArgsPlasticDeformation, intArgsPlasticDeformation,
                                                     plasticDeformationLinearSolver, plasticDeformationJacobian,
#ifdef DEBUG_MODE
                                                     plasticDeformationDEBUG,
#endif
                                                     20, 1e-9, 1e-9, 1e-4, 5, false );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMultiplierLagrangian", "Error in solution of plastic deformation" );
            result->addNext( error );
            intOuts[ 0 ] = { ( int )convergeFlag, ( int )fatalErrorFlag };
            if ( ( !fatalErrorFlag ) && ( !convergeFlag ) ){
                //Alert the Newton-Raphson solver that there is a sub-level convergence problem.
                lagrangian = -101;
                jacobian.resize( 212 );
            }
            return result;
        }

//        std::cout << "    plasticDeformationX0: "; tardigradeVectorTools::print( plasticDeformationX0 );
//        std::cout << "    currentPlasticDeformation: "; tardigradeVectorTools::print( currentPlasticDeformation );

        //Solve for the Jacobian of the plastic deformation w.r.t. the gammas
        
        tardigradeSolverTools::floatVector plasticDeformationJacobianVector = tardigradeVectorTools::appendVectors( plasticDeformationJacobian );
        Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > > 
            dPDResidualdPD( plasticDeformationJacobianVector.data(), 45, 45 );

        Eigen::Map< const Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
            dPDResidualdGammas( floatOutsPlasticDeformation[ 15 ].data(), 45, 5 );

        tardigradeSolverTools::floatVector dCurrentPlasticDeformationdGammas( 45 * 5, 0 );
        Eigen::Map< Eigen::Matrix< tardigradeSolverTools::floatType, -1, -1, Eigen::RowMajor > >
            dPDdG( dCurrentPlasticDeformationdGammas.data(), 45, 5 );

        tardigradeSolverTools::solverType linearSolver( dPDResidualdPD );

        if ( linearSolver.rank() < 45 ){
            return new errorNode( "computePlasticMultiplierLagrangian", "The plastic deformation Jacobian is not full rank" );
        }

        dPDdG = -linearSolver.solve( dPDResidualdGammas );

#ifdef DEBUG_MODE
        DEBUG.emplace( "currentPlasticDeformation", currentPlasticDeformation );
        DEBUG.emplace( "dCurrentPlasticDeformationdGammas", dCurrentPlasticDeformationdGammas );
#endif

        //Solve for the Jacobian of the stresses w.r.t. the gammas

        tardigradeSolverTools::floatMatrix dStressdGammas
            = tardigradeVectorTools::inflate( tardigradeVectorTools::matrixMultiply( floatOutsPlasticDeformation[ 11 ],
                                                                 dCurrentPlasticDeformationdGammas, 45, 45, 45, 5 ), 45, 5 );

#ifdef DEBUG_MODE
        DEBUG.emplace( "stresses", tardigradeVectorTools::appendVectors( { floatOutsPlasticDeformation[ 0 ],
                                                                 floatOutsPlasticDeformation[ 1 ],
                                                                 floatOutsPlasticDeformation[ 2 ]  } ) );
        DEBUG.emplace( "dStressdGammas", tardigradeVectorTools::appendVectors( dStressdGammas ) );
#endif

        //Construct the Jacobian of the elastic RCG w.r.t. the gammas
        variableMatrix dPlasticDeformationGradientdGammas =
            {
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 0, dCurrentPlasticDeformationdGammas.begin() + 5 * 1 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 1, dCurrentPlasticDeformationdGammas.begin() + 5 * 2 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 2, dCurrentPlasticDeformationdGammas.begin() + 5 * 3 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 3, dCurrentPlasticDeformationdGammas.begin() + 5 * 4 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 4, dCurrentPlasticDeformationdGammas.begin() + 5 * 5 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 5, dCurrentPlasticDeformationdGammas.begin() + 5 * 6 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 6, dCurrentPlasticDeformationdGammas.begin() + 5 * 7 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 7, dCurrentPlasticDeformationdGammas.begin() + 5 * 8 ),
                variableVector( dCurrentPlasticDeformationdGammas.begin() + 5 * 8, dCurrentPlasticDeformationdGammas.begin() + 5 * 9 ),
            };

        variableMatrix dElasticRightCauchyGreendGammas
            = tardigradeVectorTools::dot( tardigradeVectorTools::inflate( floatOutsPlasticDeformation[ 10 ], 9, 9 ), dPlasticDeformationGradientdGammas );

#ifdef DEBUG_MODE
        DEBUG.emplace( "currentElasticRightCauchyGreen", floatOutsPlasticDeformation[ 9 ] );
        DEBUG.emplace( "dElasticRightCauchyGreendGammas", tardigradeVectorTools::appendVectors( dElasticRightCauchyGreendGammas ) );
#endif

#ifdef DEBUG_MODE
        tardigradeSolverTools::debugMap yieldFunctionDEBUG;
#endif

        /*=============================
        | Compute the cohesion values |
        =============================*/

        //TODO: Currently the cohesion values are computed in the plastic deformation gamma solve. They don't need to be.

        variableMatrix dCohesiondGammas = tardigradeVectorTools::inflate( floatOutsPlasticDeformation[ 16 ], 5, 5 );

        /*===================================
        | Compute the yield function values |
        ===================================*/

        //Compute the yield function values
        variableVector yieldFunctionValues;
        variableType   dMacroFdMacroC, dMicroFdMicroC;
        variableVector dMacroFdPK2, dMacroFdElasticRCG, dMicroFdSigma, dMicroFdElasticRCG;
        variableMatrix dMicroGradientFdM, dMGFdMGC, dMicroGradientFdElasticRCG;

        error = evaluateYieldFunctions( floatOutsPlasticDeformation[ 0 ], floatOutsPlasticDeformation[ 1 ], //The stresses
                                        floatOutsPlasticDeformation[ 2 ],
                                        floatOutsPlasticDeformation[ 6 ][ 0 ], floatOutsPlasticDeformation[ 7 ][ 0 ], //The cohesions
                                        floatOutsPlasticDeformation[ 8 ],
                                        floatOutsPlasticDeformation[ 9 ], //The elastic right cauchy green deformation tensor
                                       *macroYieldParameters, *microYieldParameters, *microGradientYieldParameters,
                                        yieldFunctionValues,
                                        dMacroFdPK2, dMacroFdMacroC, dMacroFdElasticRCG,
                                        dMicroFdSigma, dMicroFdMicroC, dMicroFdElasticRCG,
                                        dMicroGradientFdM, dMGFdMGC, dMicroGradientFdElasticRCG
#ifdef DEBUG_MODE
                                        , yieldFunctionDEBUG
#endif
                                      );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMultiplierLagrangian", "Error in the computation of the yield functions" );
            result->addNext( error );
            fatalErrorFlag = true;
            return result;
        }

//        std::cout << "    yieldFunctionValues: "; tardigradeVectorTools::print( yieldFunctionValues );

        //Construct the Jacobians of the yield functions
        variableVector zero9( 9, 0 );
        variableVector zero27( 27, 0 );
        variableMatrix dYieldFunctionValuesdStresses =
            {
                    tardigradeVectorTools::appendVectors( { dMacroFdPK2,         zero9,                 zero27 } ),
                    tardigradeVectorTools::appendVectors( {       zero9, dMicroFdSigma,                 zero27 } ),
                    tardigradeVectorTools::appendVectors( {       zero9,         zero9, dMicroGradientFdM[ 0 ] } ),
                    tardigradeVectorTools::appendVectors( {       zero9,         zero9, dMicroGradientFdM[ 1 ] } ),
                    tardigradeVectorTools::appendVectors( {       zero9,         zero9, dMicroGradientFdM[ 2 ] } )
            };

        variableMatrix dYieldFunctionValuesdCohesion =
            {
                { dMacroFdMacroC,             0.,                 0.,                 0.,                 0. },
                {             0., dMicroFdMicroC,                 0.,                 0.,                 0. },
                {             0.,             0., dMGFdMGC[ 0 ][ 0 ], dMGFdMGC[ 0 ][ 1 ], dMGFdMGC[ 0 ][ 2 ] },
                {             0.,             0., dMGFdMGC[ 1 ][ 0 ], dMGFdMGC[ 1 ][ 1 ], dMGFdMGC[ 1 ][ 2 ] },
                {             0.,             0., dMGFdMGC[ 2 ][ 0 ], dMGFdMGC[ 2 ][ 1 ], dMGFdMGC[ 2 ][ 2 ] }
            };

        variableMatrix dYieldFunctionValuesdElasticRightCauchyGreen =
            {
                dMacroFdElasticRCG,
                dMicroFdElasticRCG,
                dMicroGradientFdElasticRCG[ 0 ],
                dMicroGradientFdElasticRCG[ 1 ],
                dMicroGradientFdElasticRCG[ 2 ]
            };

        variableMatrix dYieldFunctionValuesdGammas =
            tardigradeVectorTools::dot( dYieldFunctionValuesdStresses, dStressdGammas )
          + tardigradeVectorTools::dot( dYieldFunctionValuesdCohesion, dCohesiondGammas )
          + tardigradeVectorTools::dot( dYieldFunctionValuesdElasticRightCauchyGreen, dElasticRightCauchyGreendGammas );

#ifdef DEBUG_MODE
        DEBUG.emplace( "yieldFunctionValues", yieldFunctionValues );
        DEBUG.emplace( "dYieldFunctionValuesdGammas", tardigradeVectorTools::appendVectors( dYieldFunctionValuesdGammas ) );
#endif

        //Construct the lagrangian and the Jacobian
        lagrangian = 0;
        jacobian   = tardigradeSolverTools::floatVector( 10, 0 );

        variableType macF;
        variableType dMacYieldFunctionValuedYieldFunctionValue;
        tardigradeSolverTools::floatVector rowEye;
        tardigradeSolverTools::floatVector term1( 5 ), term2( 5 );

        for ( unsigned int i = 0; i < 5; i++ ){

            macF = tardigradeConstitutiveTools::mac( yieldFunctionValues[ i ], dMacYieldFunctionValuedYieldFunctionValue );

            lagrangian += 0.5 * macF * macF + x[ i + 5 ] * x[ i ] * yieldFunctionValues[ i ];

            rowEye = tardigradeSolverTools::floatVector( 5, 0 );
            rowEye[ i ] = 1.;

            term1 = dMacYieldFunctionValuedYieldFunctionValue * dYieldFunctionValuesdGammas[ i ] * macF
                  + x[ i + 5 ] * ( yieldFunctionValues[ i ] * rowEye + x[ i ] * dYieldFunctionValuesdGammas[ i ] );
            term2 = tardigradeSolverTools::floatVector( 5, 0 );
            term2[ i ] = x[ i ] * yieldFunctionValues[ i ];

            jacobian += tardigradeVectorTools::appendVectors( { term1, term2 } );

        }

        return NULL;
    }
}
