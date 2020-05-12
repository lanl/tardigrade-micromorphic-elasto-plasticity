/*!
 * micromorphic_elasto_plasticity.cpp
 *
 * An implementation of a elasto-plastic micromorphic constitutive model 
 * following the derivations of Farhad Shahabi in his dissertation.
 */

#include<micromorphic_elasto_plasticity.h>

namespace micromorphicElastoPlasticity{

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

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableType pressure;
        variableVector deviatoricReferenceStress;
        
        errorOut error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation",
                                             "Error in computation of second-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableType normDevStress = vectorTools::l2norm( deviatoricReferenceStress );

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

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableType pressure;
        variableVector deviatoricReferenceStress;

        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableVector dPressuredStress, dPressuredRCG;
        
        errorOut error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of second-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableType normDevStress = vectorTools::l2norm( deviatoricReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Evaluate the jacobians
        variableVector devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );

        dFdStress = vectorTools::Tdot( dDevStressdStress, devStressDirection )
                  + BAngle * dPressuredStress;

        dFdc = - AAngle;

        dFdElasticRCG = vectorTools::Tdot( dDevStressdRCG, devStressDirection )
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

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableType pressure;
        variableVector deviatoricReferenceStress;
        
        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableVector dPressuredStress, dPressuredRCG;

        variableMatrix d2DevStressdStressdRCG, d2PressuredStressdRCG;

        errorOut error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation (second order jacobian)",
                                             "Error in computation of second-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableType normDevStress = vectorTools::l2norm( deviatoricReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Evaluate the jacobians
        variableVector devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );

        dFdStress = vectorTools::Tdot( dDevStressdStress, devStressDirection )
                  + BAngle * dPressuredStress;

        dFdc = - AAngle;

        dFdElasticRCG = vectorTools::Tdot( dDevStressdRCG, devStressDirection )
                      + BAngle * dPressuredRCG;

        //Evaluate the second-order jacobians
        constantMatrix EYE = vectorTools::eye< constantType >( dim * dim );
        variableMatrix dDevStressDirectiondDevStress = ( EYE - vectorTools::dyadic( devStressDirection, devStressDirection ) ) / ( normDevStress + tol );

        d2FdStress2 = vectorTools::Tdot( dDevStressdStress, vectorTools::dot( dDevStressDirectiondDevStress, dDevStressdStress ) );

        d2FdStressdElasticRCG = vectorTools::Tdot( vectorTools::dot( dDevStressDirectiondDevStress, dDevStressdStress ), dDevStressdRCG )
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

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableVector pressure;
        variableVector deviatoricReferenceStress;
        
        errorOut error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( referenceHigherOrderStress,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation",
                                             "Error in computation of higher-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableVector normDevStress;
        error = micromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress );

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

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableVector pressure;
        variableVector deviatoricReferenceStress;

        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableMatrix dPressuredStress, dPressuredRCG;
        
        errorOut error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( referenceHigherOrderStress,
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
        error = micromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress, dNormDevStressdDevStress );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of the deviatoric higher-order stress norm" );
            result->addNext( error );
            return result;
        }

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Construct the Jacobians
        dFdStress = vectorTools::dot( dNormDevStressdDevStress, dDevStressdStress )
                  + BAngle * dPressuredStress;

        dFdc = -AAngle * vectorTools::eye< constantType >( cohesion.size() );

        dFdElasticRCG = vectorTools::dot( dNormDevStressdDevStress, dDevStressdRCG )
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

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableVector pressure;
        variableVector deviatoricReferenceStress;

        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableMatrix dPressuredStress, dPressuredRCG;

        variableMatrix d2DevStressdStressdRCG, d2PressuredStressdRCG;
        
        errorOut error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( referenceHigherOrderStress,
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
        error = micromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress,
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
        dFdStress = vectorTools::dot( dNormDevStressdDevStress, dDevStressdStress )
                  + BAngle * dPressuredStress;

        dFdc = -AAngle * vectorTools::eye< constantType >( cohesion.size() );

        dFdElasticRCG = vectorTools::dot( dNormDevStressdDevStress, dDevStressdRCG )
                      + BAngle * dPressuredRCG;

        //Construct the second-order jacobians
        d2FdStress2 = variableMatrix( dim, variableVector( dim * dim * dim * dim * dim * dim, 0 ) );
        d2FdStressdElasticRCG = vectorTools::dot( dNormDevStressdDevStress, d2DevStressdStressdRCG )
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
        inversePlasticDeformationGradient = vectorTools::inverse( plasticDeformationGradient, dim, dim );

        inversePlasticMicroDeformation = vectorTools::inverse( plasticMicroDeformation, dim, dim );

        //Assemble the elastic parts of the deformation measures
        elasticDeformationGradient = vectorTools::matrixMultiply( deformationGradient, inversePlasticDeformationGradient,
                                                                  dim, dim, dim, dim );

        elasticMicroDeformation = vectorTools::matrixMultiply( microDeformation, inversePlasticMicroDeformation,
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
        vectorTools::eye( eye );

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

        errorOut error = constitutiveTools::computeRightCauchyGreen( elasticDeformationGradient, elasticRightCauchyGreen );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = constitutiveTools::computeRightCauchyGreen( elasticMicroDeformation, elasticMicroRightCauchyGreen );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic micro right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computePsi( elasticDeformationGradient, elasticMicroDeformation, elasticPsi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic micro-deformation metric Psi" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computeGamma( elasticDeformationGradient, elasticGradientMicroDeformation, elasticGamma );

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

        errorOut error = constitutiveTools::computeRightCauchyGreen( elasticDeformationGradient, elasticRightCauchyGreen,
                                                                     dElasticRCGdElasticF );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = constitutiveTools::computeRightCauchyGreen( elasticMicroDeformation, elasticMicroRightCauchyGreen,
                                                            dElasticMicroRCGdElasticChi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic micro right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computePsi( elasticDeformationGradient, elasticMicroDeformation, elasticPsi,
                                               dElasticPsidElasticF, dElasticPsidElasticChi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic micro-deformation metric Psi" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computeGamma( elasticDeformationGradient, elasticGradientMicroDeformation, elasticGamma,
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
        vectorTools::eye( eye );

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
                                 * elasticPsi[ dim * Nb + Eb ]
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
                                 * elasticPsi[ dim * Nb + Eb ]
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
        vectorTools::eye( eye );

        //Assemble the Jacobians
        dPlasticMicroLdElasticMicroRCG = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroLdElasticPsi = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroLdMicroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int Eb = 0; Eb < dim; Eb++ ){
            for ( unsigned int Fb = 0; Fb < dim; Fb++ ){
                for ( unsigned int Ob = 0; Ob < dim; Ob++ ){
                    for ( unsigned int Pb = 0; Pb < dim; Pb++ ){
                        dPlasticMicroLdElasticPsi[ dim * Eb + Fb ][ dim * Ob + Pb ]
                            -= inverseElasticPsi[ dim * Eb + Ob ] * plasticMicroVelocityGradient[ dim * Pb + Fb ];

                        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                            dPlasticMicroLdElasticPsi[ dim * Eb + Fb ][ dim * Ob + Pb ]
                                += microGamma * inverseElasticPsi[ dim * Eb + Lb ]
                                 * microFlowDirection[ dim * Pb + Lb ]
                                 * elasticMicroRightCauchyGreen[ dim * Ob + Fb ];

                            dPlasticMicroLdMicroFlowDirection[ dim * Eb + Fb ][ dim * Ob + Pb ]
                                += microGamma * inverseElasticPsi[ dim * Eb + Pb ]
                                 * elasticPsi[ dim * Lb + Ob ]
                                 * elasticMicroRightCauchyGreen[ dim * Lb + Fb ];

                            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                                dPlasticMicroLdElasticMicroRCG[ dim * Eb + Fb ][ dim * Ob + Pb ]
                                    += microGamma * inverseElasticPsi[ dim * Eb + Lb ]
                                     * microFlowDirection[ dim * Kb + Lb ] * elasticPsi[ dim * Ob + Kb ]
                                     * eye[ dim * Fb + Pb ];
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
        vectorTools::eye( eye );

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
        vectorTools::eye( eye );

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
        variableVector inverseElasticRightCauchyGreen = vectorTools::inverse( elasticRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = vectorTools::inverse( elasticPsi, dim, dim );

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
        variableVector inverseElasticRightCauchyGreen = vectorTools::inverse( elasticRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = vectorTools::inverse( elasticPsi, dim, dim );

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

        dPlasticMicroGradientLdMicroGamma = vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
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
        variableVector inverseElasticRightCauchyGreen = vectorTools::inverse( elasticRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = vectorTools::inverse( elasticPsi, dim, dim );

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

        dPlasticMicroGradientLdMicroGamma = vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                              dPlasticMicroLdMicroGamma );

        dPlasticMicroGradientLdElasticMicroRCG = vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                                   dPlasticMicroLdElasticMicroRCG );

        dPlasticMicroGradientLdElasticPsi += vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                               dPlasticMicroLdElasticPsi );

        dPlasticMicroGradientLdMicroFlowDirection = vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
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
        vectorTools::eye( eye );

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
        currentPlasticMicroGradient = vectorTools::solveLinearSystem( LHS, RHS, rank );

        if ( rank != LHS.size() ){
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
        vectorTools::eye( eye );

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

        variableVector floatLHS = vectorTools::appendVectors( LHS );

        Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > LHSMat( floatLHS.data(), LHS.size(), LHS.size() );
        Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCDA( negdRdCurrentDtAtilde.data(), LHS.size(), dim * dim * dim );
        Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCFA( negdRdCurrentFourthA.data(), LHS.size(), dim * dim * dim * dim );

        Eigen::ColPivHouseholderQR< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > qrSolver( LHSMat );

        Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X1( vecdCurrentPlasticMicroGradientdCurrentDTAtilde.data(), LHS.size(), dim * dim * dim );
        Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X2( vecdCurrentPlasticMicroGradientdCurrentFourthA.data(), LHS.size(), dim * dim * dim * dim );

        X1 = qrSolver.solve( nDRDCDA );
        X2 = qrSolver.solve( nDRDCFA );

        variableMatrix dCurrentPlasticMicroGradientdCurrentDTAtilde = vectorTools::inflate( vecdCurrentPlasticMicroGradientdCurrentDTAtilde, dim * dim * dim, dim * dim * dim );
        variableMatrix dCurrentPlasticMicroGradientdCurrentFourthA = vectorTools::inflate( vecdCurrentPlasticMicroGradientdCurrentFourthA, dim * dim * dim, dim * dim * dim * dim );

        //Assemble the final terms of the deformation
        dCurrentPlasticMicroGradientdPlasticMicroDeformation = vectorTools::dot( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                 dCurrentDTAtildedPlasticMicroDeformation );

        dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient = vectorTools::dot( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                      dCurrentFourthAdMacroVelocityGradient );

        dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient = vectorTools::dot( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                      dCurrentFourthAdMicroVelocityGradient );

        dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient = vectorTools::dot( dCurrentPlasticMicroGradientdCurrentDTAtilde,
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

        errorOut error = constitutiveTools::evolveF( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
                                                     currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
                                                     alphaMacro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation",
                                             "Error in computation of the plastic macro deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = constitutiveTools::evolveF( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
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

        errorOut error = constitutiveTools::evolveF( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
                                                     currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
                                                     dPlasticFdPlasticMacroL, alphaMacro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation (jacobian)",
                                             "Error in computation of the plastic macro deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = constitutiveTools::evolveF( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
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

        dPlasticMicroGradientdPlasticMicroL += vectorTools::dot( dPlasticMicroGradientdPlasticMicroDeformation,
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

        //Evolve the macro-scale internal strain-like state variable
        currentMacroStrainISV = previousMacroStrainISV + Dt * (         - alphaMacro * previousMacroGamma * previousdMacroGdMacroC
                                                                - ( 1 - alphaMacro ) *  currentMacroGamma * currentdMacroGdMacroC );

        //Evolve the micro-scale internal strain-like state variable
        currentMicroStrainISV = previousMicroStrainISV + Dt * (         - alphaMicro * previousMicroGamma * previousdMicroGdMicroC
                                                                - ( 1 - alphaMicro ) *  currentMicroGamma * currentdMicroGdMicroC );

        //Evolve the micro gradient internal strain-like state variables
        currentMicroGradientStrainISV = previousMicroGradientStrainISV
        + Dt * ( 
            - alphaMicroGradient * vectorTools::Tdot( previousdMicroGradientGdMicroGradientC, previousMicroGradientGamma )
            - ( 1 - alphaMicroGradient ) * vectorTools::Tdot( currentdMicroGradientGdMicroGradientC, currentMicroGradientGamma )
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
        vectorTools::eye( eye );

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

    errorOut computePlasticDeformationResidual( const solverTools::floatVector &x, const solverTools::floatMatrix &floatArgs,
                                                const solverTools::intMatrix &intArgs, solverTools::floatVector &residual,
                                                solverTools::floatMatrix &jacobian, solverTools::floatMatrix &floatOuts,
                                                solverTools::intMatrix &intOuts
#ifdef DEBUG_MODE
                                                , solverTools::debugMap &DEBUG
#endif
                                              ){
        /*!
         * Compute the residual on the amount of plastic deformation.
         * 
         * :param solverTools::floatVector &x: The unknown vector. Organized as
         *     [ plasticDeformationGradient, plasticMicroDeformation, plasticGradientMicroDeformation,
         *       plasticMacroStrainISV, plasticMicroStrainISV, plasticMicroGradientStrainISV,
         *       currentMacroGamma, currentMicroGamma, currentMicroGradientGamma ]
         * :param const solverTools::floatMatrix &floatArgs: The floating point arguments which do not vary
         *     during the solve.
         * :param const solverTools::intMatrix &intArgs: The integer arguments which do not vary during 
         *     the solve.
         * :param solverTools::floatVector &residual: The value of the residual. This will be the 
         *     the values in the x vector - the estimated amount of plastic deformation
         * :param solverTools::floatMatrix &jacobian: The jacobian matrix
         * :param solverTools::floatMatrix &floatOuts: The floating point values that do change during the solve.
         * :param solverTools::intMatrix &intOuts: The integer values that do change during the solve.
         * :param std::map< std::string, solverTools::floatVector > &DEBUG: The debug map. Only available if
         *     DEBUG_MODE is defined.
         *
         * Ordering of floatArgs
         * floatArgs[  0 ] = Dt, The change in time
         * floatArgs[  1 ] = currentDeformationGradient, The current value of the deformation gradient
         * floatArgs[  2 ] = currentMicroDeformation, The current value of the micro-deformation
         * floatArgs[  3 ] = currentGradientMicroDeformation, The current value of the gradient of the
         *     micro deformation in the reference configuration.
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
         * floatArgs[ 25 ] = macroYieldParameters, The macro yield parameters.
         * floatArgs[ 26 ] = microYieldParameters, The micro yield parameters.
         * floatArgs[ 27 ] = microGradientYieldParameters, The micro gradient yield parameters.
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
         * floatOuts[ 0 ] = currentPK2Stress, The current value of the second Piola-Kirchoff stress
         * floatOuts[ 1 ] = currentReferenceMicroStress, The current value of the reference micro
         *     stress.
         * floatOuts[ 2 ] = currentReferenceHigherOrderStress, The current value of the reference
         *     higher order stress.
         * floatOuts[ 3 ] = dStressdx, The partial derivative of the stresses w.r.t. the 
         *     solution vector.
         * floatOuts[ 4 ] = dStressdDeformation, The partial derivative of the stresses w.r.t.
         *     the fundamental deformation measures.
         * floatOuts[ 5 ] = dResidualdDeformation, The partial derivative of the residual w.r.t.
         *     the fundamental deformation measures.
         *
         * Ordering of intOuts
         * intOuts[ 0 ] = isYielding, The flags which indicate if a yield surface is yielding or not.
         *     This activates different parts of the residual and can change during the solution if
         *     different surfaces begin yielding. Once a surface has begun yielding it is assumed to 
         *     continue yielding.
         */

        if ( x.size() != 55 ){
            return new errorNode( "computePlasticDeformationResidual",
                                  "The x vector must have a length of 55" );
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
            if ( floatOuts.size() != 6 ){
                return new errorNode( "computePlasticDeformationResidual",
                                      "The floating point output matrix floatOuts must have a length of 6" );
            }
        }
        else{ 
            if ( floatOuts.size() != 3 ){
                return new errorNode( "computePlasticDeformationResidual",
                                      "The floating point output matrix floatOuts must have a length of 3" );
            }
        }

        if ( intOuts.size() != 1 ){
            return new errorNode( "computePlasticDeformationResidual",
                                  "The integer output matrix intOuts must have a length of 1" );
        }

        /*=============================
        | Extract the incoming values |
        =============================*/

        const variableVector currentPlasticDeformationGradient( x.begin(), x.begin() + 9 );
        const variableVector currentPlasticMicroDeformation( x.begin() + 9, x.begin() + 18 );
        const variableVector currentPlasticGradientMicroDeformation( x.begin() + 18, x.begin() + 45 );
        const variableType   currentMacroStrainISV = x[ 45 ];
        const variableType   currentMicroStrainISV = x[ 46 ];
        const variableVector currentMicroGradientStrainISV( x.begin() + 47, x.begin() + 50 );
        const variableType   currentMacroGamma = x[ 50 ];
        const variableType   currentMicroGamma = x[ 51 ];
        const variableVector currentMicroGradientGamma( x.begin() + 52, x.begin() + 55 );

        //Assume the logs of the gammas are passed in ( enforcing positivity )
        currentMacroGamma = std::exp( currentMacroGamma );
        currentMicroGamma = std::exp( currentMicroGamma );
        for ( unsigned int i = 0; i < currentMicroGradientGamma.size(); i++ ){
            currentMicroGradientGamma[ i ] = std::exp( currentMicroGradientGamma[ i ] );
        }

#ifdef DEBUG_MODE
        variableVector temp = { currentMacroGamma };
        DEBUG.emplace( "currentMacroGamma", temp );

        temp = { currentMicroGamma };
        DEBUG.emplace( "currentMicroGamma", temp );

        DEBUG.emplace( "currentMicroGradientGamma", currentMicroGradientGamma );
#endif

        unsigned int ii = 0;
        const constantType    *Dt                                            = &floatArgs[ ii++ ][ 0 ];
        const variableVector  *currentDeformationGradient                    = &floatArgs[ ii++ ];
        const variableVector  *currentMicroDeformation                       = &floatArgs[ ii++ ];
        const variableVector  *currentGradientMicroDeformation               = &floatArgs[ ii++ ];
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
        const variableMatrix   previousdMicroGradientGdMicroGradientCohesion = vectorTools::inflate( floatArgs[ ii++ ], 3, 3 );
        const variableVector  *previousPlasticMacroVelocityGradient          = &floatArgs[ ii++ ];
        const variableVector  *previousPlasticMicroVelocityGradient          = &floatArgs[ ii++ ];
        const variableVector  *previousPlasticMicroGradientVelocityGradient  = &floatArgs[ ii++ ];
        const parameterVector *macroFlowParameters                           = &floatArgs[ ii++ ];
        const parameterVector *microFlowParameters                           = &floatArgs[ ii++ ];
        const parameterVector *microGradientFlowParameters                   = &floatArgs[ ii++ ];
        const parameterVector *macroHardeningParameters                      = &floatArgs[ ii++ ];
        const parameterVector *microHardeningParameters                      = &floatArgs[ ii++ ];
        const parameterVector *microGradientHardeningParameters              = &floatArgs[ ii++ ];
        const parameterVector *macroYieldParameters                          = &floatArgs[ ii++ ];
        const parameterVector *microYieldParameters                          = &floatArgs[ ii++ ];
        const parameterVector *microGradientYieldParameters                  = &floatArgs[ ii++ ];
        const parameterVector *Amatrix                                       = &floatArgs[ ii++ ];
        const parameterVector *Bmatrix                                       = &floatArgs[ ii++ ];
        const parameterVector *Cmatrix                                       = &floatArgs[ ii++ ];
        const parameterVector *Dmatrix                                       = &floatArgs[ ii++ ];
        const parameterType   *alphaMacro                                    = &floatArgs[ ii++ ][ 0 ];
        const parameterType   *alphaMicro                                    = &floatArgs[ ii++ ][ 0 ];
        const parameterType   *alphaMicroGradient                            = &floatArgs[ ii++ ][ 0 ];

        ii = 0;
        solverTools::intVector isYielding                              = intOuts[ ii++ ];

        if ( isYielding.size() != 5 ){
            return new errorNode( "computePlasticDeformationResidual",
                                  "The yielding flags vector (isYielding) must have a size of 5" );
        }

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

#ifdef DEBUG_MODE

        //Save the elastic fundamental deformation measures
        DEBUG.emplace( "currentElasticDeformationGradient", currentElasticDeformationGradient );
        DEBUG.emplace( "currentElasticMicroDeformation", currentElasticMicroDeformation );
        DEBUG.emplace( "currentElasticGradientMicroDeformation", currentElasticGradientMicroDeformation );

        //Save the flattened jacobians
        DEBUG.emplace( "dElasticDeformationGradientdDeformationGradient",
                       vectorTools::appendVectors( dElasticDeformationGradientdDeformationGradient )  );
        DEBUG.emplace( "dElasticDeformationGradientdPlasticDeformationGradient",
                       vectorTools::appendVectors( dElasticDeformationGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticMicroDeformationdMicroDeformation",
                       vectorTools::appendVectors( dElasticMicroDeformationdMicroDeformation ) );
        DEBUG.emplace( "dElasticMicroDeformationdPlasticMicroDeformation",
                       vectorTools::appendVectors( dElasticMicroDeformationdPlasticMicroDeformation ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdGradientMicroDeformation",
                       vectorTools::appendVectors( dElasticGradientMicroDeformationdGradientMicroDeformation ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dElasticGradientMicroDeformationdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdPlasticDeformationGradient",
                       vectorTools::appendVectors( dElasticGradientMicroDeformationdPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdMicroDeformation",
                       vectorTools::appendVectors( dElasticGradientMicroDeformationdMicroDeformation ) );
        DEBUG.emplace( "dElasticGradientMicroDeformationdPlasticMicroDeformation",
                       vectorTools::appendVectors( dElasticGradientMicroDeformationdPlasticMicroDeformation ) );

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

        /*!==================================================================
        | Assemble the Jacobian of the elastic derived deformation measures |
        ===================================================================*/

        //Compute the Jacobians w.r.t. the plastic deformation
        variableMatrix dElasticRightCauchyGreendPlasticDeformationGradient
            = vectorTools::dot( dElasticRightCauchyGreendElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient );

        variableMatrix dElasticMicroRightCauchyGreendPlasticMicroDeformation
            = vectorTools::dot( dElasticMicroRightCauchyGreendElasticMicroDeformation,
                                dElasticMicroDeformationdPlasticMicroDeformation );

        variableMatrix dElasticPsidPlasticDeformationGradient
            = vectorTools::dot( dElasticPsidElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient );

        variableMatrix dElasticPsidPlasticMicroDeformation
            = vectorTools::dot( dElasticPsidElasticMicroDeformation,
                                dElasticMicroDeformationdPlasticMicroDeformation );

        variableMatrix dElasticGammadPlasticDeformationGradient
            = vectorTools::dot( dElasticGammadElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient )
            + vectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticDeformationGradient );

        variableMatrix dElasticGammadPlasticMicroDeformation
            = vectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticMicroDeformation );

        variableMatrix dElasticGammadPlasticGradientMicroDeformation
            = vectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticGradientMicroDeformation );

        //Construct the full derivatives
        variableMatrix dElasticRightCauchyGreendDeformationGradient, dElasticMicroRightCauchyGreendMicroDeformation,
                       dElasticPsidDeformationGradient, dElasticPsidMicroDeformation,
                       dElasticGammadDeformationGradient, dElasticGammadMicroDeformation, dElasticGammadGradientMicroDeformation;
        if ( evaluateFullDerivatives ){
            dElasticRightCauchyGreendDeformationGradient = vectorTools::dot( dElasticRightCauchyGreendElasticDeformationGradient,
                                                                             dElasticDeformationGradientdDeformationGradient );
            dElasticMicroRightCauchyGreendMicroDeformation = vectorTools::dot( dElasticMicroRightCauchyGreendElasticMicroDeformation,
                                                                               dElasticMicroDeformationdMicroDeformation );
            dElasticPsidDeformationGradient = vectorTools::dot( dElasticPsidElasticDeformationGradient,
                                                                dElasticDeformationGradientdDeformationGradient );
            dElasticPsidMicroDeformation = vectorTools::dot( dElasticPsidElasticMicroDeformation,
                                                             dElasticMicroDeformationdMicroDeformation );
            dElasticGammadDeformationGradient = vectorTools::dot( dElasticGammadElasticDeformationGradient,
                                                                  dElasticDeformationGradientdDeformationGradient );
            dElasticGammadMicroDeformation = vectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                                               dElasticGradientMicroDeformationdMicroDeformation );
            dElasticGammadGradientMicroDeformation = vectorTools::dot( dElasticGammadElasticGradientMicroDeformation,
                                                                       dElasticGradientMicroDeformationdGradientMicroDeformation );
        }

#ifdef DEBUG_MODE

        //Save the elastic derived deformation measures
        DEBUG.emplace( "currentElasticRightCauchyGreen", currentElasticRightCauchyGreen );
        DEBUG.emplace( "currentElasticMicroRightCauchyGreen", currentElasticMicroRightCauchyGreen );
        DEBUG.emplace( "currentElasticPsi", currentElasticPsi );
        DEBUG.emplace( "currentElasticGamma", currentElasticGamma );

        //Save the Jacobians
        DEBUG.emplace( "dElasticRightCauchyGreendPlasticDeformationGradient",
                       vectorTools::appendVectors( dElasticRightCauchyGreendPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticMicroRightCauchyGreendPlasticMicroDeformation",
                       vectorTools::appendVectors( dElasticMicroRightCauchyGreendPlasticMicroDeformation ) );
        DEBUG.emplace( "dElasticPsidPlasticDeformationGradient",
                       vectorTools::appendVectors( dElasticPsidPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticPsidPlasticMicroDeformation",
                       vectorTools::appendVectors( dElasticPsidPlasticMicroDeformation ) );
        DEBUG.emplace( "dElasticGammadPlasticDeformationGradient",
                       vectorTools::appendVectors( dElasticGammadPlasticDeformationGradient ) );
        DEBUG.emplace( "dElasticGammadPlasticMicroDeformation",
                       vectorTools::appendVectors( dElasticGammadPlasticMicroDeformation ) );
        DEBUG.emplace( "dElasticGammadPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dElasticGammadPlasticGradientMicroDeformation ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dElasticRightCauchyGreendDeformationGradient",
                            vectorTools::appendVectors( dElasticRightCauchyGreendDeformationGradient ) );
            DEBUG.emplace( "dElasticMicroRightCauchyGreendMicroDeformation",
                            vectorTools::appendVectors( dElasticMicroRightCauchyGreendMicroDeformation ) );
            DEBUG.emplace( "dElasticPsidDeformationGradient",
                            vectorTools::appendVectors( dElasticPsidDeformationGradient ) );
            DEBUG.emplace( "dElasticPsidMicroDeformation",
                            vectorTools::appendVectors( dElasticPsidMicroDeformation ) );
            DEBUG.emplace( "dElasticGammadDeformationGradient",
                            vectorTools::appendVectors( dElasticGammadDeformationGradient ) );
            DEBUG.emplace( "dElasticGammadMicroDeformation",
                            vectorTools::appendVectors( dElasticGammadMicroDeformation ) );
            DEBUG.emplace( "dElasticGammadGradientMicroDeformation",
                            vectorTools::appendVectors( dElasticGammadGradientMicroDeformation ) );
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

        error = micromorphicLinearElasticity::linearElasticityReference(  currentElasticDeformationGradient,
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

        /*===============================
        | Assemble the stress Jacobians |
        ===============================*/

        //Jacobians w.r.t. the plastic deformation
        variableMatrix dPK2StressdPlasticDeformationGradient
            = vectorTools::dot( dPK2StressdElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient )
            + vectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticDeformationGradient );

        variableMatrix dPK2StressdPlasticMicroDeformation
            = vectorTools::dot( dPK2StressdElasticMicroDeformation,
                                dElasticMicroDeformationdPlasticMicroDeformation )
            + vectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticMicroDeformation );

        variableMatrix dPK2StressdPlasticGradientMicroDeformation
            = vectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticGradientMicroDeformation );

        variableMatrix dReferenceMicroStressdPlasticDeformationGradient
            = vectorTools::dot( dReferenceMicroStressdElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient )
            + vectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticDeformationGradient );

        variableMatrix dReferenceMicroStressdPlasticMicroDeformation
            = vectorTools::dot( dReferenceMicroStressdElasticMicroDeformation,
                                dElasticMicroDeformationdPlasticMicroDeformation )
            + vectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticMicroDeformation );

        variableMatrix dReferenceMicroStressdPlasticGradientMicroDeformation
            = vectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticGradientMicroDeformation );

        variableMatrix dReferenceHigherOrderStressdPlasticDeformationGradient
            = vectorTools::dot( dReferenceHigherOrderStressdElasticDeformationGradient,
                                dElasticDeformationGradientdPlasticDeformationGradient )
            + vectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticDeformationGradient );

        variableMatrix dReferenceHigherOrderStressdPlasticMicroDeformation
            = vectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticMicroDeformation );

        variableMatrix dReferenceHigherOrderStressdPlasticGradientMicroDeformation
            = vectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                dElasticGradientMicroDeformationdPlasticGradientMicroDeformation );

        //Jacobians w.r.t. the fundamental deformation measures
        variableMatrix dPK2StressdDeformationGradient, dPK2StressdMicroDeformation, dPK2StressdGradientMicroDeformation,
                       dReferenceMicroStressdDeformationGradient, dReferenceMicroStressdMicroDeformation,
                       dReferenceMicroStressdGradientMicroDeformation,
                       dReferenceHigherOrderStressdDeformationGradient, dReferenceHigherOrderStressdMicroDeformation,
                       dReferenceHigherOrderStressdGradientMicroDeformation;

        if ( evaluateFullDerivatives ){
            dPK2StressdDeformationGradient = vectorTools::dot( dPK2StressdElasticDeformationGradient,
                                                               dElasticDeformationGradientdDeformationGradient );
            dPK2StressdMicroDeformation = vectorTools::dot( dPK2StressdElasticMicroDeformation,
                                                            dElasticMicroDeformationdMicroDeformation )
                                        + vectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                                            dElasticGradientMicroDeformationdMicroDeformation );
            dPK2StressdGradientMicroDeformation = vectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            dReferenceMicroStressdDeformationGradient = vectorTools::dot( dReferenceMicroStressdElasticDeformationGradient,
                                                                          dElasticDeformationGradientdDeformationGradient );
            dReferenceMicroStressdMicroDeformation = vectorTools::dot( dReferenceMicroStressdElasticMicroDeformation,
                                                                       dElasticMicroDeformationdMicroDeformation )
                                                   + vectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                                                       dElasticGradientMicroDeformationdMicroDeformation );
            dReferenceMicroStressdGradientMicroDeformation
                = vectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            dReferenceHigherOrderStressdDeformationGradient
                = vectorTools::dot( dReferenceHigherOrderStressdElasticDeformationGradient,
                                    dElasticDeformationGradientdDeformationGradient );

            dReferenceHigherOrderStressdMicroDeformation
                = vectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdMicroDeformation );

            dReferenceHigherOrderStressdGradientMicroDeformation
                = vectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            //Assemble the jacobians into the output vector
            
            floatOuts[ 3 ] = solverTools::floatVector( 45 * 55, 0 ); //Jacobians w.r.t. the solution vector x
            floatOuts[ 4 ] = solverTools::floatVector( 45 * 45, 0 ); //Jacobians w.r.t. the fundamental deformation measures

            //Save the Jacobians of the PK2 and reference symmetric micro stresses
            for ( unsigned int i = 0; i < 9; i++ ){
                for ( unsigned int j = 0; j < 9; j++ ){
                    floatOuts[ 3 ][ 55 * i + j ]     = dPK2StressdPlasticDeformationGradient[ i ][ j ];
                    floatOuts[ 4 ][ 45 * i + j ]     = dPK2StressdDeformationGradient[ i ][ j ];
                    floatOuts[ 3 ][ 55 * i + j + 9 ] = dPK2StressdPlasticMicroDeformation[ i ][ j ];
                    floatOuts[ 4 ][ 45 * i + j + 9 ] = dPK2StressdMicroDeformation[ i ][ j ];

                    floatOuts[ 3 ][ 55 * ( i + 9 ) + j ] = dReferenceMicroStressdPlasticDeformationGradient[ i ][ j ];
                    floatOuts[ 4 ][ 45 * ( i + 9 ) + j ] = dReferenceMicroStressdDeformationGradient[ i ][ j ];
                    floatOuts[ 3 ][ 55 * ( i + 9 ) + j + 9 ] = dReferenceMicroStressdPlasticMicroDeformation[ i ][ j ];
                    floatOuts[ 4 ][ 45 * ( i + 9 ) + j + 9 ] = dReferenceMicroStressdMicroDeformation[ i ][ j ];
                }

                for ( unsigned int j = 0; j < 27; j++ ){
                    floatOuts[ 3 ][ 55 * i + j + 18 ] = dPK2StressdPlasticGradientMicroDeformation[ i ][ j ];
                    floatOuts[ 4 ][ 45 * i + j + 18 ] = dPK2StressdGradientMicroDeformation[ i ][ j ];

                    floatOuts[ 3 ][ 55 * ( i + 9 ) + j + 18 ] = dReferenceMicroStressdPlasticGradientMicroDeformation[ i ][ j ];
                    floatOuts[ 4 ][ 45 * ( i + 9 ) + j + 18 ] = dReferenceMicroStressdGradientMicroDeformation[ i ][ j ];
                }
            }

            //Save the Jacobians of the reference higher order stress
            for ( unsigned int i = 0; i < 27; i++ ){
                for ( unsigned int j = 0; j < 9; j++ ){
                    floatOuts[ 3 ][ 55 * ( i + 18 ) + j ] = dReferenceHigherOrderStressdPlasticDeformationGradient[ i ][ j ];
                    floatOuts[ 4 ][ 45 * ( i + 18 ) + j ] = dReferenceHigherOrderStressdDeformationGradient[ i ][ j ];
                    floatOuts[ 3 ][ 55 * ( i + 18 ) + j + 9 ] = dReferenceHigherOrderStressdPlasticMicroDeformation[ i ][ j ];
                    floatOuts[ 4 ][ 45 * ( i + 18 ) + j + 9 ] = dReferenceHigherOrderStressdMicroDeformation[ i ][ j ];
                }

                for ( unsigned int j = 0; j < 27; j++ ){
                    floatOuts[ 3 ][ 55 * ( i + 18 ) + j + 18 ] = dReferenceHigherOrderStressdPlasticGradientMicroDeformation[ i ][ j ];
                    floatOuts[ 4 ][ 45 * ( i + 18 ) + j + 18 ] = dReferenceHigherOrderStressdGradientMicroDeformation[ i ][ j ];
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
                       vectorTools::appendVectors( dPK2StressdPlasticDeformationGradient ) );
        DEBUG.emplace( "dPK2StressdPlasticMicroDeformation",
                       vectorTools::appendVectors( dPK2StressdPlasticMicroDeformation ) );
        DEBUG.emplace( "dPK2StressdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dPK2StressdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dReferenceMicroStressdPlasticDeformationGradient",
                       vectorTools::appendVectors( dReferenceMicroStressdPlasticDeformationGradient ) );
        DEBUG.emplace( "dReferenceMicroStressdPlasticMicroDeformation",
                       vectorTools::appendVectors( dReferenceMicroStressdPlasticMicroDeformation ) );
        DEBUG.emplace( "dReferenceMicroStressdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dReferenceMicroStressdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dReferenceHigherOrderStressdPlasticDeformationGradient",
                       vectorTools::appendVectors( dReferenceHigherOrderStressdPlasticDeformationGradient ) );
        DEBUG.emplace( "dReferenceHigherOrderStressdPlasticMicroDeformation",
                       vectorTools::appendVectors( dReferenceHigherOrderStressdPlasticMicroDeformation ) );
        DEBUG.emplace( "dReferenceHigherOrderStressdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dReferenceHigherOrderStressdPlasticGradientMicroDeformation ) );

        //Save the Jacobians w.r.t. the fundamental deformation measures
        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dPK2StressdDeformationGradient",
                            vectorTools::appendVectors( dPK2StressdDeformationGradient ) );
            DEBUG.emplace( "dPK2StressdMicroDeformation",
                            vectorTools::appendVectors( dPK2StressdMicroDeformation ) );
            DEBUG.emplace( "dPK2StressdGradientMicroDeformation",
                            vectorTools::appendVectors( dPK2StressdGradientMicroDeformation ) );
    
            DEBUG.emplace( "dReferenceMicroStressdDeformationGradient",
                            vectorTools::appendVectors( dReferenceMicroStressdDeformationGradient ) );
            DEBUG.emplace( "dReferenceMicroStressdMicroDeformation",
                            vectorTools::appendVectors( dReferenceMicroStressdMicroDeformation ) );
            DEBUG.emplace( "dReferenceMicroStressdGradientMicroDeformation",
                            vectorTools::appendVectors( dReferenceMicroStressdGradientMicroDeformation ) );
    
            DEBUG.emplace( "dReferenceHigherOrderStressdDeformationGradient",
                            vectorTools::appendVectors( dReferenceHigherOrderStressdDeformationGradient ) );
            DEBUG.emplace( "dReferenceHigherOrderStressdMicroDeformation",
                            vectorTools::appendVectors( dReferenceHigherOrderStressdMicroDeformation ) );
            DEBUG.emplace( "dReferenceHigherOrderStressdGradientMicroDeformation",
                            vectorTools::appendVectors( dReferenceHigherOrderStressdGradientMicroDeformation ) );
        }

#endif

        /*!============================
        | Compute the cohesion values |
        =============================*/

        variableType   currentMacroCohesion, currentMicroCohesion;
        variableVector currentMicroGradientCohesion;

        variableType   dMacroCohesiondMacroStrainISV, dMicroCohesiondMicroStrainISV;
        variableMatrix dMicroGradientCohesiondMicroGradientStrainISV;

        error = computeCohesion(  currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                 *macroHardeningParameters, *microHardeningParameters, *microGradientHardeningParameters,
                                  currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion,
                                  dMacroCohesiondMacroStrainISV, dMicroCohesiondMicroStrainISV,
                                  dMicroGradientCohesiondMicroGradientStrainISV );
        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the cohesion values" );
            result->addNext( error );
            return result;
        }

#ifdef DEBUG_MODE

        //Save the cohesion values
        temp = { currentMacroCohesion };

        DEBUG.emplace( "currentMacroCohesion", temp );

        temp = { currentMicroCohesion };

        DEBUG.emplace( "currentMicroCohesion", temp );

        DEBUG.emplace( "currentMicroGradientCohesion", currentMicroGradientCohesion );

        //Save the cohesion Jacobians
        temp = { dMacroCohesiondMacroStrainISV };

        DEBUG.emplace( "dMacroCohesiondMacroStrainISV", temp );

        temp = { dMicroCohesiondMicroStrainISV };

        DEBUG.emplace( "dMicroCohesiondMicroStrainISV", temp );

        DEBUG.emplace( "dMicroGradientCohesiondMicroGradientStrainISV",
                       vectorTools::appendVectors( dMicroGradientCohesiondMicroGradientStrainISV ) );

#endif


        /*!============================
        | Compute the Flow Directions |
        =============================*/

        variableVector currentMacroFlowDirection, currentMicroFlowDirection, currentMicroGradientFlowDirection;
        variableType currentdMacroGdMacroCohesion, currentdMicroGdMicroCohesion;
        variableMatrix currentdMicroGradientGdMicroGradientCohesion;

        variableMatrix dMacroFlowDirectiondPK2Stress, dMacroFlowDirectiondElasticRightCauchyGreen,
                       dMicroFlowDirectiondReferenceMicroStress, dMicroFlowDirectiondElasticRightCauchyGreen,
                       dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                       dMicroGradientFlowDirectiondElasticRightCauchyGreen;

        error = computeFlowDirections(  currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress,
                                        currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion,
                                        currentElasticRightCauchyGreen,
                                       *macroFlowParameters, *microFlowParameters, *microGradientFlowParameters,
                                        currentMacroFlowDirection, currentMicroFlowDirection, currentMicroGradientFlowDirection,
                                        currentdMacroGdMacroCohesion, currentdMicroGdMicroCohesion,
                                        currentdMicroGradientGdMicroGradientCohesion,
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

        /*!============================
        | Assemble the flow Jacobians |
        =============================*/

        //Assemble the Jacobians w.r.t. the plastic deformation
        variableMatrix dMacroFlowDirectiondPlasticDeformationGradient
            = vectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdPlasticDeformationGradient )
            + vectorTools::dot( dMacroFlowDirectiondElasticRightCauchyGreen, dElasticRightCauchyGreendPlasticDeformationGradient );

        variableMatrix dMacroFlowDirectiondPlasticMicroDeformation
            = vectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdPlasticMicroDeformation );

        variableMatrix dMacroFlowDirectiondPlasticGradientMicroDeformation
            = vectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdPlasticGradientMicroDeformation );

        variableMatrix dMicroFlowDirectiondPlasticDeformationGradient
            = vectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdPlasticDeformationGradient )
            + vectorTools::dot( dMicroFlowDirectiondElasticRightCauchyGreen, dElasticRightCauchyGreendPlasticDeformationGradient );

        variableMatrix dMicroFlowDirectiondPlasticMicroDeformation
            = vectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdPlasticMicroDeformation );

        variableMatrix dMicroFlowDirectiondPlasticGradientMicroDeformation
            = vectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdPlasticGradientMicroDeformation );

        variableMatrix dMicroGradientFlowDirectiondPlasticDeformationGradient
            = vectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                dReferenceHigherOrderStressdPlasticDeformationGradient )
            + vectorTools::dot( dMicroGradientFlowDirectiondElasticRightCauchyGreen,
                                dElasticRightCauchyGreendPlasticDeformationGradient );

        variableMatrix dMicroGradientFlowDirectiondPlasticMicroDeformation
            = vectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                dReferenceHigherOrderStressdPlasticMicroDeformation );

        variableMatrix dMicroGradientFlowDirectiondPlasticGradientMicroDeformation
            = vectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
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
                = vectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdDeformationGradient )
                + vectorTools::dot( dMacroFlowDirectiondElasticRightCauchyGreen, dElasticRightCauchyGreendDeformationGradient );

            dMacroFlowDirectiondMicroDeformation
                = vectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdMicroDeformation );

            dMacroFlowDirectiondGradientMicroDeformation
                = vectorTools::dot( dMacroFlowDirectiondPK2Stress, dPK2StressdGradientMicroDeformation );

            dMicroFlowDirectiondDeformationGradient
                = vectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdDeformationGradient )
                + vectorTools::dot( dMicroFlowDirectiondElasticRightCauchyGreen, dElasticRightCauchyGreendDeformationGradient );

            dMicroFlowDirectiondMicroDeformation
                = vectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdMicroDeformation );

            dMicroFlowDirectiondGradientMicroDeformation
                = vectorTools::dot( dMicroFlowDirectiondReferenceMicroStress, dReferenceMicroStressdGradientMicroDeformation );

            dMicroGradientFlowDirectiondDeformationGradient
                = vectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                    dReferenceHigherOrderStressdDeformationGradient )
                + vectorTools::dot( dMicroGradientFlowDirectiondElasticRightCauchyGreen,
                                    dElasticRightCauchyGreendDeformationGradient );

            dMicroGradientFlowDirectiondMicroDeformation
                = vectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                    dReferenceHigherOrderStressdMicroDeformation );

            dMicroGradientFlowDirectiondGradientMicroDeformation
                = vectorTools::dot( dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                    dReferenceHigherOrderStressdGradientMicroDeformation );
        }

#ifdef DEBUG_MODE

        //Save the flow directions
        DEBUG.emplace( "currentMacroFlowDirection", currentMacroFlowDirection );
        DEBUG.emplace( "currentMicroFlowDirection", currentMicroFlowDirection );
        DEBUG.emplace( "currentMicroGradientFlowDirection", currentMicroGradientFlowDirection );

        //Save the Jacobians of the flow directions
        DEBUG.emplace( "dMacroFlowDirectiondPlasticDeformationGradient",
                       vectorTools::appendVectors( dMacroFlowDirectiondPlasticDeformationGradient ) );
        DEBUG.emplace( "dMacroFlowDirectiondPlasticMicroDeformation",
                       vectorTools::appendVectors( dMacroFlowDirectiondPlasticMicroDeformation ) );
        DEBUG.emplace( "dMacroFlowDirectiondPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dMacroFlowDirectiondPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dMicroFlowDirectiondPlasticDeformationGradient",
                       vectorTools::appendVectors( dMicroFlowDirectiondPlasticDeformationGradient ) );
        DEBUG.emplace( "dMicroFlowDirectiondPlasticMicroDeformation",
                       vectorTools::appendVectors( dMicroFlowDirectiondPlasticMicroDeformation ) );
        DEBUG.emplace( "dMicroFlowDirectiondPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dMicroFlowDirectiondPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dMicroGradientFlowDirectiondPlasticDeformationGradient",
                       vectorTools::appendVectors( dMicroGradientFlowDirectiondPlasticDeformationGradient ) );
        DEBUG.emplace( "dMicroGradientFlowDirectiondPlasticMicroDeformation",
                       vectorTools::appendVectors( dMicroGradientFlowDirectiondPlasticMicroDeformation ) );
        DEBUG.emplace( "dMicroGradientFlowDirectiondPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dMicroGradientFlowDirectiondPlasticGradientMicroDeformation  ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dMacroFlowDirectiondDeformationGradient",
                            vectorTools::appendVectors( dMacroFlowDirectiondDeformationGradient ) );
            DEBUG.emplace( "dMacroFlowDirectiondMicroDeformation",
                            vectorTools::appendVectors( dMacroFlowDirectiondMicroDeformation ) );
            DEBUG.emplace( "dMacroFlowDirectiondGradientMicroDeformation",
                            vectorTools::appendVectors( dMacroFlowDirectiondGradientMicroDeformation ) );
            DEBUG.emplace( "dMicroFlowDirectiondDeformationGradient",
                            vectorTools::appendVectors( dMicroFlowDirectiondDeformationGradient ) );
            DEBUG.emplace( "dMicroFlowDirectiondMicroDeformation",
                            vectorTools::appendVectors( dMicroFlowDirectiondMicroDeformation ) );
            DEBUG.emplace( "dMicroFlowDirectiondGradientMicroDeformation",
                            vectorTools::appendVectors( dMicroFlowDirectiondGradientMicroDeformation ) );
            DEBUG.emplace( "dMicroGradientFlowDirectiondDeformationGradient",
                            vectorTools::appendVectors( dMicroGradientFlowDirectiondDeformationGradient ) );
            DEBUG.emplace( "dMicroGradientFlowDirectiondMicroDeformation",
                            vectorTools::appendVectors( dMicroGradientFlowDirectiondMicroDeformation ) );
            DEBUG.emplace( "dMicroGradientFlowDirectiondGradientMicroDeformation",
                            vectorTools::appendVectors( dMicroGradientFlowDirectiondGradientMicroDeformation ) );
        }

#endif

        /*!======================================
        | Compute the expected strain-like ISVs |
        =======================================*/

        variableType expectedMacroStrainISV, expectedMicroStrainISV;
        variableVector expectedMicroGradientStrainISV;

        variableType dExpectedMacroISVdMacroGamma, dExpectedMacroISVddMacroGdMacroCohesion;
        variableType dExpectedMicroISVdMicroGamma, dExpectedMicroISVddMicroGdMicroCohesion;

        variableMatrix dExpectedMicroGradientISVdMicroGradientGamma,
                       dExpectedMicroGradientISVddMicroGradientGdMicroGradientCohesion;

        error = evolveStrainStateVariables( *Dt,
                                             currentMacroGamma, currentMicroGamma, currentMicroGradientGamma,
                                             currentdMacroGdMacroCohesion, currentdMicroGdMicroCohesion,
                                             currentdMicroGradientGdMicroGradientCohesion,
                                            *previousMacroStrainISV, *previousMicroStrainISV, *previousMicroGradientStrainISV,
                                            *previousMacroGamma, *previousMicroGamma, *previousMicroGradientGamma,
                                            *previousdMacroGdMacroCohesion, *previousdMicroGdMicroCohesion,
                                             previousdMicroGradientGdMicroGradientCohesion,
                                             expectedMacroStrainISV, expectedMicroStrainISV, expectedMicroGradientStrainISV,
                                             dExpectedMacroISVdMacroGamma, dExpectedMacroISVddMacroGdMacroCohesion,
                                             dExpectedMicroISVdMicroGamma, dExpectedMicroISVddMicroGdMicroCohesion,
                                             dExpectedMicroGradientISVdMicroGradientGamma,
                                             dExpectedMicroGradientISVddMicroGradientGdMicroGradientCohesion,
                                            *alphaMacro, *alphaMicro, *alphaMicroGradient );
        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the expected strain-like ISVs" );
            result->addNext( error );
            return result;
        }

#ifdef DEBUG_MODE

        //Save the expected macro-strain values
        temp = { expectedMacroStrainISV };
        DEBUG.emplace( "expectedMacroStrainISV", temp );

        temp = { expectedMicroStrainISV };
        DEBUG.emplace( "expectedMicroStrainISV", temp );

        DEBUG.emplace( "expectedMicroGradientStrainISV", expectedMicroGradientStrainISV );

        //Save the Jacobians w.r.t. the Gammas
        temp = { dExpectedMacroISVdMacroGamma };
        DEBUG.emplace( "dExpectedMacroISVdMacroGamma", temp );

        temp = { dExpectedMicroISVdMicroGamma };
        DEBUG.emplace( "dExpectedMicroISVdMicroGamma", temp );

        DEBUG.emplace( "dExpectedMicroGradientISVdMicroGradientGamma",
                        vectorTools::appendVectors( dExpectedMicroGradientISVdMicroGradientGamma ) );

        //Save the Jacobians w.r.t. the derivative of the flow direction w.r.t. the cohesion
        temp = { dExpectedMacroISVddMacroGdMacroCohesion };
        DEBUG.emplace( "dExpectedMacroISVddMacroGdMacroCohesion", temp );

        temp = { dExpectedMicroISVddMicroGdMicroCohesion };
        DEBUG.emplace( "dExpectedMicroISVddMicroGdMicroCohesion", temp );

        DEBUG.emplace( "dExpectedMicroGradientISVddMicroGradientGdMicroGradientCohesion",
                        vectorTools::appendVectors( dExpectedMicroGradientISVddMicroGradientGdMicroGradientCohesion ) );

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

        error = computePlasticVelocityGradients( currentMacroGamma, currentMicroGamma, currentMicroGradientGamma, 
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
            = vectorTools::dot( dPlasticMacroVelocityGradientdElasticRightCauchyGreen,
                                dElasticRightCauchyGreendPlasticDeformationGradient )
            + vectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                dMacroFlowDirectiondPlasticDeformationGradient )
            + vectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticDeformationGradient );

        variableMatrix dPlasticMacroVelocityGradientdPlasticMicroDeformation
            = vectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                dMacroFlowDirectiondPlasticMicroDeformation )
            + vectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticMicroDeformation );

        variableMatrix dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation
            = vectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                dMacroFlowDirectiondPlasticGradientMicroDeformation )
            + vectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticGradientMicroDeformation );

        variableMatrix dPlasticMicroVelocityGradientdPlasticDeformationGradient
            = vectorTools::dot( dPlasticMicroVelocityGradientdElasticPsi,
                                dElasticPsidPlasticDeformationGradient )
            + vectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticDeformationGradient );

        variableMatrix dPlasticMicroVelocityGradientdPlasticMicroDeformation
            = vectorTools::dot( dPlasticMicroVelocityGradientdElasticMicroRightCauchyGreen,
                                dElasticMicroRightCauchyGreendPlasticMicroDeformation )
            + vectorTools::dot( dPlasticMicroVelocityGradientdElasticPsi,
                                dElasticPsidPlasticMicroDeformation )
            + vectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticMicroDeformation );

        variableMatrix dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation
            = vectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticGradientMicroDeformation );

        variableMatrix dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient
            = vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticPsi,
                                dElasticPsidPlasticDeformationGradient )
            + vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                dElasticGammadPlasticDeformationGradient )
            + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticDeformationGradient )
            + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                dMicroGradientFlowDirectiondPlasticDeformationGradient );

        variableMatrix dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation
            = vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticMicroRightCauchyGreen,
                                dElasticMicroRightCauchyGreendPlasticMicroDeformation )
            + vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticPsi,
                                dElasticPsidPlasticMicroDeformation )
            + vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                dElasticGammadPlasticMicroDeformation )
            + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticMicroDeformation )
            + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                dMicroGradientFlowDirectiondPlasticMicroDeformation );

        variableMatrix dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation
            = vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                dElasticGammadPlasticGradientMicroDeformation )
            + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                dMicroFlowDirectiondPlasticGradientMicroDeformation )
            + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                dMicroGradientFlowDirectiondPlasticGradientMicroDeformation );

        variableMatrix dPlasticMacroVelocityGradientdDeformationGradient, dPlasticMacroVelocityGradientdMicroDeformation,
                       dPlasticMacroVelocityGradientdGradientMicroDeformation,
                       dPlasticMicroVelocityGradientdDeformationGradient, dPlasticMicroVelocityGradientdMicroDeformation,
                       dPlasticMicroVelocityGradientdGradientMicroDeformation,
                       dPlasticGradientMicroVelocityGradientdDeformationGradient, dPlasticGradientMicroVelocityGradientdMicroDeformation,
                       dPlasticGradientMicroVelocityGradientdGradientMicroDeformation;

        if ( evaluateFullDerivatives ){
            dPlasticMacroVelocityGradientdDeformationGradient
                = vectorTools::dot( dPlasticMacroVelocityGradientdElasticRightCauchyGreen,
                                    dElasticRightCauchyGreendDeformationGradient )
                + vectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                    dMacroFlowDirectiondDeformationGradient )
                + vectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondDeformationGradient );
    
            dPlasticMacroVelocityGradientdMicroDeformation
                = vectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                    dMacroFlowDirectiondMicroDeformation )
                + vectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondMicroDeformation );
    
            dPlasticMacroVelocityGradientdGradientMicroDeformation
                = vectorTools::dot( dPlasticMacroVelocityGradientdMacroFlowDirection,
                                    dMacroFlowDirectiondGradientMicroDeformation )
                + vectorTools::dot( dPlasticMacroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondGradientMicroDeformation );
    
            dPlasticMicroVelocityGradientdDeformationGradient
                = vectorTools::dot( dPlasticMicroVelocityGradientdElasticPsi,
                                    dElasticPsidDeformationGradient )
                + vectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondDeformationGradient );
    
            dPlasticMicroVelocityGradientdMicroDeformation
                = vectorTools::dot( dPlasticMicroVelocityGradientdElasticMicroRightCauchyGreen,
                                    dElasticMicroRightCauchyGreendMicroDeformation )
                + vectorTools::dot( dPlasticMicroVelocityGradientdElasticPsi,
                                    dElasticPsidMicroDeformation )
                + vectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondMicroDeformation );
    
            dPlasticMicroVelocityGradientdGradientMicroDeformation
                = vectorTools::dot( dPlasticMicroVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondGradientMicroDeformation );
    
            dPlasticGradientMicroVelocityGradientdDeformationGradient
                = vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticPsi,
                                    dElasticPsidDeformationGradient )
                + vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                    dElasticGammadDeformationGradient )
                + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondDeformationGradient )
                + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                    dMicroGradientFlowDirectiondDeformationGradient );
    
            dPlasticGradientMicroVelocityGradientdMicroDeformation
                = vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticMicroRightCauchyGreen,
                                    dElasticMicroRightCauchyGreendMicroDeformation )
                + vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticPsi,
                                    dElasticPsidMicroDeformation )
                + vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                    dElasticGammadMicroDeformation )
                + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondMicroDeformation )
                + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
                                    dMicroGradientFlowDirectiondMicroDeformation );
    
            dPlasticGradientMicroVelocityGradientdGradientMicroDeformation
                = vectorTools::dot( dPlasticMicroGradientVelocityGradientdElasticGamma,
                                    dElasticGammadGradientMicroDeformation )
                + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroFlowDirection,
                                    dMicroFlowDirectiondGradientMicroDeformation )
                + vectorTools::dot( dPlasticMicroGradientVelocityGradientdMicroGradientFlowDirection,
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
                       vectorTools::appendVectors( dPlasticMicroGradientVelocityGradientdMicroGradientGamma ) );

        //Save the Jacobians of the velocity gradients
        DEBUG.emplace( "dPlasticMacroVelocityGradientdPlasticDeformationGradient",
                       vectorTools::appendVectors( dPlasticMacroVelocityGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dPlasticMacroVelocityGradientdPlasticMicroDeformation",
                       vectorTools::appendVectors( dPlasticMacroVelocityGradientdPlasticMicroDeformation ) );
        DEBUG.emplace( "dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dPlasticMicroVelocityGradientdPlasticDeformationGradient",
                       vectorTools::appendVectors( dPlasticMicroVelocityGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dPlasticMicroVelocityGradientdPlasticMicroDeformation",
                       vectorTools::appendVectors( dPlasticMicroVelocityGradientdPlasticMicroDeformation ) );
        DEBUG.emplace( "dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient",
                       vectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation",
                       vectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation ) );
        DEBUG.emplace( "dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dPlasticMacroVelocityGradientdDeformationGradient",
                            vectorTools::appendVectors( dPlasticMacroVelocityGradientdDeformationGradient ) );
            DEBUG.emplace( "dPlasticMacroVelocityGradientdMicroDeformation",
                            vectorTools::appendVectors( dPlasticMacroVelocityGradientdMicroDeformation ) );
            DEBUG.emplace( "dPlasticMacroVelocityGradientdGradientMicroDeformation",
                            vectorTools::appendVectors( dPlasticMacroVelocityGradientdGradientMicroDeformation ) );
            DEBUG.emplace( "dPlasticMicroVelocityGradientdDeformationGradient",
                            vectorTools::appendVectors( dPlasticMicroVelocityGradientdDeformationGradient ) );
            DEBUG.emplace( "dPlasticMicroVelocityGradientdMicroDeformation",
                            vectorTools::appendVectors( dPlasticMicroVelocityGradientdMicroDeformation ) );
            DEBUG.emplace( "dPlasticMicroVelocityGradientdGradientMicroDeformation",
                            vectorTools::appendVectors( dPlasticMicroVelocityGradientdGradientMicroDeformation ) );
            DEBUG.emplace( "dPlasticMicroGradientVelocityGradientdDeformationGradient",
                            vectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdDeformationGradient ) );
            DEBUG.emplace( "dPlasticMicroGradientVelocityGradientdMicroDeformation",
                            vectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdMicroDeformation ) );
            DEBUG.emplace( "dPlasticMicroGradientVelocityGradientdGradientMicroDeformation",
                            vectorTools::appendVectors( dPlasticGradientMicroVelocityGradientdGradientMicroDeformation ) );
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
            = vectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdMacroGamma );

        variableVector dExpectedPlasticDeformationGradientdMicroGamma
            = vectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdMicroGamma );

        variableVector dExpectedPlasticMicroDeformationdMicroGamma
            = vectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdMicroGamma );

        variableVector dExpectedPlasticGradientMicroDeformationdMacroGamma
            = vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdMacroGamma );

        variableVector dExpectedPlasticGradientMicroDeformationdMicroGamma
            = vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdMicroGamma )
            + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdMicroGamma )
            + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                dPlasticMicroGradientVelocityGradientdMicroGamma );

        variableMatrix dExpectedPlasticGradientMicroDeformationdMicroGradientGamma
            = vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                dPlasticMicroGradientVelocityGradientdMicroGradientGamma );

        //Compute the Jacobians w.r.t. the plastic deformation
        variableMatrix dExpectedPlasticDeformationGradientdPlasticDeformationGradient
            = vectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticDeformationGradient );

        variableMatrix dExpectedPlasticDeformationGradientdPlasticMicroDeformation
            = vectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticMicroDeformation );

        variableMatrix dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation
            = vectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation );

        variableMatrix dExpectedPlasticMicroDeformationdPlasticDeformationGradient
            = vectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticDeformationGradient );

        variableMatrix dExpectedPlasticMicroDeformationdPlasticMicroDeformation
            = vectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticMicroDeformation );

        variableMatrix dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation
            = vectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation );

        variableMatrix dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient
            = vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticDeformationGradient )
            + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticDeformationGradient )
            + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient );

        variableMatrix dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation
            = vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticMicroDeformation )
            + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticMicroDeformation )
            + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation );

        variableMatrix dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation
            = vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation )
            + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation )
            + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
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
                = vectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdDeformationGradient );
    
            dExpectedPlasticDeformationGradientdMicroDeformation
                = vectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdMicroDeformation );
    
            dExpectedPlasticDeformationGradientdGradientMicroDeformation
                = vectorTools::dot( dExpectedPlasticDeformationGradientdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdGradientMicroDeformation );
    
            dExpectedPlasticMicroDeformationdDeformationGradient
                = vectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdDeformationGradient );
    
            dExpectedPlasticMicroDeformationdMicroDeformation
                = vectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdMicroDeformation );
    
            dExpectedPlasticMicroDeformationdGradientMicroDeformation
                = vectorTools::dot( dExpectedPlasticMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdGradientMicroDeformation );
    
            dExpectedPlasticGradientMicroDeformationdDeformationGradient
                = vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdDeformationGradient )
                + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdDeformationGradient )
                + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                    dPlasticGradientMicroVelocityGradientdDeformationGradient );
    
            dExpectedPlasticGradientMicroDeformationdMicroDeformation
                = vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdMicroDeformation )
                + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdMicroDeformation )
                + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
                                    dPlasticGradientMicroVelocityGradientdMicroDeformation );
    
            dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation
                = vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMacroVelocityGradient,
                                    dPlasticMacroVelocityGradientdGradientMicroDeformation )
                + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticMicroVelocityGradient,
                                    dPlasticMicroVelocityGradientdGradientMicroDeformation )
                + vectorTools::dot( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroVelocityGradient,
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
                        vectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdMicroGradientGamma ) );

        //Save the Jacobians of the velocity gradients
        DEBUG.emplace( "dExpectedPlasticDeformationGradientdPlasticDeformationGradient",
                       vectorTools::appendVectors( dExpectedPlasticDeformationGradientdPlasticDeformationGradient ) );
        DEBUG.emplace( "dExpectedPlasticDeformationGradientdPlasticMicroDeformation",
                       vectorTools::appendVectors( dExpectedPlasticDeformationGradientdPlasticMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticMicroDeformationdPlasticDeformationGradient",
                       vectorTools::appendVectors( dExpectedPlasticMicroDeformationdPlasticDeformationGradient ) );
        DEBUG.emplace( "dExpectedPlasticMicroDeformationdPlasticMicroDeformation",
                       vectorTools::appendVectors( dExpectedPlasticMicroDeformationdPlasticMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient",
                       vectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient ) );
        DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation",
                       vectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation ) );
        DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation",
                       vectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dExpectedPlasticDeformationGradientdDeformationGradient",
                            vectorTools::appendVectors( dExpectedPlasticDeformationGradientdDeformationGradient ) );
            DEBUG.emplace( "dExpectedPlasticDeformationGradientdMicroDeformation",
                            vectorTools::appendVectors( dExpectedPlasticDeformationGradientdMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticDeformationGradientdGradientMicroDeformation",
                            vectorTools::appendVectors( dExpectedPlasticDeformationGradientdGradientMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticMicroDeformationdDeformationGradient",
                            vectorTools::appendVectors( dExpectedPlasticMicroDeformationdDeformationGradient ) );
            DEBUG.emplace( "dExpectedPlasticMicroDeformationdMicroDeformation",
                            vectorTools::appendVectors( dExpectedPlasticMicroDeformationdMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticMicroDeformationdGradientMicroDeformation",
                            vectorTools::appendVectors( dExpectedPlasticMicroDeformationdGradientMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdDeformationGradient",
                            vectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdDeformationGradient ) );
            DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdMicroDeformation",
                            vectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdMicroDeformation ) );
            DEBUG.emplace( "dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation",
                            vectorTools::appendVectors( dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation ) );
        }

#endif

        /*!=============================
        | Evaluate the yield equations |
        ==============================*/

        variableVector yieldFunctionValues;
        
        variableType   dMacroYielddMacroCohesion, dMicroYielddMicroCohesion;
        variableVector dMacroYielddPK2Stress, dMacroYielddElasticRightCauchyGreen;
        variableVector dMicroYielddReferenceMicroStress, dMicroYielddElasticRightCauchyGreen;

        variableMatrix dMicroGradientYielddReferenceHigherOrderStress, dMicroGradientYielddMicroGradientCohesion,
                       dMicroGradientYielddElasticRightCauchyGreen;

        error = evaluateYieldFunctions( currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress,
                                        currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion,
                                        currentElasticRightCauchyGreen,
                                       *macroYieldParameters, *microYieldParameters, *microGradientYieldParameters,
                                        yieldFunctionValues,
                                        dMacroYielddPK2Stress, dMacroYielddMacroCohesion, dMacroYielddElasticRightCauchyGreen,
                                        dMicroYielddReferenceMicroStress, dMicroYielddMicroCohesion, dMicroYielddElasticRightCauchyGreen,
                                        dMicroGradientYielddReferenceHigherOrderStress, dMicroGradientYielddMicroGradientCohesion,
                                        dMicroGradientYielddElasticRightCauchyGreen
#ifdef DEBUG_MODE
                                     , DEBUG
#endif
                                    );

        if ( error ){
            errorOut result = new errorNode( "computePlasticDeformationResidual",
                                             "Error in the computation of the yield functions\n" );
            result->addNext( error );
            return result;
        }

        /*!===============================================
        | Construct the Jacobians of the yield equations |
        ================================================*/

        //Construct the Jacobians w.r.t. the plastic deformation measures
        variableVector dMacroYielddPlasticDeformationGradient
            = vectorTools::Tdot( dPK2StressdPlasticDeformationGradient, dMacroYielddPK2Stress )
            + vectorTools::Tdot( dElasticRightCauchyGreendPlasticDeformationGradient, dMacroYielddElasticRightCauchyGreen );

        variableVector dMacroYielddPlasticMicroDeformation
            = vectorTools::Tdot( dPK2StressdPlasticMicroDeformation, dMacroYielddPK2Stress );

        variableVector dMacroYielddPlasticGradientMicroDeformation
            = vectorTools::Tdot( dPK2StressdPlasticGradientMicroDeformation, dMacroYielddPK2Stress );

        variableVector dMicroYielddPlasticDeformationGradient
            = vectorTools::Tdot( dReferenceMicroStressdPlasticDeformationGradient, dMicroYielddReferenceMicroStress )
            + vectorTools::Tdot( dElasticRightCauchyGreendPlasticDeformationGradient, dMicroYielddElasticRightCauchyGreen );

        variableVector dMicroYielddPlasticMicroDeformation
            = vectorTools::Tdot( dReferenceMicroStressdPlasticMicroDeformation, dMicroYielddReferenceMicroStress );

        variableVector dMicroYielddPlasticGradientMicroDeformation
            = vectorTools::Tdot( dReferenceMicroStressdPlasticGradientMicroDeformation, dMicroYielddReferenceMicroStress );

        variableMatrix dMicroGradientYielddPlasticDeformationGradient
            = vectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress, dReferenceHigherOrderStressdPlasticDeformationGradient )
            + vectorTools::dot( dMicroGradientYielddElasticRightCauchyGreen, dElasticRightCauchyGreendPlasticDeformationGradient );

        variableMatrix dMicroGradientYielddPlasticMicroDeformation
            = vectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress, dReferenceHigherOrderStressdPlasticMicroDeformation );

        variableMatrix dMicroGradientYielddPlasticGradientMicroDeformation
            = vectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress,
                                dReferenceHigherOrderStressdPlasticGradientMicroDeformation );

        //Construct the Jacobians w.r.t. the strain-like ISV
        variableType dMacroYielddMacroStrainISV = dMacroYielddMacroCohesion * dMacroCohesiondMacroStrainISV;
        variableType dMicroYielddMicroStrainISV = dMicroYielddMicroCohesion * dMicroCohesiondMicroStrainISV;
        variableMatrix dMicroGradientYielddMicroGradientStrainISV = vectorTools::dot( dMicroGradientYielddMicroGradientCohesion,
                                                                                      dMicroGradientCohesiondMicroGradientStrainISV );

        variableVector dMacroYielddDeformationGradient, dMacroYielddMicroDeformation, dMacroYielddGradientMicroDeformation,
                       dMicroYielddDeformationGradient, dMicroYielddMicroDeformation, dMicroYielddGradientMicroDeformation;
        variableMatrix dMicroGradientYielddDeformationGradient, dMicroGradientYielddMicroDeformation,
                       dMicroGradientYielddGradientMicroDeformation;

        if ( evaluateFullDerivatives ){
            dMacroYielddDeformationGradient
                = vectorTools::Tdot( dPK2StressdDeformationGradient, dMacroYielddPK2Stress )
                + vectorTools::Tdot( dElasticRightCauchyGreendDeformationGradient, dMacroYielddElasticRightCauchyGreen );
    
            dMacroYielddMicroDeformation
                = vectorTools::Tdot( dPK2StressdMicroDeformation, dMacroYielddPK2Stress );
    
            dMacroYielddGradientMicroDeformation
                = vectorTools::Tdot( dPK2StressdGradientMicroDeformation, dMacroYielddPK2Stress );
    
            dMicroYielddDeformationGradient
                = vectorTools::Tdot( dReferenceMicroStressdDeformationGradient, dMicroYielddReferenceMicroStress )
                + vectorTools::Tdot( dElasticRightCauchyGreendDeformationGradient, dMicroYielddElasticRightCauchyGreen );
    
            dMicroYielddMicroDeformation
                = vectorTools::Tdot( dReferenceMicroStressdMicroDeformation, dMicroYielddReferenceMicroStress );
    
            dMicroYielddGradientMicroDeformation
                = vectorTools::Tdot( dReferenceMicroStressdGradientMicroDeformation, dMicroYielddReferenceMicroStress );
    
            dMicroGradientYielddDeformationGradient
                = vectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress, dReferenceHigherOrderStressdDeformationGradient )
                + vectorTools::dot( dMicroGradientYielddElasticRightCauchyGreen, dElasticRightCauchyGreendDeformationGradient );
    
            dMicroGradientYielddMicroDeformation
                = vectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress, dReferenceHigherOrderStressdMicroDeformation );
    
            dMicroGradientYielddGradientMicroDeformation
                = vectorTools::dot( dMicroGradientYielddReferenceHigherOrderStress,
                                    dReferenceHigherOrderStressdGradientMicroDeformation );
        }

#ifdef DEBUG_MODE

        //Save the values to the debug map
        temp = { yieldFunctionValues[ 0 ] };
        DEBUG.emplace( "macroYieldFunction", temp );
 
        temp = { yieldFunctionValues[ 1 ] };
        DEBUG.emplace( "microYieldFunction", temp );
 
        DEBUG.emplace( "microGradientYieldFunction",
                        variableVector( yieldFunctionValues.begin() + 2, yieldFunctionValues.begin() + 5 ) );
 
        //Save the jacobians w.r.t. the plastic deformation
        DEBUG.emplace( "dMacroYielddPlasticDeformationGradient", dMacroYielddPlasticDeformationGradient );
        DEBUG.emplace( "dMacroYielddPlasticMicroDeformation", dMacroYielddPlasticMicroDeformation );
        DEBUG.emplace( "dMacroYielddPlasticGradientMicroDeformation", dMacroYielddPlasticGradientMicroDeformation );
 
        DEBUG.emplace( "dMicroYielddPlasticDeformationGradient", dMicroYielddPlasticDeformationGradient );
        DEBUG.emplace( "dMicroYielddPlasticMicroDeformation", dMicroYielddPlasticMicroDeformation );
        DEBUG.emplace( "dMicroYielddPlasticGradientMicroDeformation", dMicroYielddPlasticGradientMicroDeformation );
 
        DEBUG.emplace( "dMicroGradientYielddPlasticDeformationGradient",
                        vectorTools::appendVectors( dMicroGradientYielddPlasticDeformationGradient ) );
        DEBUG.emplace( "dMicroGradientYielddPlasticMicroDeformation",
                        vectorTools::appendVectors( dMicroGradientYielddPlasticMicroDeformation ) );
        DEBUG.emplace( "dMicroGradientYielddPlasticGradientMicroDeformation",
                        vectorTools::appendVectors( dMicroGradientYielddPlasticGradientMicroDeformation ) );
 
        //Save the Jacobians w.r.t. the strain-like ISVs
        temp = { dMacroYielddMacroStrainISV };
        DEBUG.emplace( "dMacroYielddMacroStrainISV", temp );
        temp = { dMicroYielddMicroStrainISV };
        DEBUG.emplace( "dMicroYielddMicroStrainISV", temp );
        DEBUG.emplace( "dMicroGradientYielddMicroGradientStrainISV",
                        vectorTools::appendVectors( dMicroGradientYielddMicroGradientStrainISV ) );

        if ( evaluateFullDerivatives ){
            DEBUG.emplace( "dMacroYielddDeformationGradient", dMacroYielddDeformationGradient );
            DEBUG.emplace( "dMacroYielddMicroDeformation", dMacroYielddMicroDeformation );
            DEBUG.emplace( "dMacroYielddGradientMicroDeformation", dMacroYielddGradientMicroDeformation );
            DEBUG.emplace( "dMicroYielddDeformationGradient", dMicroYielddDeformationGradient );
            DEBUG.emplace( "dMicroYielddMicroDeformation", dMicroYielddMicroDeformation );
            DEBUG.emplace( "dMicroYielddGradientMicroDeformation", dMicroYielddGradientMicroDeformation );
            DEBUG.emplace( "dMicroGradientYielddDeformationGradient",
                            vectorTools::appendVectors( dMicroGradientYielddDeformationGradient ) );
            DEBUG.emplace( "dMicroGradientYielddMicroDeformation",
                            vectorTools::appendVectors( dMicroGradientYielddMicroDeformation ) );
            DEBUG.emplace( "dMicroGradientYielddGradientMicroDeformation",
                            vectorTools::appendVectors( dMicroGradientYielddGradientMicroDeformation ) );
        }

#endif

        
        /*!===============================================
        | Compute the residual equation and the Jacobian |
        ================================================*/

        residual = solverTools::floatVector( x.size(), 0 );
        jacobian = vectorTools::eye< solverTools::floatType >( x.size() );

        if ( evaluateFullDerivatives ){
            floatOuts[ 5 ] = solverTools::floatVector( x.size() * 45, 0 );
        }

        //Compute the residuals and the jacobians for the plastic deformations
        for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
            residual[ i ] = currentPlasticDeformationGradient[ i ] - expectedPlasticDeformationGradient[ i ];

            //The plastic deformation residuals
            for ( unsigned int j = 0; j < currentPlasticDeformationGradient.size(); j++ ){
                jacobian[ i ][ j ] -= dExpectedPlasticDeformationGradientdPlasticDeformationGradient[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ 45 * i + j ] -= dExpectedPlasticDeformationGradientdDeformationGradient[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticMicroDeformation.size(); j++ ){
                jacobian[ i ][ j + 9 ] -= dExpectedPlasticDeformationGradientdPlasticMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ 45 * i + j + 9 ] -= dExpectedPlasticDeformationGradientdMicroDeformation[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticGradientMicroDeformation.size(); j++ ){
                jacobian[ i ][ j + 18 ] -= dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ 45 * i + j + 18 ] -= dExpectedPlasticDeformationGradientdGradientMicroDeformation[ i ][ j ];
                }
            }

            //The strain-like ISV residuals
            
            //The plastic multiplier residuals
            jacobian[ i ][ 50 ] -= dExpectedPlasticDeformationGradientdMacroGamma[ i ];
            jacobian[ i ][ 51 ] -= dExpectedPlasticDeformationGradientdMicroGamma[ i ];

        }

        for ( unsigned int i = 0; i < currentPlasticMicroDeformation.size(); i++ ){
            residual[ i + 9 ] = currentPlasticMicroDeformation[ i ] - expectedPlasticMicroDeformation[ i ];

            for ( unsigned int j = 0; j < currentPlasticDeformationGradient.size(); j++ ){
                jacobian[ i + 9 ][ j ] -= dExpectedPlasticMicroDeformationdPlasticDeformationGradient[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ 45 * ( i + 9 ) + j ] -= dExpectedPlasticMicroDeformationdDeformationGradient[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticMicroDeformation.size(); j++ ){
                jacobian[ i + 9 ][ j + 9 ] -= dExpectedPlasticMicroDeformationdPlasticMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ 45 * ( i + 9 ) + j + 9 ] -= dExpectedPlasticMicroDeformationdMicroDeformation[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticGradientMicroDeformation.size(); j++ ){
                jacobian[ i + 9 ][ j + 18 ] -= dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ 45 * ( i + 9 ) + j + 18 ] -= dExpectedPlasticMicroDeformationdGradientMicroDeformation[ i ][ j ];
                }
            }

            //The strain-like ISV residuals
            
            //The plastic multiplier residuals
            jacobian[ i + 9 ][ 51 ] -= dExpectedPlasticMicroDeformationdMicroGamma[ i ];
        }

        for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
            residual[ i + 18 ] = currentPlasticGradientMicroDeformation[ i ] - expectedPlasticGradientMicroDeformation[ i ];

            for ( unsigned int j = 0; j < currentPlasticDeformationGradient.size(); j++ ){
                jacobian[ i + 18 ][ j ] -= dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ 45 * ( i + 18 ) + j ] -= dExpectedPlasticGradientMicroDeformationdDeformationGradient[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticMicroDeformation.size(); j++ ){
                jacobian[ i + 18 ][ j + 9 ] -= dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ 45 * ( i + 18 ) + j + 9 ] -= dExpectedPlasticGradientMicroDeformationdMicroDeformation[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticGradientMicroDeformation.size(); j++ ){
                jacobian[ i + 18 ][ j + 18 ] -= dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ 45 * ( i + 18 ) + j + 18 ] -= dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation[ i ][ j ];
                }
            }

            //The strain-like ISV residuals are all zero

            //The plastic multiplier residuals
            jacobian[ i + 18 ][ 50 ] -= dExpectedPlasticGradientMicroDeformationdMacroGamma[ i ]; 
            jacobian[ i + 18 ][ 51 ] -= dExpectedPlasticGradientMicroDeformationdMicroGamma[ i ]; 
            for ( unsigned int j = 0; j < currentMicroGradientGamma.size(); j++ ){
                jacobian[ i + 18 ][ 52 + j ] -= dExpectedPlasticGradientMicroDeformationdMicroGradientGamma[ i ][ j ];
            }
        }

        //Compute the residuals and the Jacobians for the plastic strain-like ISVs
        residual[ 45 ] = currentMacroStrainISV - expectedMacroStrainISV;
        jacobian[ 45 ][ 50 ] = -dExpectedMacroISVdMacroGamma;

        residual[ 46 ] = currentMicroStrainISV - expectedMicroStrainISV;
        jacobian[ 46 ][ 51 ] = -dExpectedMicroISVdMicroGamma;

        for ( unsigned int i = 0; i < currentMicroGradientStrainISV.size(); i++ ){
            residual[ 47 + i ] = currentMicroGradientStrainISV[ i ] - expectedMicroGradientStrainISV[ i ];
            
            for ( unsigned int j = 0; j < currentMicroGradientGamma.size(); j++ ){
                jacobian[ 47 + i ][ 52 + j ] = -dExpectedMicroGradientISVdMicroGradientGamma[ i ][ j ];
            }
        }

        //Compute the residuals and the Jacobians for the Kuhn-Tucker and yield conditions

        //Set the residuals and Jacobians for the macro-scale terms
        variableType dMacMacroYielddMacroYield;
        residual[ 50 ] = currentMacroGamma * yieldFunctionValues[ 0 ]
                       + constitutiveTools::mac( yieldFunctionValues[ 0 ], dMacMacroYielddMacroYield );

        //The Jacobian terms w.r.t. the plastic deformation
        for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
            jacobian[ 50 ][ i ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddPlasticDeformationGradient[ i ];

            if ( evaluateFullDerivatives ){
                floatOuts[ 5 ][ 50 * 45 + i ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddDeformationGradient[ i ];
            }
        }

        for ( unsigned int i = 0; i < currentPlasticMicroDeformation.size(); i++ ){
            jacobian[ 50 ][ i + 9 ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddPlasticMicroDeformation[ i ];

            if ( evaluateFullDerivatives ){
                floatOuts[ 5 ][ 50 * 45 + i + 9 ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddMicroDeformation[ i ];
            }
        }

        for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
            jacobian[ 50 ][ i + 18 ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddPlasticGradientMicroDeformation[ i ];

            if ( evaluateFullDerivatives ){
                floatOuts[ 5 ][ 50 * 45 + i + 18 ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddGradientMicroDeformation[ i ];
            }
        }

        //The Jacobian terms w.r.t. the strain-like ISVs
        jacobian[ 50 ][ 45 ] = ( currentMacroGamma + dMacMacroYielddMacroYield ) * dMacroYielddMacroStrainISV;

        //The Jacobian terms w.r.t. the plastic multipliers
        jacobian[ 50 ][ 50 ] = yieldFunctionValues[ 0 ];

        //Set the residuals and Jacobians for the micro-scale terms
        variableType dMacMicroYielddMicroYield;
        residual[ 51 ] = currentMicroGamma * yieldFunctionValues[ 1 ]
                       + constitutiveTools::mac( yieldFunctionValues[ 1 ], dMacMicroYielddMicroYield );

        //The Jacobian terms w.r.t. the plastic deformation
        for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
            jacobian[ 51 ][ i ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddPlasticDeformationGradient[ i ];

            if ( evaluateFullDerivatives ){
                floatOuts[ 5 ][ 51 * 45 + i ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddDeformationGradient[ i ];
            }
        }

        for ( unsigned int i = 0; i < currentPlasticMicroDeformation.size(); i++ ){
            jacobian[ 51 ][ i + 9 ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddPlasticMicroDeformation[ i ];

            if ( evaluateFullDerivatives ){
                floatOuts[ 5 ][ 51 * 45 + i + 9 ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddMicroDeformation[ i ];
            }
        }

        for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
            jacobian[ 51 ][ i + 18 ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddPlasticGradientMicroDeformation[ i ];

            if ( evaluateFullDerivatives ){
                floatOuts[ 5 ][ 51 * 45 + i + 18 ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddGradientMicroDeformation[ i ];
            }
        }

        //The Jacobian terms w.r.t. the strain-like ISVs
        jacobian[ 51 ][ 46 ] = ( currentMicroGamma + dMacMicroYielddMicroYield ) * dMicroYielddMicroStrainISV;

        //The Jacobian terms w.r.t. the plastic multipliers
        jacobian[ 51 ][ 51 ] = yieldFunctionValues[ 1 ];

        //Set the residuals and jacobians for the micro gradient terms
        variableType dMacMicroGradientYielddMicroGradientYield;
        for ( unsigned int i = 0; i < 3; i++ ){
            residual[ 52 + i ] = currentMicroGradientGamma[ i ] * yieldFunctionValues[ 2 + i ]
                               + constitutiveTools::mac( yieldFunctionValues[ 2 + i ], dMacMicroGradientYielddMicroGradientYield );

            //The Jacobian terms w.r.t. the plastic deformation
            for ( unsigned int j = 0; j < currentPlasticDeformationGradient.size(); j++ ){
                jacobian[ 52 + i ][ j ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddPlasticDeformationGradient[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ ( 52 + i ) * 45 + j ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddDeformationGradient[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticMicroDeformation.size(); j++ ){
                jacobian[ 52 + i ][ j + 9 ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddPlasticMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ ( 52 + i ) * 45 + j + 9 ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddMicroDeformation[ i ][ j ];
                }
            }

            for ( unsigned int j = 0; j < currentPlasticGradientMicroDeformation.size(); j++ ){
                jacobian[ 52 + i ][ j + 18 ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddPlasticGradientMicroDeformation[ i ][ j ];

                if ( evaluateFullDerivatives ){
                    floatOuts[ 5 ][ ( 52 + i ) * 45 + j + 18 ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddGradientMicroDeformation[ i ][ j ];
                }
            }

            //The Jacobian terms w.r.t. the strain-like ISVs
            for ( unsigned int j = 0; j < currentMicroGradientStrainISV.size(); j++ ){
                jacobian[ 52 + i ][ 47 + j ] = ( currentMicroGradientGamma[ i ] + dMacMicroGradientYielddMicroGradientYield ) * dMicroGradientYielddMicroGradientStrainISV[ i ][ j ];
            }

            //The Jacobian terms w.r.t. the plastic multipliers
            jacobian[ 52 + i ][ 52 + i ] = yieldFunctionValues[ 2 + i ];
        }

        //Save the stresses
        floatOuts[ 0 ] = currentPK2Stress;
        floatOuts[ 1 ] = currentReferenceMicroStress;
        floatOuts[ 2 ] = currentReferenceHigherOrderStress;

        //Save whether the function is yielding
        intOuts[ 0 ] = isYielding;

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
                        , solverTools::homotopyMap &DEBUG
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
         * :param solverTools::homotopyMap DEBUG: The debugging object ( only available if DEBUG_MODE is defined )
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        //Assume 3D
        unsigned int dim = 3;

#ifdef DEBUG_MODE
        solverTools::debugMap tempDEBUG;
        solverTools::iterationMap tempITERATION;
#endif

        //Construct identity matrix
        constantVector eye( dim * dim, 0 );
        vectorTools::eye< constantType >( eye );

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
        variableType previousdMacroGdMacroCohesion, previousdMicroGdMicroCohesion;
        variableMatrix previousdMicroGradientGdMicroGradientCohesion( previousMicroGradientStrainISV.size(),
                                                                      variableVector( previousMicroGradientStrainISV.size() ) );

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
        solverTools::floatVector tmp = { previousMacroCohesion };
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

//            std::cout << "  previous yielding detected\n";
            //Compute the previous stress
            error = micromorphicLinearElasticity::linearElasticityReference( previousElasticDeformationGradient,
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
                               vectorTools::appendVectors( previousdMicroGradientGdMicroGradientCohesion ) );
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

        //Compute the right Cauchy-Green deformation gradient
        variableVector currentElasticRightCauchyGreen;

        error = constitutiveTools::computeRightCauchyGreen( currentElasticDeformationGradient, currentElasticRightCauchyGreen );

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

        error = micromorphicLinearElasticity::linearElasticityReference( currentElasticDeformationGradient,
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

#ifdef DEBUG_MODE
        tempDEBUG.emplace( "intermediatePK2Stress", currentPK2Stress );
        tempDEBUG.emplace( "intermediateReferenceMicroStress", currentReferenceMicroStress );
        tempDEBUG.emplace( "intermediateReferenceHigherOrderStress", currentReferenceHigherOrderStress );
#endif

//        std::cout << "current stress measures\n";
//        std::cout << "currentPK2Stress:\n"; vectorTools::print( currentPK2Stress );
//        std::cout << "currentReferenceMicroStress:\n"; vectorTools::print( currentReferenceMicroStress );
//        std::cout << "currentReferenceHigherOrderStress:\n"; vectorTools::print( currentReferenceHigherOrderStress );

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

//        std::cout << "initial yield function values:\n"; vectorTools::print( currentYieldFunctionValues );

        /*============================
        | Begin the non-linear solve |
        ============================*/
//        std::cout << "beginning the nonlinear solve\n";

        //Check if any of the surfaces are yielding and begin the non-linear solver if they are
        solverTools::floatVector solutionVector;
        solverTools::intVector activePlasticity( currentYieldFunctionValues.size(), 0 );
        bool convergenceFlag = true;

        for ( unsigned int i = 0; i < currentYieldFunctionValues[ i ]; i++ ){
            if ( currentYieldFunctionValues[ i ] > absoluteTolerance ){

                convergenceFlag = false;
                
                activePlasticity[ i ] = 1;
            }
        }
//        std::cout << "convergenceFlag: " << convergenceFlag << "\n";
//        std::cout << "activePlasticity: "; vectorTools::print( activePlasticity );

        if ( !convergenceFlag ){
            solverTools::floatMatrix floatArgs =
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
                    vectorTools::appendVectors( previousdMicroGradientGdMicroGradientCohesion ),
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
        
            solverTools::floatMatrix floatOuts =
                {
                    currentPK2Stress,
                    currentReferenceMicroStress,
                    currentReferenceHigherOrderStress,
                };

            solverTools::intMatrix intOuts = { activePlasticity };

            solverTools::stdFncNLFJ func
                = static_cast< solverTools::NonLinearFunctionWithJacobian >( computePlasticDeformationResidual );

            solverTools::floatVector x0 = { currentMacroStrainISV, currentMicroStrainISV,
                                            currentMicroGradientStrainISV[ 0 ],
                                            currentMicroGradientStrainISV[ 1 ],
                                            currentMicroGradientStrainISV[ 2 ],
                                            currentMacroGamma, currentMicroGamma, 
                                            currentMicroGradientGamma[ 0 ],
                                            currentMicroGradientGamma[ 1 ],
                                            currentMicroGradientGamma[ 2 ]
                                          };

            x0 = vectorTools::appendVectors( { currentPlasticDeformationGradient,
                                               currentPlasticMicroDeformation,
                                               currentPlasticGradientMicroDeformation, x0 } );

            solverTools::floatVector solutionVector;

            bool convergeFlag, fatalErrorFlag;
           
            solverTools::intMatrix intArgs = { { 0 } };
//            std::cout << "entering homotopy solver\n"; 
//            error = solverTools::newtonRaphson( func, x0, solutionVector, convergeFlag, fatalErrorFlag,
//                                                floatOuts, intOuts, floatArgs, intArgs,
//                                                20, relativeTolerance, absoluteTolerance );

            error = solverTools::homotopySolver( func, x0, solutionVector, convergeFlag, fatalErrorFlag,
                                                 floatOuts, intOuts, floatArgs, intArgs,
#ifdef DEBUG_MODE
                                                 DEBUG,
#endif
                                                 20, relativeTolerance, absoluteTolerance,
                                                 1e-4, 5, 1.0, .1 );

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
                return 1;
            }

            //Extract the deformation measures
//            std::cout << "solutionVector:\n"; vectorTools::print( solutionVector );
            currentPlasticDeformationGradient      = variableVector( solutionVector.begin() +  0, solutionVector.begin() +  9 );
            currentPlasticMicroDeformation         = variableVector( solutionVector.begin() +  9, solutionVector.begin() + 18 );
            currentPlasticGradientMicroDeformation = variableVector( solutionVector.begin() + 18, solutionVector.begin() + 45 );
            currentMacroStrainISV                  = solutionVector[ 45 ];
            currentMicroStrainISV                  = solutionVector[ 46 ];
            currentMicroGradientStrainISV          = variableVector( solutionVector.begin() + 47, solutionVector.begin() + 50 );
            currentMacroGamma                      = solutionVector[ 50 ];
            currentMicroGamma                      = solutionVector[ 51 ];
            currentMicroGradientGamma              = variableVector( solutionVector.begin() + 52, solutionVector.begin() + 55 );

//            std::cout << "currentPlasticDeformationGradient:\n"; vectorTools::print( currentPlasticDeformationGradient );
//            std::cout << "currentPlasticMicroDeformation:\n"; vectorTools::print( currentPlasticMicroDeformation );
//            std::cout << "currentPlasticGradientMicroDeformation:\n"; vectorTools::print( currentPlasticGradientMicroDeformation );
//            std::cout << "currentMacroStrainISV:\n" << currentMacroStrainISV << "\n";
//            std::cout << "currentMicroStrainISV:\n" << currentMicroStrainISV << "\n";
//            std::cout << "currentMicroGradientStrainISV:\n"; vectorTools::print( currentMicroGradientStrainISV );
//            std::cout << "currentMacroGamma:\n" << currentMacroGamma << "\n";
//            std::cout << "currentMicroStrainISV:\n" << currentMicroGamma << "\n";
//            std::cout << "currentMicroGradientStrainISV:\n"; vectorTools::print( currentMicroGradientGamma );

            //Extract the stresses
            currentPK2Stress                       = floatOuts[ 0 ];
            currentReferenceMicroStress            = floatOuts[ 1 ];
            currentReferenceHigherOrderStress      = floatOuts[ 2 ];

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "convergedPlasticDeformationGradient", currentPlasticDeformationGradient );
            tempDEBUG.emplace( "convergedPlasticMicroDeformation", currentPlasticMicroDeformation );
            tempDEBUG.emplace( "convergedPlasticGradientMicroDeformation", currentPlasticGradientMicroDeformation );

            solverTools::floatVector temp = { currentMacroStrainISV };
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
        currentStrainISVS = vectorTools::appendVectors( { currentStrainISVS, currentMicroGradientGamma } );

        std::vector< double > currentGammas = { currentMacroGamma, currentMicroGamma };
        currentGammas = vectorTools::appendVectors( { currentGammas, currentMicroGradientGamma } );

        SDVS = vectorTools::appendVectors( { currentStrainISVS, currentGammas, currentPlasticDeformationGradient - eye,
                                             currentPlasticMicroDeformation - eye, currentPlasticGradientMicroDeformation } );

        //Output the stresses mapped to the true reference configuration
        error = micromorphicTools::pullBackCauchyStress( currentPK2Stress, currentPlasticDeformationGradient, current_PK2 );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in pullback operation on the PK2 stress" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        error = micromorphicTools::pullBackMicroStress( currentReferenceMicroStress, currentPlasticDeformationGradient, current_SIGMA );

        if ( error ){
            errorOut result = new errorNode( "evaluate_model",
                                             "Error in pullback operation on the reference symmetric micro-stress" );
            result->addNext( error );
            result->print();           //Print the error message
            output_message = buffer.str(); //Save the output to enable message passing
            return 2;
        }

        error = micromorphicTools::pullBackHigherOrderStress( currentReferenceHigherOrderStress, 
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

//        std::cout << "SDVS:\n"; vectorTools::print( SDVS );
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
                        , solverTools::homotopyMap &DEBUG
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
         *     terms w.r.t. the deformation ( unused )
         * :param std::string &output_message: The output message string.
         * :param solverTools::homotopyMap DEBUG: The debugging map ( only available if DEBUG_MODE is defined )
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        //Assume 3D
        unsigned int dim = 3;

#ifdef DEBUG_MODE
        solverTools::debugMap tempDEBUG;
        solverTools::iterationMap tempITERATION;
#endif

        //Construct identity matrix
        constantVector eye( dim * dim, 0 );
        vectorTools::eye< constantType >( eye );

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
                           vectorTools::appendVectors( dDeformationGradientdGradientMacroDisplacement ) );
        tempDEBUG.emplace( "totaldMicroDeformationdMicroDisplacement",
                           vectorTools::appendVectors( dMicroDeformationdMicroDisplacement ) );
        tempDEBUG.emplace( "totaldGradientMicroDeformationdGradientMicroDisplacement",
                           vectorTools::appendVectors( dGradientMicroDeformationdGradientMicroDisplacement ) );
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
        variableType previousdMacroGdMacroCohesion, previousdMicroGdMicroCohesion;
        variableMatrix previousdMicroGradientGdMicroGradientCohesion( previousMicroGradientStrainISV.size(),
                                                                      variableVector( previousMicroGradientStrainISV.size() ) );

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
        solverTools::floatVector tmp = { previousMacroCohesion };
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
            error = micromorphicLinearElasticity::linearElasticityReference( previousElasticDeformationGradient,
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
                               vectorTools::appendVectors( previousdMicroGradientGdMicroGradientCohesion ) );
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

        error = constitutiveTools::computeRightCauchyGreen( currentElasticDeformationGradient, currentElasticRightCauchyGreen );

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

        error = micromorphicLinearElasticity::linearElasticityReference( currentElasticDeformationGradient,
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
        solverTools::floatVector solutionVector;
        solverTools::intVector activePlasticity( currentYieldFunctionValues.size(), 0 );
        bool convergenceFlag = true;

        for ( unsigned int i = 0; i < currentYieldFunctionValues[ i ]; i++ ){
            if ( currentYieldFunctionValues[ i ] > absoluteTolerance ){

                convergenceFlag = false;
                
                activePlasticity[ i ] = 1;
            }
        }

        if ( !convergenceFlag ){
            solverTools::floatMatrix floatArgs =
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
                    vectorTools::appendVectors( previousdMicroGradientGdMicroGradientCohesion ),
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
        
            variableVector dStressdSolutionVector( 45 * 55, 0 );
            variableVector dStressdDeformation( 45 * 45, 0 );
            variableVector dResidualdDeformation( 55 * 45, 0 );
            solverTools::floatMatrix floatOuts =
                {
                    currentPK2Stress,
                    currentReferenceMicroStress,
                    currentReferenceHigherOrderStress,
                    dStressdSolutionVector,
                    dStressdDeformation,
                    dResidualdDeformation
                };

            solverTools::intMatrix intOuts = { activePlasticity };

            solverTools::stdFncNLFJ func
                = static_cast< solverTools::NonLinearFunctionWithJacobian >( computePlasticDeformationResidual );

            solverTools::floatVector x0 = { currentMacroStrainISV, currentMicroStrainISV,
                                            currentMicroGradientStrainISV[ 0 ],
                                            currentMicroGradientStrainISV[ 1 ],
                                            currentMicroGradientStrainISV[ 2 ],
                                            currentMacroGamma, currentMicroGamma, 
                                            currentMicroGradientGamma[ 0 ],
                                            currentMicroGradientGamma[ 1 ],
                                            currentMicroGradientGamma[ 2 ]
                                          };

            x0 = vectorTools::appendVectors( { currentPlasticDeformationGradient,
                                               currentPlasticMicroDeformation,
                                               currentPlasticGradientMicroDeformation, x0 } );

            solverTools::floatVector solutionVector;

            bool convergeFlag, fatalErrorFlag;
           
            solverTools::intMatrix intArgs = { { 1 } };
            //TODO: I would like to compute the derivative of the residual w.r.t. the fundamental deformation measures
            //      more efficiently in the future rather than every time.

            solverTools::solverType finalSolver;
            solverTools::floatMatrix J;
            error = solverTools::homotopySolver( func, x0, solutionVector, convergeFlag, fatalErrorFlag,
                                                 floatOuts, intOuts, floatArgs, intArgs, finalSolver, J,
#ifdef DEBUG_MODE
                                                 DEBUG,
#endif
                                                 20, relativeTolerance, absoluteTolerance,
                                                 1e-4, 5, 1.0, .1 );

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
                return 1;
            }

            for ( unsigned int i = 0; i < 5; i++ ){
                output_message += std::to_string( solutionVector[ 50 + i ] ) + " ";
            }

            //Extract the deformation measures
            currentPlasticDeformationGradient      = variableVector( solutionVector.begin() +  0, solutionVector.begin() +  9 );
            currentPlasticMicroDeformation         = variableVector( solutionVector.begin() +  9, solutionVector.begin() + 18 );
            currentPlasticGradientMicroDeformation = variableVector( solutionVector.begin() + 18, solutionVector.begin() + 45 );
            currentMacroStrainISV                  = solutionVector[ 45 ];
            currentMicroStrainISV                  = solutionVector[ 46 ];
            currentMicroGradientStrainISV          = variableVector( solutionVector.begin() + 47, solutionVector.begin() + 50 );
            currentMacroGamma                      = solutionVector[ 50 ];
            currentMicroGamma                      = solutionVector[ 51 ];
            currentMicroGradientGamma              = variableVector( solutionVector.begin() + 52, solutionVector.begin() + 55 );

            //Extract the stresses and Jacobians
            currentPK2Stress                       = floatOuts[ 0 ];
            currentReferenceMicroStress            = floatOuts[ 1 ];
            currentReferenceHigherOrderStress      = floatOuts[ 2 ];
            dStressdSolutionVector                 = floatOuts[ 3 ];
            dStressdDeformation                    = floatOuts[ 4 ];
            dResidualdDeformation                  = floatOuts[ 5 ];

#ifdef DEBUG_MODE
            tempDEBUG.emplace( "convergedPlasticDeformationGradient", currentPlasticDeformationGradient );
            tempDEBUG.emplace( "convergedPlasticMicroDeformation", currentPlasticMicroDeformation );
            tempDEBUG.emplace( "convergedPlasticGradientMicroDeformation", currentPlasticGradientMicroDeformation );

            solverTools::floatVector temp = { currentMacroStrainISV };
            tempDEBUG.emplace( "convergedMacroStrainISV", temp );

            temp = { currentMicroStrainISV };
            tempDEBUG.emplace( "convergedMicroStrainISV", temp );

            tempDEBUG.emplace( "convergedMicroGradientStrainISV", currentMicroGradientStrainISV );

            temp = { currentMacroGamma };
            tempDEBUG.emplace( "convergedMacroGamma", temp );

            temp = { currentMicroGamma };
            tempDEBUG.emplace( "convergedMicroGamma", temp );

            tempDEBUG.emplace( "convergedMicroGradientGamma", currentMicroGradientGamma );

            tempDEBUG.emplace( "intermediatedPK2Stress", currentPK2Stress );
            tempDEBUG.emplace( "intermediateReferenceMicroStress", currentReferenceMicroStress );
            tempDEBUG.emplace( "intermediateReferenceHigherOrderStress", currentReferenceHigherOrderStress );
#endif

            //Form the Jacobian of the stresses w.r.t. the fundamental deformation measures

            //Compute decomposition of the final Jacobian
            solverTools::floatVector dResidualdSolutionVector = vectorTools::appendVectors( J );
            Eigen::Map< const Eigen::Matrix< solverTools::floatType, -1, -1, Eigen::RowMajor > >
                dRdX( dResidualdSolutionVector.data(), 55, 55 );

            solverTools::solverType linearSolver( dRdX );

            Eigen::Map< const Eigen::Matrix< solverTools::floatType, -1, -1, Eigen::RowMajor > >
                dRdD( dResidualdDeformation.data(), 55, 45 );

            Eigen::Map< const Eigen::Matrix< solverTools::floatType, -1, -1, Eigen::RowMajor > >
                dSdX( dStressdSolutionVector.data(), 45, 55 );

            Eigen::Map< const Eigen::Matrix< solverTools::floatType, -1, -1, Eigen::RowMajor > >
                dSdD( dStressdDeformation.data(), 45, 45 );

            //Solve for the total derivative of the residual w.r.t. 
            solverTools::floatVector DSolutionVectorDDeformation( 55 * 45, 0 );
            Eigen::Map< Eigen::Matrix< solverTools::floatType, -1, -1, Eigen::RowMajor > >
                DXDD( DSolutionVectorDDeformation.data(), 55, 45 );

            //Perform the linear solve
            DXDD = -linearSolver.solve( dRdD );

            solverTools::floatVector DStressDDeformation( 45 * 45, 0 );
            Eigen::Map< Eigen::Matrix< solverTools::floatType, -1, -1, Eigen::RowMajor > >
                DSDD( DStressDDeformation.data(), 45, 45 );

            DSDD = dSdD + dSdX * DXDD;

            //Extract the total derivatives of the plastic deformation measures and stresses
            for ( unsigned int i = 0; i < 9; i++ ){
                for ( unsigned int j = 0; j < 9; j++ ){

                    //The plastic deformation measures
                    dPlasticDeformationGradientdDeformationGradient[ i ][ j ] = DSolutionVectorDDeformation[ 45 * i + j ];
                    dPlasticDeformationGradientdMicroDeformation[ i ][ j ]    = DSolutionVectorDDeformation[ 45 * i + j + 9 ];

                    dPlasticMicroDeformationdDeformationGradient[ i ][ j ]    = DSolutionVectorDDeformation[ 45 * ( i + 9 ) + j ];
                    dPlasticMicroDeformationdMicroDeformation[ i ][ j ]       = DSolutionVectorDDeformation[ 45 * ( i + 9 ) + j + 9 ];

                    //The stress measures
                    dPK2StressdDeformationGradient[ i ][ j ]            = DStressDDeformation[ 45 * i + j ];
                    dPK2StressdMicroDeformation[ i ][ j ]               = DStressDDeformation[ 45 * i + j + 9 ];

                    dReferenceMicroStressdDeformationGradient[ i ][ j ] = DStressDDeformation[ 45 * ( i + 9 ) + j ];
                    dReferenceMicroStressdMicroDeformation[ i ][ j ]    = DStressDDeformation[ 45 * ( i + 9 ) + j + 9 ];
                }

                for ( unsigned int j = 0; j < 27; j++ ){

                    //The plastic deformation measures
                    dPlasticDeformationGradientdGradientMicroDeformation[ i ][ j ]
                        = DSolutionVectorDDeformation[ 45 * i + j + 18 ];
                    dPlasticMicroDeformationdGradientMicroDeformation[ i ][ j ]
                        = DSolutionVectorDDeformation[ 45 * ( i + 9 ) + j + 18 ];

                    //The stress measures
                    dPK2StressdGradientMicroDeformation[ i ][ j ]            = DStressDDeformation[ 45 * i + j + 18 ];
                    dReferenceMicroStressdGradientMicroDeformation[ i ][ j ] = DStressDDeformation[ 45 * ( i + 9 ) + j + 18 ];
                }
            }

            for ( unsigned int i = 0; i < 27; i++ ){
                for ( unsigned int j = 0; j < 9; j++ ){

                    //The plastic deformation measures
                    dPlasticGradientMicroDeformationdDeformationGradient[ i ][ j ]
                        = DSolutionVectorDDeformation[ 45 * ( i + 18 ) + j ];
                    dPlasticGradientMicroDeformationdMicroDeformation[ i ][ j ]
                        = DSolutionVectorDDeformation[ 45 * ( i + 18 ) + j + 9 ];

                    //The stress measures
                    dReferenceHigherOrderStressdDeformationGradient[ i ][ j ]
                        = DStressDDeformation[ 45 * ( i + 18 ) + j ];
                    dReferenceHigherOrderStressdMicroDeformation[ i ][ j ]
                        = DStressDDeformation[ 45 * ( i + 18 ) + j + 9 ];
                }

                for ( unsigned int j = 0; j < 27; j++ ){

                    //The plastic deformation measures
                    dPlasticGradientMicroDeformationdGradientMicroDeformation[ i ][ j ]
                        = DSolutionVectorDDeformation[ 45 * ( i + 18 ) + j + 18 ];

                    //The stress measures
                    dReferenceHigherOrderStressdGradientMicroDeformation[ i ][ j ]
                        = DStressDDeformation[ 45 * ( i + 18 ) + j + 18 ];
                }
            }
        }
        else{

            dPK2StressdDeformationGradient = vectorTools::dot( dPK2StressdElasticDeformationGradient,
                                                               dElasticDeformationGradientdDeformationGradient );
            dPK2StressdMicroDeformation = vectorTools::dot( dPK2StressdElasticMicroDeformation,
                                                            dElasticMicroDeformationdMicroDeformation )
                                        + vectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                                            dElasticGradientMicroDeformationdMicroDeformation );
            dPK2StressdGradientMicroDeformation = vectorTools::dot( dPK2StressdElasticGradientMicroDeformation,
                                                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            dReferenceMicroStressdDeformationGradient = vectorTools::dot( dReferenceMicroStressdElasticDeformationGradient,
                                                                          dElasticDeformationGradientdDeformationGradient );
            dReferenceMicroStressdMicroDeformation = vectorTools::dot( dReferenceMicroStressdElasticMicroDeformation,
                                                                       dElasticMicroDeformationdMicroDeformation )
                                                   + vectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                                                       dElasticGradientMicroDeformationdMicroDeformation );
            dReferenceMicroStressdGradientMicroDeformation
                = vectorTools::dot( dReferenceMicroStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdGradientMicroDeformation );

            dReferenceHigherOrderStressdDeformationGradient = vectorTools::dot( dReferenceHigherOrderStressdElasticDeformationGradient,
                                                                                dElasticDeformationGradientdDeformationGradient );
            dReferenceHigherOrderStressdMicroDeformation = vectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                                                             dElasticGradientMicroDeformationdMicroDeformation );
            dReferenceHigherOrderStressdGradientMicroDeformation
                = vectorTools::dot( dReferenceHigherOrderStressdElasticGradientMicroDeformation,
                                    dElasticGradientMicroDeformationdGradientMicroDeformation );


        }

#ifdef DEBUG_MODE
        //Save the total derivatives of the plastic deformation measures w.r.t. the fundamental deformation measures
        tempDEBUG.emplace( "totaldPlasticDeformationGradientdDeformationGradient",
                            vectorTools::appendVectors( dPlasticDeformationGradientdDeformationGradient ) );
        tempDEBUG.emplace( "totaldPlasticDeformationGradientdMicroDeformation",
                            vectorTools::appendVectors( dPlasticDeformationGradientdMicroDeformation ) );
        tempDEBUG.emplace( "totaldPlasticDeformationGradientdGradientMicroDeformation",
                            vectorTools::appendVectors( dPlasticDeformationGradientdGradientMicroDeformation ) );

        tempDEBUG.emplace( "totaldPlasticMicroDeformationdDeformationGradient",
                            vectorTools::appendVectors( dPlasticMicroDeformationdDeformationGradient ) );
        tempDEBUG.emplace( "totaldPlasticMicroDeformationdMicroDeformation",
                            vectorTools::appendVectors( dPlasticMicroDeformationdMicroDeformation ) );
        tempDEBUG.emplace( "totaldPlasticMicroDeformationdGradientMicroDeformation",
                            vectorTools::appendVectors( dPlasticMicroDeformationdGradientMicroDeformation ) );

        tempDEBUG.emplace( "totaldPlasticGradientMicroDeformationdDeformationGradient",
                            vectorTools::appendVectors( dPlasticGradientMicroDeformationdDeformationGradient ) );
        tempDEBUG.emplace( "totaldPlasticGradientMicroDeformationdMicroDeformation",
                            vectorTools::appendVectors( dPlasticGradientMicroDeformationdMicroDeformation ) );
        tempDEBUG.emplace( "totaldPlasticGradientMicroDeformationdGradientMicroDeformation",
                            vectorTools::appendVectors( dPlasticGradientMicroDeformationdGradientMicroDeformation ) );

        //Save the total derivatives of the intermediate stresses w.r.t. the fundamental deformation measures
        tempDEBUG.emplace( "totaldPK2StressdDeformationGradient",
                            vectorTools::appendVectors( dPK2StressdDeformationGradient ) );
        tempDEBUG.emplace( "totaldPK2StressdMicroDeformation",
                            vectorTools::appendVectors( dPK2StressdMicroDeformation ) );
        tempDEBUG.emplace( "totaldPK2StressdGradientMicroDeformation",
                            vectorTools::appendVectors( dPK2StressdGradientMicroDeformation ) );

        tempDEBUG.emplace( "totaldReferenceMicroStressdDeformationGradient",
                            vectorTools::appendVectors( dReferenceMicroStressdDeformationGradient ) );
        tempDEBUG.emplace( "totaldReferenceMicroStressdMicroDeformation",
                            vectorTools::appendVectors( dReferenceMicroStressdMicroDeformation ) );
        tempDEBUG.emplace( "totaldReferenceMicroStressdGradientMicroDeformation",
                            vectorTools::appendVectors( dReferenceMicroStressdGradientMicroDeformation ) );

        tempDEBUG.emplace( "totaldReferenceHigherOrderStressdDeformationGradient",
                           vectorTools::appendVectors( dReferenceHigherOrderStressdDeformationGradient ) );
        tempDEBUG.emplace( "totaldReferenceHigherOrderStressdMicroDeformation",
                           vectorTools::appendVectors( dReferenceHigherOrderStressdMicroDeformation ) );
        tempDEBUG.emplace( "totaldReferenceHigherOrderStressdGradientMicroDeformation",
                           vectorTools::appendVectors( dReferenceHigherOrderStressdGradientMicroDeformation ) );

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
        currentStrainISVS = vectorTools::appendVectors( { currentStrainISVS, currentMicroGradientGamma } );

        std::vector< double > currentGammas = { currentMacroGamma, currentMicroGamma };
        currentGammas = vectorTools::appendVectors( { currentGammas, currentMicroGradientGamma } );

        SDVS = vectorTools::appendVectors( { currentStrainISVS, currentGammas, currentPlasticDeformationGradient - eye,
                                             currentPlasticMicroDeformation - eye, currentPlasticGradientMicroDeformation } );

        //Output the stresses mapped to the true reference configuration
        variableMatrix dPK2dIntermediatePK2, dPK2dPlasticDeformationGradient;
        error = micromorphicTools::pullBackCauchyStress( currentPK2Stress, currentPlasticDeformationGradient, current_PK2,
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
            = vectorTools::dot( dPK2dIntermediatePK2, dPK2StressdDeformationGradient )
            + vectorTools::dot( dPK2dPlasticDeformationGradient, dPlasticDeformationGradientdDeformationGradient );

        variableMatrix dReferencePK2StressdMicroDeformation
            = vectorTools::dot( dPK2dIntermediatePK2, dPK2StressdMicroDeformation )
            + vectorTools::dot( dPK2dPlasticDeformationGradient, dPlasticDeformationGradientdMicroDeformation );

        variableMatrix dReferencePK2StressdGradientMicroDeformation
            = vectorTools::dot( dPK2dIntermediatePK2, dPK2StressdGradientMicroDeformation )
            + vectorTools::dot( dPK2dPlasticDeformationGradient, dPlasticDeformationGradientdGradientMicroDeformation );

        variableMatrix dSIGMAdIntermediateSIGMA, dSIGMAdPlasticDeformationGradient;
        error = micromorphicTools::pullBackMicroStress( currentReferenceMicroStress, currentPlasticDeformationGradient, current_SIGMA,
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
            = vectorTools::dot( dSIGMAdIntermediateSIGMA, dReferenceMicroStressdDeformationGradient )
            + vectorTools::dot( dSIGMAdPlasticDeformationGradient, dPlasticDeformationGradientdDeformationGradient );

        variableMatrix dSIGMAStressdMicroDeformation
            = vectorTools::dot( dSIGMAdIntermediateSIGMA, dReferenceMicroStressdMicroDeformation )
            + vectorTools::dot( dSIGMAdPlasticDeformationGradient, dPlasticDeformationGradientdMicroDeformation );

        variableMatrix dSIGMAStressdGradientMicroDeformation
            = vectorTools::dot( dSIGMAdIntermediateSIGMA, dReferenceMicroStressdGradientMicroDeformation )
            + vectorTools::dot( dSIGMAdPlasticDeformationGradient, dPlasticDeformationGradientdGradientMicroDeformation );

        variableMatrix dMdIntermediateM, dMdPlasticDeformationGradient, dMdPlasticMicroDeformation;
        error = micromorphicTools::pullBackHigherOrderStress( currentReferenceHigherOrderStress, 
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
            = vectorTools::dot( dMdIntermediateM, dReferenceHigherOrderStressdDeformationGradient )
            + vectorTools::dot( dMdPlasticDeformationGradient, dPlasticDeformationGradientdDeformationGradient )
            + vectorTools::dot( dMdPlasticMicroDeformation, dPlasticMicroDeformationdDeformationGradient );

        variableMatrix dMStressdMicroDeformation
            = vectorTools::dot( dMdIntermediateM, dReferenceHigherOrderStressdMicroDeformation )
            + vectorTools::dot( dMdPlasticDeformationGradient, dPlasticDeformationGradientdMicroDeformation )
            + vectorTools::dot( dMdPlasticMicroDeformation, dPlasticMicroDeformationdMicroDeformation );

        variableMatrix dMStressdGradientMicroDeformation
            = vectorTools::dot( dMdIntermediateM, dReferenceHigherOrderStressdGradientMicroDeformation )
            + vectorTools::dot( dMdPlasticDeformationGradient, dPlasticDeformationGradientdGradientMicroDeformation )
            + vectorTools::dot( dMdPlasticMicroDeformation, dPlasticMicroDeformationdGradientMicroDeformation );

        //Assemble the Jacobians w.r.t. the degrees of freedom
        DPK2Dgrad_u   = vectorTools::dot( dReferencePK2StressdDeformationGradient,
                                          dDeformationGradientdGradientMacroDisplacement );
        DPK2Dphi      = vectorTools::dot( dReferencePK2StressdMicroDeformation,
                                          dMicroDeformationdMicroDisplacement );
        DPK2Dgrad_phi = vectorTools::dot( dReferencePK2StressdGradientMicroDeformation,
                                          dGradientMicroDeformationdGradientMicroDisplacement );

        DSIGMADgrad_u   = vectorTools::dot( dSIGMAStressdDeformationGradient, dDeformationGradientdGradientMacroDisplacement );
        DSIGMADphi      = vectorTools::dot( dSIGMAStressdMicroDeformation, dMicroDeformationdMicroDisplacement );
        DSIGMADgrad_phi = vectorTools::dot( dSIGMAStressdGradientMicroDeformation, dGradientMicroDeformationdGradientMicroDisplacement );

        DMDgrad_u   = vectorTools::dot( dMStressdDeformationGradient, dDeformationGradientdGradientMacroDisplacement );
        DMDphi      = vectorTools::dot( dMStressdMicroDeformation, dMicroDeformationdMicroDisplacement );
        DMDgrad_phi = vectorTools::dot( dMStressdGradientMicroDeformation, dGradientMicroDeformationdGradientMicroDisplacement );

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
            error = micromorphicLinearElasticity::formIsotropicA( outputs[ 9 ][ 0 ], outputs[ 9 ][ 1 ], Amatrix );
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
            error = micromorphicLinearElasticity::formIsotropicB( outputs[ 10 ][ 0 ], outputs[ 10 ][ 1 ], outputs[ 10 ][ 2 ],
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
            error = micromorphicLinearElasticity::formIsotropicC( outputs[ 11 ], Cmatrix );
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
            error = micromorphicLinearElasticity::formIsotropicD( outputs[ 12 ][ 0 ], outputs[ 12 ][ 1 ], Dmatrix );
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
        vectorTools::eye( eye );

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

        errorOut error = micromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures",
                                             "Error in assembly of the deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures",
                                             "Error in assembly of the micro deformation" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation );

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

        errorOut error = micromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient, dFdGradU );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures (jacobian)",
                                             "Error in assembly of the deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation, dChidPhi );

        if ( error ){
            errorOut result = new errorNode( "assembleFundamentalDeformationMeasures (jacobian)",
                                             "Error in assembly of the micro deformation" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation,
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
                                     , solverTools::debugMap &DEBUG
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
         * :param std::map< std::string, solverTools::floatVector > &DEBUG: The debug map. Only output when in debug mode.
         */

        if ( macroYieldParameters.size() != 2 ){
            std::cout << "macroYieldParameters: "; vectorTools::print( macroYieldParameters );
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
        solverTools::floatVector tmp = { yieldFunctionValues[ 0 ] };
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
                                     , solverTools::debugMap &DEBUG
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
         * :param std::map< std::string, solverTools::floatVector > &DEBUG: The debug map. Only output when in debug mode.
         */

        if ( macroYieldParameters.size() != 2 ){
            std::cout << "macroYieldParameters: "; vectorTools::print( macroYieldParameters );
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
        solverTools::floatVector tmp = { yieldFunctionValues[ 0 ] };
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
        if ( vectorTools::dot( PK2Stress, PK2Stress ) < 1e-9 ){
            dMacroFdPK2        = variableVector( PK2Stress.size(), 0 );
            dMacroFdElasticRCG = variableVector( elasticRightCauchyGreen.size(), 0 );
        }

        if ( vectorTools::dot( referenceMicroStress, referenceMicroStress ) < 1e-9 ){
            dMicroFdSigma      = variableVector( referenceMicroStress.size(), 0 );
            dMicroFdElasticRCG = variableVector( elasticRightCauchyGreen.size(), 0 );
        }

        if ( vectorTools::dot( referenceHigherOrderStress, referenceHigherOrderStress ) < 1e-9 ){
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
        constantMatrix eye = vectorTools::eye< constantType >( microGradientStrainISV.size() );

        dMacroCdMacroStrainISV = macroHardeningParameters[ 1 ];
        dMicroCdMicroStrainISV = microHardeningParameters[ 1 ];
        dMicroGradientCdMicroGradientStrainISV = microGradientHardeningParameters[ 1 ] * eye;

        return NULL;
    }
}
