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
                                                           variableVector &dFdElasticRCG ){
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
         */

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
        variableVector devStressDirection = deviatoricReferenceStress / normDevStress;

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
                                                           variableMatrix &d2FdStressdElasticRCG ){
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
         */
        //Assume 3D
        unsigned int dim = 3;

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
        variableVector devStressDirection = deviatoricReferenceStress / normDevStress;

        dFdStress = vectorTools::Tdot( dDevStressdStress, devStressDirection )
                  + BAngle * dPressuredStress;

        dFdc = - AAngle;

        dFdElasticRCG = vectorTools::Tdot( dDevStressdRCG, devStressDirection )
                      + BAngle * dPressuredRCG;

        //Evaluate the second-order jacobians
        constantMatrix EYE = vectorTools::eye< constantType >( dim * dim );
        variableMatrix dDevStressDirectiondDevStress = ( EYE - vectorTools::dyadic( devStressDirection, devStressDirection ) ) / normDevStress;

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


    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma,
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableMatrix &microGradientFlowDirection,
                                              variableVector &plasticMacroVelocityGradient, variableVector &plasticMicroVelocityGradient,
                                              variableVector &plasticMicroGradientVelocityGradient ){
        /*!
         * Compute the plastic velocity gradients in the intermediate configuration.
         *
         * \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         * \bar{ L }_{ \bar{B} \bar{K} }^{ \chi, p } = \dot{ \bar{ \gamma } }^{\chi} \Psi_{ \bar{E} \bar{L} }^{e, -1} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } } \bar{ \Psi }_{ \bar{N} \bar{K} }^{e} \bar{C}_{ \bar{N} \bar{F} }^{\chi, p, -1 }
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
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
         * :param const variableMatrix &microGradientFlowDirection: The flow direction of the micro gradient plasticity.
         *     Note: This is a matrix because it is computed as the gradient of the flow potential which is a vector.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &plasticMicroVelocityGradient: The plastic velocity gradient for the micro plastic deformation.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic velocity gradient for the micro gradient 
         *     plastic deformation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( microGradientGamma.size() != dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The micro gradient gamma must be of dimension 3" );
        }

        if ( elasticRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The elastic right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( elasticMicroRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The elastic micro right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The elastic micro deformation metric Psi must be 3D" );
        }

        if ( elasticGamma.size() != dim * dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The elastic higher order deformation metric Gamma must be 3D" );
        }

        if ( macroFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The macro flow direction must be 3D" );
        }

        if ( microFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The micro flow direction must be 3D" );
        }

        if ( microGradientFlowDirection.size() != dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The micro gradient flow direction must be 3D" );
        }

        for ( unsigned int i = 0; i < dim; i++ ){
            if ( microGradientFlowDirection[ i ].size() != dim * dim * dim ){
                return new errorNode( "computePlasticVelocityGradients",
                                      "The rows of the micro gradient flow direction must be of length 27" );
            }
        }

        //Compute the required inverses of the deformation metrics
        variableVector inverseElasticRightCauchyGreen = vectorTools::inverse( elasticRightCauchyGreen, dim, dim );
//        variableVector inverseElasticMicroRightCauchyGreen = vectorTools::inverse( elasticMicroRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = vectorTools::inverse( elasticPsi, dim, dim );

        //Compute the macro-scale velocity gradient

        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients",
                                             "Error in computation of plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        error = computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen, elasticPsi, inverseElasticPsi,
                                                     microFlowDirection, plasticMicroVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients",
                                             "Error in computation of plastic micro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Assemble the 'skew' term
        variableVector skewTerm = variableVector( dim * dim * dim, 0 );

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
                                 * ( microGradientGamma[ Ib ] * microGradientFlowDirection[ Ib ][ dim * dim * Kb + dim * Lb + Mb ]
                                 +   elasticPsi[ dim * Lb + Ib ] * skewTerm[ dim * dim * Ib + dim * Mb + Kb ] );
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
                                              const variableVector &microFlowDirection, const variableMatrix &microGradientFlowDirection,
                                              variableVector &macroPlasticVelocityGradient, variableVector &microPlasticVelocityGradient,
                                              variableVector &microGradientPlasticVelocityGradient, variableVector &dMacroLpdMacroGamma,                                               variableVector &dMacroLpdMicroGamma, variableVector &dMicroLpdMicroGamma,
                                              variableVector &dMicroGradientLpdMacroGamma, variableVector &dMicroGradientLpdMicroGamma,
                                              variableMatrix &dMicroGradientLpdMicroGradientGamma ){

    }
}
