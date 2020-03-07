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

        variableVector elasticRightCauchyGreenInverse;
        return computeSecondOrderDruckerPragerYieldEquation( referenceStressMeasure, cohesion, elasticRightCauchyGreen, frictionAngle,
                                                             beta, yieldValue, elasticRightCauchyGreenInverse );
    }

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &elasticRightCauchyGreenInverse ){
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
         * :param variableVector &elasticRightCauchyGreenInverse: The inverse of the Right Cauchy-Green
         *     elastic deformation tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Compute the parameters
        parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the pressure
        variableType pressure;
        errorOut error = micromorphicTools::computeReferenceSecondOrderStressPressure( referenceStressMeasure, 
                             elasticRightCauchyGreen, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation",
                                             "Error in computation of second-order stress pressure" );
            result->addNext( error );
            return result;
        }

        //Invert the elastic Right Cauchy-Green deformation tensor
        elasticRightCauchyGreenInverse = vectorTools::inverse( elasticRightCauchyGreen, dim, dim );

        //Compute the deviatoric stress
        variableVector deviatoricReferenceStress = 
            referenceStressMeasure - pressure * elasticRightCauchyGreenInverse;

        //Compute the l2norm of the deviatoric stress
        variableType normDev = vectorTools::l2norm( deviatoricReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDev - ( AAngle * cohesion - BAngle * pressure );

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
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         *  Also compute the Jacobians
         *  \frac{ \partial F }{ \partial stressMeasure_{IJ} } = \frac{ dev ( stressMeasure_{IJ} ) }{ || dev ( stressMeasure ) || + B^{\phi} \frac{ \partial \bar{p} }{ \partial stressMeasure_{IJ} }
         *  \frac{ \partial F }{ \partial RCG_{IJ} } = B^{\phi} \frac{ \partial \bar{p} }{ \partial RCG_{IJ} }
         *  \frac{ \partial F }{ cohesion } = -A^{\phi}
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
         */

        variableVector dNormDevdDevReferenceStress;
        variableMatrix dDevStressdStress, dDevStressdElasticRCG;
        variableType normDevStress, Bphi;

        return computeSecondOrderDruckerPragerYieldEquation( referenceStressMeasure, cohesion, elasticRightCauchyGreen,
                                                             frictionAngle, beta, yieldValue, dFdStress, dFdc,
                                                             dFdElasticRCG, dNormDevdDevReferenceStress, 
                                                             dDevStressdStress, dDevStressdElasticRCG, normDevStress, Bphi );

    }

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdElasticRCG, variableVector &dNormDevdDevReferenceStress, 
                                                           variableMatrix &dDevStressdReferenceStress,
                                                           variableMatrix &dDevStressdElasticRCG,
                                                           variableType &normDevStress, parameterType &BAngle ){
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
         *  \frac{ \partial F }{ \partial stressMeasure_{IJ} } = \frac{ dev ( stressMeasure_{IJ} ) }{ || dev ( stressMeasure ) || + B^{\phi} \frac{ \partial \bar{p} }{ \partial stressMeasure_{IJ} }
         *  \frac{ \partial F }{ \partial RCG_{IJ} } = B^{\phi} \frac{ \partial \bar{p} }{ \partial RCG_{IJ} }
         *  \frac{ \partial F }{ cohesion } = -A^{\phi}
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
         * :param variableVector &dNormDevdDevReferenceStress: The Jacobian of the norm of the deviatoric stress w.r.t.
         *     the reference stress.
         * :param variableMatrix &dDevStressdReferenceStress: The Jacobian of the deviatoric stress w.r.t. the reference stress.
         * :param variableMatrix &dDevStressdElasticRCG: The Jacobian of the deviatoric stress w.r.t. the elastic part of the right 
         *     cauchy green deformation tensor.
         * :param variableType &normDevStress: The norm of the deviatoric part of the stress.
         * :param parameterType &BAngle: The BAngle parameter. 
         */

        //Assume 3D
        unsigned int dim = 3;

        //Compute the parameters
        parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the pressure
        variableType pressure;
        variableVector dpdStress, dpdElasticRCG;
        errorOut error = micromorphicTools::computeReferenceSecondOrderStressPressure( referenceStressMeasure,
                             elasticRightCauchyGreen, pressure, dpdStress, dpdElasticRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation",
                                             "Error in computation of second-order stress pressure" );
            result->addNext( error );
            return result;
        }

        //Invert the elastic Right Cauchy-Green deformation tensor
        variableVector elasticRightCauchyGreenInverse = vectorTools::inverse( elasticRightCauchyGreen, dim, dim );

        variableMatrix dElasticRCGInversedElasticRCG( elasticRightCauchyGreenInverse.size(),
                           variableVector( elasticRightCauchyGreen.size(), 0 ) );

        //Compute the Jacobian of the inverse of the Right Cauchy-Green deformation tensor with 
        //respect to itself
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        dElasticRCGInversedElasticRCG[ dim * I + J ][ dim * K + L ] -= elasticRightCauchyGreenInverse[ dim * I + K ]
                                                                                     * elasticRightCauchyGreenInverse[ dim * L + J ];
                    }
                }
            }
        }

        //Compute the deviatoric stress
        variableVector deviatoricReferenceStress = 
            referenceStressMeasure - pressure * elasticRightCauchyGreenInverse;

        //Compute the Jacobian of the deviatoric stress
        constantMatrix EYE = vectorTools::eye< variableType >( dim * dim );
        dDevStressdReferenceStress = EYE - vectorTools::dyadic( elasticRightCauchyGreenInverse, dpdStress );
        dDevStressdElasticRCG = - ( vectorTools::dyadic( elasticRightCauchyGreenInverse, dpdElasticRCG )
                                  + pressure * dElasticRCGInversedElasticRCG );

        //Compute the l2norm of the deviatoric stress
        normDevStress = vectorTools::l2norm( deviatoricReferenceStress );

        //Compute the Jacobian of the l2norm of the deviatoric stress
        dNormDevdDevReferenceStress = deviatoricReferenceStress / normDevStress;
        variableVector dNormDevdReferenceStress = vectorTools::Tdot( dDevStressdReferenceStress, dNormDevdDevReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Compute the jacobian
        dFdStress = dNormDevdReferenceStress + BAngle * dpdStress;
        dFdc = -AAngle;
        dFdElasticRCG = vectorTools::Tdot( dDevStressdElasticRCG, dNormDevdDevReferenceStress ) + BAngle * dpdElasticRCG;

        return NULL;
    }

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdElasticRCG, variableMatrix &d2FdStress2,
                                                           variableMatrix d2FdStressdElasticRCG ){
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
         *  \frac{ \partial F }{ \partial stressMeasure_{IJ} } = \frac{ dev ( stressMeasure_{IJ} ) }{ || dev ( stressMeasure ) || + B^{\phi} \frac{ \partial \bar{p} }{ \partial stressMeasure_{IJ} }
         *  \frac{ \partial F }{ \partial RCG_{IJ} } = B^{\phi} \frac{ \partial \bar{p} }{ \partial RCG_{IJ} }
         *  \frac{ \partial F }{ cohesion } = -A^{\phi}
         *
         *  The second deriatives of \frac{ \partial F }{ \partial \stressMeasure_{IJ} } are
         *  \frac{ \partial^2 F }{ \partial stressMeasure_{IJ} \partial stressMeasure_{KL} } = \frac{1}{ || dev ( stressMeasure ) || } \left ( \delta_{IM} \delta{JN} - \frac{ dev( stressMeasure_{IJ} ) }{ || dev ( stressMeasure ) || } \frac{ stressMeasure_{MN} }{ || dev ( stressMeasure ) || } \right ) \frac{ \partial dev ( stressMeasure )_{MN} }{ \partial stressMeasure_{KL} }
         *  \frac{ \partial^2 F }{ \partial stressMeasure_{IJ} \partial RCG_{KL} } = \frac{1}{ || dev ( stressMeasure ) || } \left ( \delta_{IM} \delta{JN} - \frac{ dev( stressMeasure_{IJ} ) }{ || dev ( stressMeasure ) || } \frac{ stressMeasure_{MN} }{ || dev ( stressMeasure ) || } \right ) \frac{ \partial dev ( stressMeasure )_{MN} }{ \partial RCG_{KL} } + B^{\phi} \delta_{IK} \delta_{JL}

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

        variableVector dNormDevdDevReferenceStress;
        variableMatrix dDevStressdStress, dDevStressdElasticRCG;
        variableType normDevStress;
        parameterType BAngle;
        errorOut error = computeSecondOrderDruckerPragerYieldEquation( referenceStressMeasure, cohesion,
                                                                       elasticRightCauchyGreen, frictionAngle, beta,
                                                                       yieldValue, dFdStress, dFdc, dFdElasticRCG,
                                                                       dNormDevdDevReferenceStress, dDevStressdStress,
                                                                       dDevStressdElasticRCG, normDevStress, BAngle );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation (flow direction jacobian)",
                                             "Error in computation of the jacobians of the yield surface" );
            result->addNext( error );
            return result;
        }

        constantMatrix EYE = vectorTools::eye< constantType >( dim * dim );
        variableMatrix d2FdDeviatoricStress2 = ( EYE - vectorTools::dyadic( dNormDevdDevReferenceStress, dNormDevdDevReferenceStress ) ) / normDevStress;

        d2FdStress2 = vectorTools::dot( d2FdDeviatoricStress2, dDevStressdStress );
        d2FdStressdElasticRCG = vectorTools::dot( d2FdDeviatoricStress2, dDevStressdElasticRCG ) + BAngle * EYE / 3;

        return NULL;
    }
}
