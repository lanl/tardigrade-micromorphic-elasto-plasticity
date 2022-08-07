//Tests for constitutive_tools

#include<micromorphic_elasto_plasticity.h>
#include<sstream>
#include<fstream>
#include<iostream>
#include<iomanip>

#define BOOST_TEST_MODULE test_micromorphic_linear_elasticity
#include <boost/test/included/unit_test.hpp>

typedef micromorphicTools::constantType constantType;
typedef micromorphicTools::constantVector constantVector;
typedef micromorphicTools::constantMatrix constantMatrix;

typedef micromorphicTools::parameterType parameterType;
typedef micromorphicTools::parameterVector parameterVector;
typedef micromorphicTools::parameterMatrix parameterMatrix;

typedef micromorphicTools::variableType variableType;
typedef micromorphicTools::variableVector variableVector;
typedef micromorphicTools::variableMatrix variableMatrix;

typedef micromorphicTools::errorNode errorNode;
typedef micromorphicTools::errorOut errorOut;

BOOST_AUTO_TEST_CASE( testComputeSecondOrderDruckerPragerYieldEquation ){
    /*!
     * Test the computation of the second order stress Drucker-Prager yield equation.
     *
     */

    variableVector S = { 0.77189588, -0.84417528,  0.95929231,
                        -0.50465708,  0.50576944,  0.05335127,
                         0.81510751,  0.76814059, -0.82146208 };

    variableVector C = { 0.03468919, -0.31275742, -0.57541261,
                        -0.27865312, -0.45844965,  0.52325004,
                        -0.0439162 , -0.80201065, -0.44921044 };

    variableType cohesion = 0.51;
    parameterType frictionAngle = 0.46;
    parameterType beta = 0.34;

    constantType answer = 3.236807543277929;
    variableType result;

    errorOut error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, C,
                                                                                                 frictionAngle, beta, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    //Test the Jacobian
    variableType resultJ, dFdCohesion;
    variableVector dFdS, dFdC;

    error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, C,
                                                                                        frictionAngle, beta, resultJ,
                                                                                        dFdS, dFdCohesion, dFdC );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultJ, answer ) );

    //Test the Jacobian
    variableType resultJ2, dFdcohesionJ2;
    variableVector dFdSJ2, dFdCJ2;
    variableMatrix d2FdStress2J2;
    variableMatrix d2FdStressdElasticRCGJ2;
    variableMatrix d2FdS2J2, d2FdSdCJ2;

    error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, C, 
                                                                                        frictionAngle, beta, resultJ2,
                                                                                        dFdSJ2, dFdcohesionJ2, dFdCJ2, d2FdS2J2, 
                                                                                        d2FdSdCJ2);

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultJ2, answer ) );

    //Test dFdStress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S + delta, cohesion, C,
                                                                                            frictionAngle, beta, resultJ );

        BOOST_CHECK( !error );

        constantType gradCol = ( resultJ - result ) / delta[i];

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dFdS[i] ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dFdSJ2[i] ) );
    }

    //Test dFdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, C + delta,
                                                                                            frictionAngle, beta, resultJ );

        BOOST_CHECK( !error );

        constantType gradCol = ( resultJ - result ) / delta[i];

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dFdC[i], 1e-4 ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dFdCJ2[i], 1e-4 ) );
    }

    //Test dFdcohesion
    constantType deltas = eps * fabs( cohesion ) + eps;

    error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion + deltas, C, 
                                                                                        frictionAngle, beta, resultJ );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( ( resultJ - result ) / deltas, dFdcohesionJ2 ) );

    //Test d2FdStress2
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        variableVector dFdSp, dFdSm, dFdCp, dFdCm;
        variableType dFdcohesionp, dFdcohesionm;

        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S + delta, cohesion, C,
                                                                                            frictionAngle, beta, resultJ,
                                                                                            dFdSp, dFdcohesionp, dFdCp );

        BOOST_CHECK( !error );
                
        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S - delta, cohesion, C,
                                                                                            frictionAngle, beta, resultJ,
                                                                                            dFdSm, dFdcohesionm, dFdCm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( dFdSp - dFdSm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], d2FdS2J2[j][i] ) );
        }
    }

    //Test d2FdStressdElasticRCG
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableVector dFdSp, dFdSm, dFdCp, dFdCm;
        variableType dFdcohesionp, dFdcohesionm;

        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, C + delta,
                                                                                            frictionAngle, beta, resultJ,
                                                                                            dFdSp, dFdcohesionp, dFdCp );

        BOOST_CHECK( !error );
                
        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, C - delta,
                                                                                            frictionAngle, beta, resultJ,
                                                                                            dFdSm, dFdcohesionm, dFdCm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( dFdSp - dFdSm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], d2FdSdCJ2[j][i] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputeHigherOrderDruckerPragerYieldEquation ){
    /*!
     * Test the computation of the higher order stress Drucker-Prager yield equation.
     *
     */

    variableVector M = { 0.80732114,  0.79202055,  0.17990022,  0.97454675,  0.703207  ,
                        -0.58236697,  0.53324571, -0.93438873, -0.40650796,  0.14071918,
                         0.66933708, -0.67854069, -0.30317772, -0.93821882,  0.97270622,
                         0.00295302, -0.12441126,  0.30539971, -0.0580227 ,  0.89696105,
                         0.17567709, -0.9592962 ,  0.63535407,  0.95437804, -0.64531877,
                         0.69978907,  0.81327586 };

    variableVector C = { 0.03468919, -0.31275742, -0.57541261,
                        -0.27865312, -0.45844965,  0.52325004,
                        -0.0439162 , -0.80201065, -0.44921044 };

    variableVector cohesion = { 0.51, 0.71, .14 };
    parameterType frictionAngle = 0.46;
    parameterType beta = 0.34;

    constantVector answer = { 3.90579607, 0.86300970, 4.63923709 };
    variableVector result;

    errorOut error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C,
                                                                                                 frictionAngle, beta, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    //Test the Jacobians

    variableVector resultJ;
    variableMatrix dFdStress, dFdc, dFdRCG;

    error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C,
                                                                                        frictionAngle, beta, resultJ,
                                                                                        dFdStress, dFdc, dFdRCG );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultJ, answer ) );

    variableVector resultJ2;
    variableMatrix dFdStressJ2, dFdcJ2, dFdRCGJ2;
    variableMatrix d2FdStress2J2, d2FdStressdRCGJ2;

    error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C,
                                                                                        frictionAngle, beta, resultJ2,
                                                                                        dFdStressJ2, dFdcJ2, dFdRCGJ2,
                                                                                        d2FdStress2J2, d2FdStressdRCGJ2 );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultJ2, answer ) );
    
    //Test derivatives w.r.t stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M + delta, cohesion, C,
                                                                                            frictionAngle, beta, resultP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M - delta, cohesion, C,
                                                                                            frictionAngle, beta, resultM );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dFdStress[j][i] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dFdStressJ2[j][i] ) );
        }

        variableMatrix dFdStressP, dFdcP, dFdRCGP;
        variableMatrix dFdStressM, dFdcM, dFdRCGM;

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M + delta, cohesion, C,
                                                                                            frictionAngle, beta, resultP,
                                                                                            dFdStressP, dFdcP, dFdRCGP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M - delta, cohesion, C,
                                                                                            frictionAngle, beta, resultM,
                                                                                            dFdStressM, dFdcM, dFdRCGM );

        BOOST_CHECK( !error );


        constantMatrix gradMat = ( dFdStressP - dFdStressM ) / ( 2 * delta[i] );

        unsigned int n, o, p;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        n = ( int )( i / 9 );
                        o = ( int )( (i - 9 * n ) / 3 );
                        p = ( i - 9 * n - 3 * o ) % 3;
                        BOOST_CHECK( vectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                        d2FdStress2J2[ j ][ 243 * k + 81 * l + 27 * m + 9 * n + 3 * o + p ] ) );
                    }
                }
            }
        }
    }

    //Test derivatives w.r.t. the cohesion
    for ( unsigned int i = 0; i < cohesion.size(); i++ ){
        constantVector delta( cohesion.size(), 0 );
        delta[i] = eps * fabs( cohesion[i] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion + delta, C,
                                                                                            frictionAngle, beta, resultP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion - delta, C,
                                                                                            frictionAngle, beta, resultM );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dFdc[j][i] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dFdcJ2[j][i] ) );
        }
    }

    //Test derivatives w.r.t. the right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C + delta,
                                                                                            frictionAngle, beta, resultP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C - delta,
                                                                                            frictionAngle, beta, resultM );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dFdRCG[j][i] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dFdRCG[j][i] ) );
        }

        variableMatrix dFdStressP, dFdcP, dFdRCGP;
        variableMatrix dFdStressM, dFdcM, dFdRCGM;

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C + delta,
                                                                                            frictionAngle, beta, resultP,
                                                                                            dFdStressP, dFdcP, dFdRCGP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C - delta,
                                                                                            frictionAngle, beta, resultM,
                                                                                            dFdStressM, dFdcM, dFdRCGM );

        BOOST_CHECK( !error );

        constantMatrix gradMat = ( dFdStressP - dFdStressM ) / ( 2 * delta[i] );

        unsigned int n, o;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        n = ( int )( i / 3 );
                        o = ( i % 3 );
                        BOOST_CHECK( vectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                        d2FdStressdRCGJ2[ j ][ 81 * k + 27 * l + 9 * m + 3 * n + o ] ) );
                    }
                }
            }
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputeElasticPartOfDeformation ){
    /*!
     * Test of the computation of the elastic part of the various deformation measures.
     *
     */

    //Define the full deformation

    variableVector F = { -1.08831037, -0.66333427, -0.48239487,
                         -0.904554  ,  1.28942848, -0.02156112,
                         -0.08464824, -0.07730218,  0.86415668 };

    variableVector chi = { 0.27548762,  0.71704571, -1.33727026,
                          -0.42283117, -1.65591182, -0.22625105,
                           1.13780238, -0.51375234, -0.62058393 };

    variableVector gradChi = { 0.99727481, 0.19819654, 0.20010236, 0.17532525, 0.06633508,
                               0.72458529, 0.75509811, 0.25826068, 0.49724274, 0.92781534,
                               0.11869579, 0.8604985 , 0.17654448, 0.33682062, 0.17328738,
                               0.89505449, 0.93584041, 0.07746146, 0.96392792, 0.41350595,
                               0.50880671, 0.3957155 , 0.29743433, 0.65458412, 0.85345184,
                               0.09451065, 0.60621985 };

    variableVector Fp = { 0.28777732,  0.09144675,  0.22352836,
                         -0.27852895,  0.88746541, -0.77925405,
                          0.08833233,  0.02544955, -0.29848719 };

    variableVector chip = { -0.49993961, -0.36238477, -0.24525394,
                            -0.00566907,  0.66797545,  0.43135092,
                            -0.20196189,  0.04922572, -0.32425703 };

    variableVector gradChip = { 0.66070743, 0.53000315, 0.5279699 , 0.6826509 , 0.6746987 ,
                                0.38417229, 0.6754847 , 0.57437005, 0.50388402, 0.29566708,
                                0.41339237, 0.48189968, 0.29916895, 0.28835666, 0.76680497,
                                0.68779831, 0.44623234, 0.8414494 , 0.72621575, 0.50702875,
                                0.23052973, 0.35436878, 0.89915565, 0.87194386, 0.22637193,
                                0.97189461, 0.15759071 };

    variableVector answerFe = { -3.94817649, -0.32662858, -0.48781965,
                                -0.26623179,  1.60410514, -4.31494123,
                                 0.39966071, -0.05820434, -2.44387444 };

    variableVector answerChie = { -2.61886025, -0.72602181,  5.13908937,
                                   1.8405542 , -1.3016982 , -2.42597922,
                                  -2.61135485, -2.25166639,  0.89364486 };

    variableVector answerGradChie = { 1.14852526e+00,  4.82577470e-01,  3.17687002e-01,  3.15279224e+00,
                                     -3.98511181e+00, -5.72845791e+00, -8.03005606e+00,  4.03942034e+00,
                                     -1.19752568e+01, -4.58373637e+00,  4.69392821e-01, -2.55663946e-01,
                                     -2.90717158e-01,  2.96657140e+00,  2.65874647e+00, -8.50814019e+00,
                                     -4.27501281e+00,  8.35453206e-03, -6.64635005e+00, -3.65333789e+00,
                                     -2.08699660e+00,  3.09942546e+00,  7.38196309e-01,  5.96234908e-01,
                                     -8.64387705e+00, -7.45765574e-01, -8.42198541e+00 };

    variableVector resultFe, resultChie, resultGradChie;

    errorOut error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi, 
                                                                                    Fp, chip, gradChip,
                                                                                    resultFe, resultChie,
                                                                                    resultGradChie );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultFe, answerFe ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultChie, answerChie ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultGradChie, answerGradChie ) );

    variableVector resultFe2, resultChie2, resultGradChie2, resultInvFp, resultInvChip;

    error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi, 
                                                                            Fp, chip, gradChip,
                                                                            resultInvFp, resultInvChip,
                                                                            resultFe2, resultChie2,
                                                                            resultGradChie2 );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultFe2, answerFe ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultChie2, answerChie ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultGradChie2, answerGradChie ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultInvFp, vectorTools::inverse( Fp, 3, 3 ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultInvChip, vectorTools::inverse( chip, 3, 3 ) ) );

    //Test the computation of the jacobians

    variableVector resultFeJ, resultChieJ, resultGradChieJ;
    variableMatrix dFedF, dFedFp, dChiedChi, dChiedChip, dGradChiedFp, 
                   dGradChiedGradChi, dGradChiedGradChip, dGradChiedChi, dGradChiedChip;

    error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi, 
                                                                            Fp, chip, gradChip,
                                                                            resultFeJ, resultChieJ,
                                                                            resultGradChieJ,
                                                                            dFedF, dFedFp, dChiedChi,
                                                                            dChiedChip, dGradChiedGradChi,
                                                                            dGradChiedGradChip, dGradChiedFp,
                                                                            dGradChiedChi, dGradChiedChip );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultFeJ, answerFe ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultChieJ, answerChie ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultGradChieJ, answerGradChie ) );

    //Test the Jacobians w.r.t. the deformation gradient
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < F.size(); i++ ){
        constantVector delta( F.size(), 0 );
        delta[i] = eps * fabs( F[ i ] ) + eps;

        variableVector resultFeP,       resultFeM;
        variableVector resultChieP,     resultChieM;
        variableVector resultGradChieP, resultGradChieM;

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F + delta,  chi,  gradChi,
                                                                                Fp, chip, gradChip,
                                                                                resultFeP, resultChieP,
                                                                                resultGradChieP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F - delta,  chi,  gradChi,
                                                                                Fp, chip, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        BOOST_CHECK( !error );

        //Test dFedF
        variableVector gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dFedF[ j ][ i ] ) );
        }

        //Test dChiedF
        gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        //Test dGradChiedF
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );
    }

    //Test the Jacobians w.r.t. the plastic deformation gradient
    for ( unsigned int i = 0; i < Fp.size(); i++ ){
        constantVector delta( Fp.size(), 0 );
        delta[i] = eps * fabs( Fp[ i ] ) + eps;

        variableVector resultFeP,       resultFeM;
        variableVector resultChieP,     resultChieM;
        variableVector resultGradChieP, resultGradChieM;

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi,
                                                                                Fp + delta, chip, gradChip,
                                                                                resultFeP, resultChieP,
                                                                                resultGradChieP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi,
                                                                                Fp - delta, chip, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dFedFp[ j ][ i ] ) );
        }

        //Test dChiedFp
        gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        //Test dGradChiedFp
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dGradChiedFp[j][i] ) );
        }
    }

    //Test the Jacobians w.r.t. the micro deformation
    for ( unsigned int i = 0; i < chi.size(); i++ ){
        constantVector delta( chi.size(), 0 );
        delta[i] = eps * fabs( chi[ i ] ) + eps;

        variableVector resultFeP,       resultFeM;
        variableVector resultChieP,     resultChieM;
        variableVector resultGradChieP, resultGradChieM;

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi + delta,  gradChi,
                                                                                Fp, chip, gradChip,
                                                                                resultFeP, resultChieP,
                                                                                resultGradChieP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi - delta,  gradChi,
                                                                                Fp, chip, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        BOOST_CHECK( !error );

        //Test dChiedChi
        variableVector gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dChiedChi[ j ][ i ] ) );
        }

        //Test dFedChi
        gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        //Test dGradChiedChi
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dGradChiedChi[j][i] ) );
        }
    }

    //Test the Jacobians w.r.t. the plastic micro-deformation
    for ( unsigned int i = 0; i < chip.size(); i++ ){
        constantVector delta( chip.size(), 0 );
        delta[i] = eps * fabs( chip[ i ] ) + eps;

        variableVector resultFeP,       resultFeM;
        variableVector resultChieP,     resultChieM;
        variableVector resultGradChieP, resultGradChieM;

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi,
                                                                                Fp, chip + delta, gradChip,
                                                                                resultFeP, resultChieP,
                                                                                resultGradChieP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi,
                                                                                Fp, chip - delta, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dChiedChip[ j ][ i ] ) );
        }

        //Test dFedFp
        gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        //Test dGradChiedChip
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dGradChiedChip[j][i] ) );
        }
    }

    //Test the Jacobians w.r.t. the gradient of chi
    for ( unsigned int i = 0; i < gradChi.size(); i++ ){
        constantVector delta( gradChi.size(), 0 );
        delta[i] = eps * fabs( gradChi[ i ] ) + eps;

        variableVector resultFeP,       resultFeM;
        variableVector resultChieP,     resultChieM;
        variableVector resultGradChieP, resultGradChieM;

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi + delta,
                                                                                Fp, chip, gradChip,
                                                                                resultFeP, resultChieP,
                                                                                resultGradChieP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi - delta,
                                                                                Fp, chip, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        BOOST_CHECK( !error );

        //Test dFedGradChi
        variableVector gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        //Test dChiedGradChi
        gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        //Test dGradChiedGradChi
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dGradChiedGradChi[j][i] ) );
        }
    }

    //Test the Jacobians w.r.t. the plastic part of the gradient of chi
    for ( unsigned int i = 0; i < gradChi.size(); i++ ){
        constantVector delta( gradChi.size(), 0 );
        delta[i] = eps * fabs( gradChi[ i ] ) + eps;

        variableVector resultFeP,       resultFeM;
        variableVector resultChieP,     resultChieM;
        variableVector resultGradChieP, resultGradChieM;

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi,
                                                                                Fp, chip, gradChip + delta,
                                                                                resultFeP, resultChieP,
                                                                                resultGradChieP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi,
                                                                                Fp, chip, gradChip - delta,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        BOOST_CHECK( !error );

        //Test dFedGradChi
        variableVector gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        //Test dChiedGradChi
        gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        //Test dGradChiedGradChi
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dGradChiedGradChip[j][i] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputeElasticDeformationMeasures ){
    /*!
     * Test the computation of the elastic deformation measures 
     * required for the computation of the plasticity.
     *
     */

    variableVector Fe = { -3.94817649, -0.32662858, -0.48781965,
                          -0.26623179,  1.60410514, -4.31494123,
                           0.39966071, -0.05820434, -2.44387444 };

    variableVector chie = { -2.61886025, -0.72602181,  5.13908937,
                             1.8405542 , -1.3016982 , -2.42597922,
                            -2.61135485, -2.25166639,  0.89364486 };

    variableVector gradChie = { -1.27940191, -0.15466701, -0.38184898, -0.30173923,  0.05962743,
                                 0.88277412, -1.76241432, -0.60016475, -0.07033724, -0.98094822,
                                 0.57398185, -1.77481299, -0.10849973,  0.9656491 , -0.7146935 ,
                                -2.16271183, -2.03566316,  0.15276376, -1.26613437, -0.94886439,
                                -0.82649425,  0.02632722, -0.09189365,  0.56763039, -1.63935109,
                                 0.3039677 , -0.48933706 };

    variableVector answerCe = { 15.81870565,  0.8392615 ,  2.09805203,
                                 0.8392615 ,  2.68322729, -6.62003948,
                                 2.09805203, -6.62003948, 24.82920808 };

    variableVector answerCChie = { 17.06524293,   5.38540351, -20.25732698,
                                    5.38540351,   7.29152741,  -2.58538828,
                                  -20.25732698,  -2.58538828,  33.09421588 };

    variableVector answerPsie = { 8.80605252,   2.3131131 , -19.28700431,
                                  3.95982925,  -1.71986455,  -5.62211322,
                                 -0.28252834,  11.47370888,   5.77705312 };

    variableVector answerGammae = { 4.80644001,  0.07861661,  1.64980155,  1.23072776, -0.5292324 ,
                                   -3.06821432,  6.87892124,  3.03299854,  0.04146446, -1.08196034,
                                    1.02647393, -2.6741583 , -0.07702067,  1.53487528, -1.46782133,
                                   -2.79814493, -3.08707902,  0.29650483,  7.95112472, -0.0823429 ,
                                    9.86433536,  0.55102384, -3.97123001,  1.26600849, 14.19808301,
                                    8.33368016,  0.57102355 };

    variableVector resultCe, resultCChie, resultPsie, resultGammae;

    errorOut error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie, gradChie, resultCe,
                                                                                      resultCChie, resultPsie,
                                                                                      resultGammae );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultCe, answerCe ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultCChie, answerCChie ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultPsie, answerPsie ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultGammae, answerGammae ) );

    //Tests of Jacobians
    variableVector resultCeJ, resultCChieJ, resultPsieJ, resultGammaeJ;
    variableMatrix dRCGedFe, dMicroRCGedChie, dPsiedFe, dPsiedChie, dGammaedFe, dGammaedGradChie;

    error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie, gradChie, resultCeJ,
                                                                             resultCChieJ, resultPsieJ,
                                                                             resultGammaeJ, dRCGedFe,
                                                                             dMicroRCGedChie, dPsiedFe, dPsiedChie,
                                                                             dGammaedFe, dGammaedGradChie );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultCeJ, answerCe ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultCChieJ, answerCChie ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultPsieJ, answerPsie ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultGammaeJ, answerGammae ) );

    //Test jacobians w.r.t. the elastic deformation meaure
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < Fe.size(); i++ ){
        constantVector delta( Fe.size(), 0 );
        delta[i] = eps * fabs( Fe[ i ] ) + eps;

        variableVector resultCeP, resultCChieP, resultPsieP, resultGammaeP;
        variableVector resultCeM, resultCChieM, resultPsieM, resultGammaeM;

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe + delta, chie, gradChie,
                                                                                 resultCeP, resultCChieP, resultPsieP,
                                                                                 resultGammaeP );
        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe - delta, chie, gradChie,
                                                                                 resultCeM, resultCChieM, resultPsieM,
                                                                                 resultGammaeM );
        BOOST_CHECK( !error );

        variableVector gradCol = ( resultCeP - resultCeM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dRCGedFe[j][i] ) );
        }

        gradCol = ( resultCChieP - resultCChieM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultPsieP - resultPsieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dPsiedFe[j][i] ) );
        }

        gradCol = ( resultGammaeP - resultGammaeM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dGammaedFe[j][i] ) );
        }
    }

    //Test jacobians w.r.t. the elastic micro-deformation
    for ( unsigned int i = 0; i < chie.size(); i++ ){
        constantVector delta( chie.size(), 0 );
        delta[i] = eps * fabs( chie[ i ] ) + eps;

        variableVector resultCeP, resultCChieP, resultPsieP, resultGammaeP;
        variableVector resultCeM, resultCChieM, resultPsieM, resultGammaeM;

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie + delta, gradChie,
                                                                                 resultCeP, resultCChieP, resultPsieP,
                                                                                 resultGammaeP );
        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie - delta, gradChie,
                                                                                 resultCeM, resultCChieM, resultPsieM,
                                                                                 resultGammaeM );
        BOOST_CHECK( !error );

        variableVector gradCol = ( resultCeP - resultCeM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultCChieP - resultCChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dMicroRCGedChie[j][i] ) );
        }

        gradCol = ( resultPsieP - resultPsieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dPsiedChie[j][i] ) );
        }

        gradCol = ( resultGammaeP - resultGammaeM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );
    }

    for ( unsigned int i = 0; i < gradChie.size(); i++ ){
        constantVector delta( gradChie.size(), 0 );
        delta[i] = eps * fabs( gradChie[ i ] ) + eps;

        variableVector resultCeP, resultCChieP, resultPsieP, resultGammaeP;
        variableVector resultCeM, resultCChieM, resultPsieM, resultGammaeM;

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie, gradChie + delta,
                                                                                 resultCeP, resultCChieP, resultPsieP,
                                                                                 resultGammaeP );
        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie, gradChie - delta,
                                                                                 resultCeM, resultCChieM, resultPsieM,
                                                                                 resultGammaeM );
        BOOST_CHECK( !error );

        variableVector gradCol = ( resultCeP - resultCeM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultCChieP - resultCChieM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultPsieP - resultPsieM ) / ( 2 * delta[i] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultGammaeP - resultGammaeM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dGammaedGradChie[j][i] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticVelocityGradients ){
    /*!
     * Test the computation of the plastic velocity gradients.
     *
     */

    variableType macroGamma = 0.08166694603978908;
    variableType microGamma = 0.8652174130049269;
    variableVector microGradientGamma = { 0.17245016, 0.92420274, 0.28114459 };

    variableVector Ce = { 15.81870565,  0.8392615 ,  2.09805203,
                           0.8392615 ,  2.68322729, -6.62003948,
                           2.09805203, -6.62003948, 24.82920808 };

    variableVector microCe = { 17.06524293,   5.38540351, -20.25732698,
                                5.38540351,   7.29152741,  -2.58538828,
                              -20.25732698,  -2.58538828,  33.09421588 };

    variableVector Psie = {  8.80605252,   2.3131131 , -19.28700431,
                             3.95982925,  -1.71986455,  -5.62211322,
                            -0.28252834,  11.47370888,   5.77705312 };

    variableVector Gammae = { 4.80644001,  0.07861661,  1.64980155,  1.23072776, -0.5292324 ,
                             -3.06821432,  6.87892124,  3.03299854,  0.04146446, -1.08196034,
                              1.02647393, -2.6741583 , -0.07702067,  1.53487528, -1.46782133,
                             -2.79814493, -3.08707902,  0.29650483,  7.95112472, -0.0823429 ,
                              9.86433536,  0.55102384, -3.97123001,  1.26600849, 14.19808301,
                              8.33368016,  0.57102355 };

    variableVector macroFlowDirection = { 0.78884638, 0.19639211, 0.15523073,
                                          0.47307595, 0.28241451, 0.66404732,
                                          0.1634089 , 0.92452471, 0.77390742 };

    variableVector microFlowDirection = { 0.86300151, 0.95736394, 0.61329255,
                                          0.27780339, 0.26054793, 0.33313753,
                                          0.34289169, 0.57971261, 0.51536929 };

    variableVector microGradientFlowDirection = { 0.38320117, 0.00147635, 0.22526135, 0.24857347, 0.44904944,
                                                  0.39175461, 0.94088825, 0.04088633, 0.95042374, 0.44676197,
                                                  0.33100061, 0.79806506, 0.05883935, 0.20924962, 0.83681153,
                                                  0.12116776, 0.39737069, 0.07417313, 0.5859491 , 0.28899583,
                                                  0.91967175, 0.413024  , 0.97723212, 0.81694258, 0.92037483,
                                                  0.84857389, 0.74623422, 0.65442987, 0.15706966, 0.03580793,
                                                  0.98977654, 0.6414159 , 0.03345668, 0.73436727, 0.25417675,
                                                  0.594925  , 0.4871345 , 0.27395216, 0.23644903, 0.42902409,
                                                  0.24760169, 0.16207352, 0.68475097, 0.86214768, 0.9734798 ,
                                                  0.86141159, 0.98250926, 0.25056881, 0.8315578 , 0.95970017,
                                                  0.62180382, 0.52207192, 0.66811873, 0.06083854, 0.59855098,
                                                  0.41784728, 0.41193658, 0.3161969 , 0.75697096, 0.2172361 ,
                                                  0.5170385 , 0.52482239, 0.55849978, 0.60039656, 0.38358062,
                                                  0.66214191, 0.22829067, 0.10781315, 0.40531347, 0.25340843,
                                                  0.89016033, 0.85477638, 0.43630125, 0.35699992, 0.3784267 ,
                                                  0.12262464, 0.38383612, 0.12695384, 0.74207569, 0.58531619,
                                                  0.08294492 };

    variableVector answerLp = { -0.05489573, -0.01980382, -0.06060589,
                                 1.1610081 ,  0.4002548 ,  0.86866858,
                                 0.33607202,  0.12218348,  0.25723268 };

    variableVector answerMicroLp = { -1.14307645,  1.41916637,  2.9980591 ,
                                      0.03471568, -0.06540477, -0.10626039,
                                     -0.42684692,  0.51754693,  1.10823355 };

    variableVector answerMicroGradientLp = {  0.47438161,  0.25861926,  1.1761207 ,  1.36722161, -0.47991786,
                                              1.87711396,  3.53355655,  1.0243844 ,  1.95422675,  1.55466612,
                                              0.34239846,  1.05904209, -2.01539758, -0.30114924, -1.40333162,
                                             -4.18953464, -0.80091983, -2.94449224, -0.13359547,  0.05591392,
                                              0.39411415,  0.90630249, -0.13965495,  0.72105693,  2.11543948,
                                              0.44788947,  0.783895 }; 

    variableVector resultLp, resultMicroLp, resultMicroGradientLp;

    errorOut error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                                    Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                                    microFlowDirection, microGradientFlowDirection,
                                                                                    resultLp, resultMicroLp, resultMicroGradientLp );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultLp, answerLp ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroLp, answerMicroLp ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroGradientLp, answerMicroGradientLp ) );

    //Tests of the Jacobians
    variableVector resultLpJ, resultMicroLpJ, resultMicroGradientLpJ;

    variableVector dMacroLpdMacroGammaJ, dMacroLpdMicroGammaJ;
    variableVector dMicroLpdMicroGammaJ;
    variableVector dMicroGradientLpdMicroGammaJ;
    variableMatrix dMicroGradientLpdMicroGradientGammaJ;

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpJ, resultMicroLpJ, resultMicroGradientLpJ,
                                                                           dMacroLpdMacroGammaJ, dMacroLpdMicroGammaJ,
                                                                           dMicroLpdMicroGammaJ,
                                                                           dMicroGradientLpdMicroGammaJ,
                                                                           dMicroGradientLpdMicroGradientGammaJ );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultLpJ, answerLp ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroLpJ, answerMicroLp ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroGradientLpJ, answerMicroGradientLp ) );

    variableVector resultLpJ2, resultMicroLpJ2, resultMicroGradientLpJ2;

    variableVector dMacroLpdMacroGammaJ2, dMacroLpdMicroGammaJ2;
    variableVector dMicroLpdMicroGammaJ2;
    variableVector dMicroGradientLpdMicroGammaJ2;
    variableMatrix dMicroGradientLpdMicroGradientGammaJ2;
    variableMatrix dPlasticMacroLdElasticRCGJ2,
                   dPlasticMacroLdMacroFlowDirectionJ2,
                   dPlasticMacroLdMicroFlowDirectionJ2,
                   dPlasticMicroLdMicroElasticRCGJ2,
                   dPlasticMicroLdElasticPsiJ2,
                   dPlasticMicroLdMicroFlowDirectionJ2,
                   dPlasticMicroGradientLdMicroElasticRCGJ2,
                   dPlasticMicroGradientLdElasticPsiJ2,
                   dPlasticMicroGradientLdElasticGammaJ2,
                   dPlasticMicroGradientLdMicroFlowDirectionJ2,
                   dPlasticMicroGradientLdMicroGradientFlowDirectionJ2;

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpJ2, resultMicroLpJ2, resultMicroGradientLpJ2,
                                                                           dMacroLpdMacroGammaJ2, dMacroLpdMicroGammaJ2,
                                                                           dMicroLpdMicroGammaJ2,
                                                                           dMicroGradientLpdMicroGammaJ2,
                                                                           dMicroGradientLpdMicroGradientGammaJ2,
                                                                           dPlasticMacroLdElasticRCGJ2,
                                                                           dPlasticMacroLdMacroFlowDirectionJ2,
                                                                           dPlasticMacroLdMicroFlowDirectionJ2,
                                                                           dPlasticMicroLdMicroElasticRCGJ2,
                                                                           dPlasticMicroLdElasticPsiJ2,
                                                                           dPlasticMicroLdMicroFlowDirectionJ2,
                                                                           dPlasticMicroGradientLdMicroElasticRCGJ2,
                                                                           dPlasticMicroGradientLdElasticPsiJ2,
                                                                           dPlasticMicroGradientLdElasticGammaJ2,
                                                                           dPlasticMicroGradientLdMicroFlowDirectionJ2,
                                                                           dPlasticMicroGradientLdMicroGradientFlowDirectionJ2 );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultLpJ2, answerLp ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroLpJ2, answerMicroLp ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroGradientLpJ2, answerMicroGradientLp ) );

    //Tests of Jacobians w.r.t. macroGamma
    constantType eps = 1e-6;
    constantType scalarDelta = eps * fabs( macroGamma ) + eps;

    variableVector resultLpP, resultMicroLpP, resultMicroGradientLpP;
    variableVector resultLpM, resultMicroLpM, resultMicroGradientLpM;

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma + scalarDelta, microGamma, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpP, resultMicroLpP, resultMicroGradientLpP );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma - scalarDelta, microGamma, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpM, resultMicroLpM, resultMicroGradientLpM );

    BOOST_CHECK( !error );

    variableVector gradCol = ( resultLpP - resultLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMacroLpdMacroGammaJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMacroLpdMacroGammaJ2 ) );

    gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

    gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

    //Test of Jacobians w.r.t. microGamma
    scalarDelta = eps * fabs( microGamma ) + eps;

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma + scalarDelta, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpP, resultMicroLpP, resultMicroGradientLpP );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma - scalarDelta, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpM, resultMicroLpM, resultMicroGradientLpM );

    BOOST_CHECK( !error );

    gradCol = ( resultLpP - resultLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol,  dMacroLpdMicroGammaJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol,  dMacroLpdMicroGammaJ2 ) );

    gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGammaJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGammaJ2 ) );

    gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol,  dMicroGradientLpdMicroGammaJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol,  dMicroGradientLpdMicroGammaJ2 ) );

    //Test of Jacobians w.r.t. the microGradientGamma
    for ( unsigned int i = 0; i < microGradientGamma.size(); i++ ){
        constantVector delta( microGradientGamma.size(), 0 );
        delta[i] = eps * fabs( microGradientGamma[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma + delta,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma - delta,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientLpdMicroGradientGammaJ[ j ][ i ] ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientLpdMicroGradientGammaJ2[ j ][ i ] ) );
        }
    }

    //Test of Jacobians w.r.t. the macro right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < Ce.size(); i++ ){
        constantVector delta( Ce.size(), 0 );
        delta[i] = eps * fabs( Ce[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce + delta, microCe, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce - delta, microCe, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMacroLdElasticRCGJ2[ j ][ i ] ) );
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );
    }

    //Test of Jacobians w.r.t. the micro right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < microCe.size(); i++ ){
        constantVector delta( microCe.size(), 0 );
        delta[i] = eps * fabs( microCe[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe + delta, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe - delta, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroLdMicroElasticRCGJ2[ j ][ i ] ) );
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroElasticRCGJ2[ j ][ i ] ) );
        }
    }

    //Test of Jacobians w.r.t. the micro deformation tensor Psi
    for ( unsigned int i = 0; i < Psie.size(); i++ ){
        constantVector delta( Psie.size(), 0 );
        delta[i] = eps * fabs( Psie[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie + delta, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie - delta, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroLdElasticPsiJ2[ j ][ i ] ) );
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticPsiJ2[ j ][ i ] ) );
        }
    }

    //Test of Jacobians w.r.t. the higher order deformation tensor
    for ( unsigned int i = 0; i < Gammae.size(); i++ ){
        constantVector delta( Gammae.size(), 0 );
        delta[i] = eps * fabs( Gammae[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae + delta, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae - delta, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticGammaJ2[ j ][ i ] ) );
        }
    }

    //Test of Jacobians w.r.t. the macro flow direction
    for ( unsigned int i = 0; i < macroFlowDirection.size(); i++ ){
        constantVector delta( macroFlowDirection.size(), 0 );
        delta[i] = eps * fabs( macroFlowDirection[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection + delta, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection - delta, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMacroLdMacroFlowDirectionJ2[ j ][ i ] ) );
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    //Test of Jacobians w.r.t. the micro flow direction
    for ( unsigned int i = 0; i < microFlowDirection.size(); i++ ){
        constantVector delta( microFlowDirection.size(), 0 );
        delta[i] = eps * fabs( microFlowDirection[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection,
                                                                               microFlowDirection + delta, microGradientFlowDirection,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection,
                                                                               microFlowDirection - delta, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMacroLdMicroFlowDirectionJ2[ j ][ i ] ) );
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroLdMicroFlowDirectionJ2[ j ][ i ] ) );
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroFlowDirectionJ2[ j ][ i ] ) );
        }
    }

    //Test of Jacobians w.r.t. the micro gradient flow direction
    for ( unsigned int i = 0; i < microGradientFlowDirection.size(); i++ ){
        constantVector delta( microGradientFlowDirection.size(), 0 );
        delta[ i ] = eps * fabs( microGradientFlowDirection[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection,
                                                                               microFlowDirection, microGradientFlowDirection + delta,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection,
                                                                               microFlowDirection, microGradientFlowDirection - delta,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientFlowDirectionJ2[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticMacroVelocityGradient ){
    /*!
     * Test the computation of the plastic macro velocity gradient.
     *
     */

    variableType macroGamma = 0.08166694603978908;
    variableType microGamma = 0.8652174130049269;

    variableVector Ce = { 15.81870565,  0.8392615 ,  2.09805203,
                           0.8392615 ,  2.68322729, -6.62003948,
                           2.09805203, -6.62003948, 24.82920808 };

    variableVector inverseCe = vectorTools::inverse( Ce, 3, 3 );

    variableVector macroFlowDirection = { 0.78884638, 0.19639211, 0.15523073,
                                          0.47307595, 0.28241451, 0.66404732,
                                          0.1634089 , 0.92452471, 0.77390742 };

    variableVector microFlowDirection = { 0.86300151, 0.95736394, 0.61329255,
                                          0.27780339, 0.26054793, 0.33313753,
                                          0.34289169, 0.57971261, 0.51536929 };

    variableVector answerMacroLp = { -0.05489573, -0.01980382, -0.06060589,
                                      1.1610081 ,  0.4002548 ,  0.86866858,
                                      0.33607202,  0.12218348,  0.25723268 };

    variableVector resultMacroLp;

    errorOut error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                                        macroFlowDirection, microFlowDirection,
                                                                                        resultMacroLp );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMacroLp, resultMacroLp ) );

    //Tests of the Jacobians
    variableVector resultMacroLpJ;
    variableVector dMacroLdMacroGammaJ, dMacroLdMicroGammaJ;

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                               macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpJ, dMacroLdMacroGammaJ,
                                                                               dMacroLdMicroGammaJ );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMacroLp, resultMacroLpJ ) );

    variableVector resultMacroLpJ2;
    variableVector dMacroLdMacroGammaJ2, dMacroLdMicroGammaJ2;
    variableMatrix dMacroLdElasticRCG, dMacroLdMacroFlowDirection, dMacroLdMicroFlowDirection;

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                               macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpJ2, dMacroLdMacroGammaJ2,
                                                                               dMacroLdMicroGammaJ2, dMacroLdElasticRCG,
                                                                               dMacroLdMacroFlowDirection,
                                                                               dMacroLdMicroFlowDirection );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMacroLp, resultMacroLpJ2 ) );

    //Tests Jacobians w.r.t. macroGamma
    constantType eps = 1e-6;
    constantType scalarDelta = eps * fabs( macroGamma) + eps;

    variableVector resultMacroLpP, resultMacroLpM;

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma + scalarDelta, microGamma,
                                                                               inverseCe, macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpP );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma - scalarDelta, microGamma,
                                                                               inverseCe, macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpM );

    BOOST_CHECK( !error );

    variableVector gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMacroLdMacroGammaJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMacroLdMacroGammaJ2 ) );

    //Test Jacobians w.r.t. microGamma
    scalarDelta = eps * fabs( microGamma) + eps;

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma + scalarDelta,
                                                                               inverseCe, macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpP );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma - scalarDelta,
                                                                               inverseCe, macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpM );

    BOOST_CHECK( !error );

    gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMacroLdMicroGammaJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMacroLdMicroGammaJ2 ) );

    //Test Jacobians w.r.t. the right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < Ce.size(); i++ ){
        constantVector delta( Ce.size(), 0 );
        delta[i] = eps * fabs( Ce[i] ) + eps;

        variableVector inverseCeTemp = vectorTools::inverse( Ce + delta, 3, 3 );

        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCeTemp, macroFlowDirection,
                                                                                   microFlowDirection, resultMacroLpP );
    
        BOOST_CHECK( !error );

        inverseCeTemp = vectorTools::inverse( Ce - delta, 3, 3 );
    
        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCeTemp, macroFlowDirection,
                                                                                   microFlowDirection, resultMacroLpM );
    
        BOOST_CHECK( !error );
    
        gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMacroLdElasticRCG[ j ][ i ] ) );
        }
    }

    //Test Jacobians w.r.t. the macro flow direction
    for ( unsigned int i = 0; i < macroFlowDirection.size(); i++ ){
        constantVector delta( macroFlowDirection.size(), 0 );
        delta[i] = eps * fabs( macroFlowDirection[i] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCe, macroFlowDirection + delta,
                                                                                   microFlowDirection, resultMacroLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCe, macroFlowDirection - delta,
                                                                                   microFlowDirection, resultMacroLpM );
    
        BOOST_CHECK( !error );
    
        gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMacroLdMacroFlowDirection[ j ][ i ] ) );
        }
    }

    //Test Jacobians w.r.t. the micro flow direction
    for ( unsigned int i = 0; i < microFlowDirection.size(); i++ ){
        constantVector delta( microFlowDirection.size(), 0 );
        delta[i] = eps * fabs( microFlowDirection[i] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCe, macroFlowDirection,
                                                                                   microFlowDirection + delta, resultMacroLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCe, macroFlowDirection,
                                                                                   microFlowDirection - delta, resultMacroLpM );
    
        BOOST_CHECK( !error );
    
        gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMacroLdMicroFlowDirection[ j ][ i ] ) );
        }
    }

}


BOOST_AUTO_TEST_CASE( testComputePlasticMicroVelocityGradient ){
    /*!
     * Test the computation of the plastic micro velocity gradient.
     *
     */

    variableType microGamma = 0.8652174130049269;

    variableVector Ce = { 17.06524293,   5.38540351, -20.25732698,
                           5.38540351,   7.29152741,  -2.58538828,
                         -20.25732698,  -2.58538828,  33.09421588 };

    variableVector Psie = { 8.80605252,   2.3131131 , -19.28700431,
                            3.95982925,  -1.71986455,  -5.62211322,
                           -0.28252834,  11.47370888,   5.77705312 };

    variableVector invPsie = vectorTools::inverse( Psie, 3, 3 );

    variableVector microFlowDirection = { 0.86300151, 0.95736394, 0.61329255,
                                          0.27780339, 0.26054793, 0.33313753,
                                          0.34289169, 0.57971261, 0.51536929 };

    variableVector answerMicroLp = { -1.14307645,  1.41916637,  2.9980591 ,
                                      0.03471568, -0.06540477, -0.10626039,
                                     -0.42684692,  0.51754693,  1.10823355 };

    variableVector resultMicroLp;

    errorOut error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                        microFlowDirection, resultMicroLp );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroLp, resultMicroLp ) );

    //Test the Jacobians
    variableVector resultMicroLpJ;
    variableVector dMicroLpdMicroGamma;

    error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                               microFlowDirection, resultMicroLpJ,
                                                                               dMicroLpdMicroGamma );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroLp, resultMicroLpJ ) );

    variableVector resultMicroLpJ2;
    variableVector dMicroLpdMicroGamma2;
    variableMatrix dMicroLpdMicroRCG, dMicroLpdPsie, dMicroLpdMicroFlowDirection;

    error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                               microFlowDirection, resultMicroLpJ2,
                                                                               dMicroLpdMicroGamma2,
                                                                               dMicroLpdMicroRCG, dMicroLpdPsie,
                                                                               dMicroLpdMicroFlowDirection );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroLp, resultMicroLpJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroLp, resultMicroLpJ2 ) );

    constantType eps = 1e-6;
    constantType scalarDelta = eps * fabs( microGamma ) + eps;

    variableVector resultMicroLpP, resultMicroLpM;

    error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma + scalarDelta, Ce, Psie, invPsie,
                                                                               microFlowDirection, resultMicroLpP );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma - scalarDelta, Ce, Psie, invPsie,
                                                                               microFlowDirection, resultMicroLpM );

    BOOST_CHECK( !error );

    variableVector gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGamma ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGamma2 ) );

    //Test Jacobian w.r.t. the elastic micro right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < Ce.size(); i++ ){
        constantVector delta( Ce.size(), 0 );
        delta[i] = eps * fabs( Ce[i] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce + delta, Psie, invPsie,
                                                                                   microFlowDirection, resultMicroLpP );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce - delta, Psie, invPsie,
                                                                                   microFlowDirection, resultMicroLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dMicroLpdMicroRCG[ j ][ i ] ) );
        }
    }

    //Test Jacobian w.r.t. the elastic micro deformation measure Psi
    for ( unsigned int i = 0; i < Psie.size(); i++ ){
        constantVector delta( Psie.size(), 0 );
        delta[i] = eps * fabs( Psie[i] ) + eps;

        variableVector PsieTemp = Psie + delta;
        variableVector invPsieTemp = vectorTools::inverse( PsieTemp, 3, 3 );

        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, PsieTemp, invPsieTemp,
                                                                                   microFlowDirection, resultMicroLpP );
    
        BOOST_CHECK( !error );
    
        PsieTemp = Psie - delta;
        invPsieTemp = vectorTools::inverse( PsieTemp, 3, 3 );

        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, PsieTemp, invPsieTemp,
                                                                                   microFlowDirection, resultMicroLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dMicroLpdPsie[ j ][ i ] ) );
        }
    }

    for ( unsigned int i = 0; i < microFlowDirection.size(); i++ ){
        constantVector delta( microFlowDirection.size(), 0 );
        delta[i] = eps * fabs( microFlowDirection[i] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                   microFlowDirection + delta, resultMicroLpP );
    
        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                   microFlowDirection - delta, resultMicroLpM );
    
        BOOST_CHECK( !error );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[j], dMicroLpdMicroFlowDirection[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticMicroGradientVelocityGradient ){
    /*!
     * Test the computation of the plastic micro gradient velocity gradient.
     *
     */

    variableVector microGradientGamma = { 0.17245016, 0.92420274, 0.28114459 };

    variableVector Psie = { 8.80605252,   2.3131131 , -19.28700431,
                            3.95982925,  -1.71986455,  -5.62211322,
                           -0.28252834,  11.47370888,   5.77705312 };

    variableVector invPsie = vectorTools::inverse( Psie, 3, 3 );

    variableVector elasticGamma = { 4.80644001,  0.07861661,  1.64980155,  1.23072776, -0.5292324 ,
                                   -3.06821432,  6.87892124,  3.03299854,  0.04146446, -1.08196034,
                                    1.02647393, -2.6741583 , -0.07702067,  1.53487528, -1.46782133,
                                   -2.79814493, -3.08707902,  0.29650483,  7.95112472, -0.0823429 ,
                                    9.86433536,  0.55102384, -3.97123001,  1.26600849, 14.19808301,
                                    8.33368016,  0.57102355 };

    variableVector microGradientFlowDirection = { 0.38320117, 0.00147635, 0.22526135, 0.24857347, 0.44904944,
                                                  0.39175461, 0.94088825, 0.04088633, 0.95042374, 0.44676197,
                                                  0.33100061, 0.79806506, 0.05883935, 0.20924962, 0.83681153,
                                                  0.12116776, 0.39737069, 0.07417313, 0.5859491 , 0.28899583,
                                                  0.91967175, 0.413024  , 0.97723212, 0.81694258, 0.92037483,
                                                  0.84857389, 0.74623422, 0.65442987, 0.15706966, 0.03580793,
                                                  0.98977654, 0.6414159 , 0.03345668, 0.73436727, 0.25417675,
                                                  0.594925  , 0.4871345 , 0.27395216, 0.23644903, 0.42902409,
                                                  0.24760169, 0.16207352, 0.68475097, 0.86214768, 0.9734798 ,
                                                  0.86141159, 0.98250926, 0.25056881, 0.8315578 , 0.95970017,
                                                  0.62180382, 0.52207192, 0.66811873, 0.06083854, 0.59855098,
                                                  0.41784728, 0.41193658, 0.3161969 , 0.75697096, 0.2172361 ,
                                                  0.5170385 , 0.52482239, 0.55849978, 0.60039656, 0.38358062,
                                                  0.66214191, 0.22829067, 0.10781315, 0.40531347, 0.25340843,
                                                  0.89016033, 0.85477638, 0.43630125, 0.35699992, 0.3784267 ,
                                                  0.12262464, 0.38383612, 0.12695384, 0.74207569, 0.58531619,
                                                  0.08294492 };

    variableVector microLp = { -85.67983387, -16.91839826, 127.3318347 ,
                                 0.65035144,   0.1459189 ,  -0.71988301,
                               -36.05794838,  -7.86041652,  52.33737079 };

    variableVector answerMicroGradLp = { -83.50670143,  -11.14831106,  -39.09529065,  -16.58901227,
                                          -8.71869628,   19.76266969,   64.74989223,    3.55092183,
                                          59.30524632,  126.06867513,   26.68572242,   83.03661657,
                                          25.93560443,    6.13134257,   16.20263835, -185.220318  ,
                                         -38.17442863, -123.4797927 ,  -56.62958098,   -7.64862742,
                                         -16.00279814,  -11.26027965,   -4.19458173,    8.34519958,
                                          58.1455205 ,    5.42232797,   23.44933638 };

    variableVector answerSkewTerm = { -84.06988022,  -11.39160049,  -39.5285613 ,  -17.07306624,
                                       -8.94361366,   19.17011533,   64.59515902,    3.2768447 ,
                                       58.94068694,  126.07596462,   26.65609256,   83.02552879,
                                       25.99081849,    6.07289457,   16.21498441, -185.26655691,
                                      -38.22947657, -123.43856309,  -56.84233212,   -7.7271727 ,
                                      -16.14907393,  -11.46103628,   -4.28260582,    8.13100026,
                                       58.07906172,    5.31870592,   23.31357689 };

    variableVector resultMicroGradLp;

    errorOut error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie, 
                                                                                                elasticGamma, microGradientFlowDirection,
                                                                                                microLp, resultMicroGradLp );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLp ) );

    variableVector resultMicroGradLp1;
    variableVector resultSkewTerm;

    error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie, 
                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                       microLp, resultMicroGradLp1,
                                                                                       resultSkewTerm );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLp1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerSkewTerm, resultSkewTerm ) );

    //Test the Jacobians
    variableVector resultMicroGradLpJ;
    variableMatrix dPlasticMicroGradientLdMicroGradientGamma;
    variableMatrix dPlasticMicroGradientLdPlasticMicroL;

    error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                       microLp, resultMicroGradLpJ,
                                                                                       dPlasticMicroGradientLdMicroGradientGamma,
                                                                                       dPlasticMicroGradientLdPlasticMicroL );
    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLpJ ) );

    variableVector resultMicroGradLpJ2, resultSkewTerm2;
    variableMatrix dPlasticMicroGradientLdMicroGradientGamma2;
    variableMatrix dPlasticMicroGradientLdPlasticMicroL2;

    error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                       microLp, resultMicroGradLpJ2, resultSkewTerm2,
                                                                                       dPlasticMicroGradientLdMicroGradientGamma2,
                                                                                       dPlasticMicroGradientLdPlasticMicroL2 );
    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLpJ2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerSkewTerm, resultSkewTerm2 ) );

    variableVector resultMicroGradLpJ3;
    variableMatrix dPlasticMicroGradientLdMicroGradientGamma3;
    variableMatrix dPlasticMicroGradientLdPlasticMicroL3;
    variableMatrix dPlasticMicroGradientLdElasticPsi;
    variableMatrix dPlasticMicroGradientLdElasticGamma;
    variableMatrix dPlasticMicroGradientLdMicroGradientFlowDirection;

    error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                       microLp, resultMicroGradLpJ3,
                                                                                       dPlasticMicroGradientLdMicroGradientGamma3,
                                                                                       dPlasticMicroGradientLdPlasticMicroL3,
                                                                                       dPlasticMicroGradientLdElasticPsi,
                                                                                       dPlasticMicroGradientLdElasticGamma,
                                                                                       dPlasticMicroGradientLdMicroGradientFlowDirection );
    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLpJ3 ) );

    //Test computation of Jacobians w.r.t. microGradientGamma
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < microGradientGamma.size(); i++ ){
        constantVector delta( microGradientGamma.size(), 0 );
        delta[i] = eps * fabs( microGradientGamma[i] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma + delta, Psie, invPsie,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma - delta, Psie, invPsie,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientGamma[ j ][ i ] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientGamma2[ j ][ i ] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientGamma3[ j ][ i ] ) );
        }
    }

    //Test computation of Jacobians w.r.t. the plastic micro velocity gradient
    for ( unsigned int i = 0; i < microLp.size(); i++ ){
        constantVector delta( microLp.size(), 0 );
        delta[i] = eps * fabs( microLp[i] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp + delta, resultMicroGradLpP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp - delta, resultMicroGradLpM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdPlasticMicroL[ j ][ i ] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdPlasticMicroL2[ j ][ i ] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdPlasticMicroL3[ j ][ i ] ) );
        }
    }

    //Test computation of Jacobian w.r.t. the micro deformation measure Psi
    for ( unsigned int i = 0; i < Psie.size(); i++ ){
        constantVector delta( Psie.size(), 0 );
        delta[i] = eps * fabs( Psie[i] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        variableVector PsieTemp = Psie + delta;
        variableVector invPsieTemp = vectorTools::inverse( PsieTemp, 3, 3 );

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, PsieTemp, invPsieTemp,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpP );

        BOOST_CHECK( !error );

        PsieTemp = Psie - delta;
        invPsieTemp = vectorTools::inverse( PsieTemp, 3, 3 );

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, PsieTemp, invPsieTemp,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticPsi[ j ][ i ] ) );
        }
    }

    //Test computation of Jacobian w.r.t. the elastic higher order deformation metric Gamma
    for ( unsigned int i = 0; i < elasticGamma.size(); i++ ){
        constantVector delta( elasticGamma.size(), 0 );
        delta[i] = eps * fabs( elasticGamma[i] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                           elasticGamma + delta,
                                                                                           microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                           elasticGamma - delta,
                                                                                           microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticGamma[ j ][ i ] ) );
        }
    }

    //Test computation of Jacobian w.r.t. the elastic higher order deformation metric Gamma
    for ( unsigned int i = 0; i < microGradientFlowDirection.size(); i++ ){
        constantVector delta( microGradientFlowDirection.size(), 0 );
        delta[ i ] = eps * fabs( microGradientFlowDirection[ i ] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                           elasticGamma,
                                                                                           microGradientFlowDirection + delta,
                                                                                           microLp, resultMicroGradLpP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                           elasticGamma,
                                                                                           microGradientFlowDirection - delta,
                                                                                           microLp, resultMicroGradLpM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientFlowDirection[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testEvolvePlasticMicroGradChi ){
    /*!
     * Test the evolution of the plastic micro gradient deformation.
     *
     */

    variableType Dt = 7.888463751831797;

    variableVector currentPlasticMicroDeformation = { -0.49993961, -0.36238477, -0.24525394,
                                                      -0.00566907,  0.66797545,  0.43135092,
                                                      -0.20196189,  0.04922572, -0.32425703 };

    variableVector currentPlasticMacroVelocityGradient = { -0.05489573, -0.01980382, -0.06060589,
                                                            1.1610081 ,  0.4002548 ,  0.86866858,
                                                            0.33607202,  0.12218348,  0.25723268 };

    variableVector currentPlasticMicroVelocityGradient = { -21.41995847,  -4.22959956,  31.83295868,
                                                             0.16258786,   0.03647972,  -0.17997075,
                                                            -9.01448709,  -1.96510413,  13.0843427 };

    variableVector currentPlasticMicroGradientVelocityGradient = { -70.66715632,   28.05094192,  -12.32029375,    9.17810714,
                                                                   -19.71784064,  -49.77471891,   40.06703071,  -18.75257585,
                                                                   -87.64798388,  216.34003077,   11.38331058,   52.06095799,
                                                                    44.55771056,    2.81701824,   11.53413479, -317.46326014,
                                                                   -16.12786945,  -75.01963462,  -38.71407841,   26.25281724,
                                                                   -31.89136164,    2.12480823,   -4.79956514,  -26.12979314,
                                                                    27.25829277,  -30.0572412 ,    1.8995271 };

    variableVector previousInversePlasticMicroDeformation = { -0.62976501, -0.65543815,  0.07568244,
                                                               0.25069488,  0.04845195,  0.85682533,
                                                              -0.69408256,  0.68092037,  0.06724845 };

    variableVector previousPlasticMicroGradient = { 0.02280625, 0.22512007, 0.85575746, 0.1825644 , 0.97734329,
                                                    0.15265694, 0.0984977 , 0.73878709, 0.14943404, 0.50190449,
                                                    0.1394924 , 0.49979228, 0.69297073, 0.03931335, 0.61839353,
                                                    0.32686324, 0.67661231, 0.54568142, 0.60788606, 0.46255347,
                                                    0.94481369, 0.6944252 , 0.29844176, 0.55288832, 0.33853861,
                                                    0.03626558, 0.56699219 };

    variableVector previousPlasticMacroVelocityGradient = { 0.63999789, 0.35571483, 0.63708287,
                                                            1.82175239, 1.47320642, 2.26307223,
                                                            0.9862676 , 1.15883656, 1.36908512 };

    variableVector previousPlasticMicroVelocityGradient = { -0.10429558, -0.17299031, -0.13991628,
                                                            -0.1431904 ,  0.18911688,  0.20668886,
                                                             0.01882534,  0.20244689,  0.16117068 };

    variableVector previousPlasticMicroGradientVelocityGradient = { -1.28938531, -0.78106843, -0.16682568, -0.98659316, -0.33046864,
                                                                     1.00840067, -0.86124472,  0.16817579,  0.26768819,  2.42835457,
                                                                     2.54123143,  0.66070946,  2.88282966,  3.30378025,  0.12668134,
                                                                     3.52709363,  3.17568822,  0.4061731 ,  2.3619489 ,  1.92411292,
                                                                     0.39428284,  0.84845341,  0.29862915,  2.21825071,  1.58372838,
                                                                     0.01050663,  2.0526354 };

    parameterType alpha = 0.19639211333133877;

    variableVector answerCurrentPlasticMicroGradient = {  201.71721607,   -8.0663384 ,   63.93960976,  299.19611487,
                                                          -15.75827887,   98.74480552, -319.49770439,   13.5458325 ,
                                                         -102.33587287, -195.29766492,   -5.08537627,  -39.70978333,
                                                         -291.02098952,   -3.24570068,  -64.71587095,  314.13870087,
                                                            5.09285553,   66.08381791,  109.96443867,   -5.79910293,
                                                           37.22865464,  163.04783477,  -10.35151403,   59.19030845,
                                                         -175.10575052,   10.13754922,  -60.7404024 };

    variableMatrix answerLHS = {
        {1.3643808e+02, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0., 0., 0., 0., 0.},
        {-1.2554098e-01, 1.3932339e+02, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0., 0., 0., 0.},
        {-3.8419476e-01, 5.5066914e+00, 1.3841674e+02, 0., 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 1.3643808e+02, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0., 0.},
        {0., 0., 0., -1.2554098e-01, 1.3932339e+02, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0.},
        {0., 0., 0., -3.8419476e-01, 5.5066914e+00, 1.3841674e+02, 0., 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 1.3643808e+02, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0.},
        {0., 0., 0., 0., 0., 0., -1.2554098e-01, 1.3932339e+02, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0.},
        {0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, 1.3841674e+02, 0., 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02},
        {-1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., 0., 4.2074983e-01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, 3.3060545e+00, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, 2.3994041e+00, 0., 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., 0., 4.2074983e-01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, 3.3060545e+00, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, 2.3994041e+00, 0., 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., 0., 4.2074983e-01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, 3.3060545e+00, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, 2.3994041e+00, 0., 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00},
        {5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., 0., -8.2292677e+01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0.},
        {0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, -7.9407372e+01, 7.7454940e-01, 0., 0., 0., 0., 0., 0.},
        {0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, -8.0314022e+01, 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., 0., -8.2292677e+01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0.},
        {0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, -7.9407372e+01, 7.7454940e-01, 0., 0., 0.},
        {0., 0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, -8.0314022e+01, 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., 0., -8.2292677e+01, 7.3598994e+00, 2.1304384e+00},
        {0., 0., 0., 0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, -7.9407372e+01, 7.7454940e-01},
        {0., 0., 0., 0., 0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, -8.0314022e+01}
    };

    variableVector resultCurrentPlasticMicroGradient;

    errorOut error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                              currentPlasticMacroVelocityGradient,
                                                                              currentPlasticMicroVelocityGradient,
                                                                              currentPlasticMicroGradientVelocityGradient,
                                                                              previousInversePlasticMicroDeformation,
                                                                              previousPlasticMicroGradient,
                                                                              previousPlasticMacroVelocityGradient,
                                                                              previousPlasticMicroVelocityGradient,
                                                                              previousPlasticMicroGradientVelocityGradient,
                                                                              resultCurrentPlasticMicroGradient, alpha );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerCurrentPlasticMicroGradient, resultCurrentPlasticMicroGradient ) );

    variableVector resultCurrentPlasticMicroGradient2;
    variableMatrix LHS2;

    error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                     currentPlasticMacroVelocityGradient,
                                                                     currentPlasticMicroVelocityGradient,
                                                                     currentPlasticMicroGradientVelocityGradient,
                                                                     previousInversePlasticMicroDeformation,
                                                                     previousPlasticMicroGradient,
                                                                     previousPlasticMacroVelocityGradient,
                                                                     previousPlasticMicroVelocityGradient,
                                                                     previousPlasticMicroGradientVelocityGradient,
                                                                     resultCurrentPlasticMicroGradient2, LHS2, alpha );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerCurrentPlasticMicroGradient, resultCurrentPlasticMicroGradient2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerLHS, LHS2 ) );

    //Test the Jacobians
    variableVector resultCurrentPlasticMicroGradientJ;

    variableMatrix dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                   dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                   dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                   dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient;

    error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, 
                                                                     currentPlasticMacroVelocityGradient,
                                                                     currentPlasticMicroVelocityGradient,
                                                                     currentPlasticMicroGradientVelocityGradient,
                                                                     previousInversePlasticMicroDeformation,
                                                                     previousPlasticMicroGradient,
                                                                     previousPlasticMacroVelocityGradient,
                                                                     previousPlasticMicroVelocityGradient,
                                                                     previousPlasticMicroGradientVelocityGradient,
                                                                     resultCurrentPlasticMicroGradientJ,
                                                                     dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                                                     dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                                                     dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                                                     dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                                                     alpha );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerCurrentPlasticMicroGradient, resultCurrentPlasticMicroGradientJ ) );

    //Test the jacobian w.r.t. the current plastic macro deformation
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < currentPlasticMicroDeformation.size(); i++ ){
        constantVector delta( currentPlasticMicroDeformation.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroDeformation[ i ] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation + delta,
                                                                         currentPlasticMacroVelocityGradient,
                                                                         currentPlasticMicroVelocityGradient,
                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                         previousInversePlasticMicroDeformation,
                                                                         previousPlasticMicroGradient,
                                                                         previousPlasticMacroVelocityGradient,
                                                                         previousPlasticMicroVelocityGradient,
                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                         resultP, alpha );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation - delta,
                                                                         currentPlasticMacroVelocityGradient,
                                                                         currentPlasticMicroVelocityGradient,
                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                         previousInversePlasticMicroDeformation,
                                                                         previousPlasticMicroGradient,
                                                                         previousPlasticMacroVelocityGradient,
                                                                         previousPlasticMicroVelocityGradient,
                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                         resultM, alpha );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMicroDeformation[ j ][ i ] ) );
        }
    }

    //Test the jacobian w.r.t. the current plastic macro velocity gradient
    for ( unsigned int i = 0; i < currentPlasticMacroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMacroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMacroVelocityGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                         currentPlasticMacroVelocityGradient + delta,
                                                                         currentPlasticMicroVelocityGradient,
                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                         previousInversePlasticMicroDeformation,
                                                                         previousPlasticMicroGradient,
                                                                         previousPlasticMacroVelocityGradient,
                                                                         previousPlasticMicroVelocityGradient,
                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                         resultP, alpha );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                         currentPlasticMacroVelocityGradient - delta,
                                                                         currentPlasticMicroVelocityGradient,
                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                         previousInversePlasticMicroDeformation,
                                                                         previousPlasticMicroGradient,
                                                                         previousPlasticMacroVelocityGradient,
                                                                         previousPlasticMicroVelocityGradient,
                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                         resultM, alpha );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient[ j ][ i ] ) );
        }
    }

    //Test the jacobian w.r.t. the current plastic micro velocity gradient
    for ( unsigned int i = 0; i < currentPlasticMicroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroVelocityGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                         currentPlasticMacroVelocityGradient,
                                                                         currentPlasticMicroVelocityGradient + delta,
                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                         previousInversePlasticMicroDeformation,
                                                                         previousPlasticMicroGradient,
                                                                         previousPlasticMacroVelocityGradient,
                                                                         previousPlasticMicroVelocityGradient,
                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                         resultP, alpha );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                         currentPlasticMacroVelocityGradient,
                                                                         currentPlasticMicroVelocityGradient - delta,
                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                         previousInversePlasticMicroDeformation,
                                                                         previousPlasticMicroGradient,
                                                                         previousPlasticMacroVelocityGradient,
                                                                         previousPlasticMicroVelocityGradient,
                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                         resultM, alpha );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient[ j ][ i ] ) );
        }
    }

    //Test the jacobian w.r.t. the current plastic micro gradient velocity gradient
    for ( unsigned int i = 0; i < currentPlasticMicroGradientVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroGradientVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroGradientVelocityGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                         currentPlasticMacroVelocityGradient,
                                                                         currentPlasticMicroVelocityGradient,
                                                                         currentPlasticMicroGradientVelocityGradient + delta,
                                                                         previousInversePlasticMicroDeformation,
                                                                         previousPlasticMicroGradient,
                                                                         previousPlasticMacroVelocityGradient,
                                                                         previousPlasticMicroVelocityGradient,
                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                         resultP, alpha );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                         currentPlasticMacroVelocityGradient,
                                                                         currentPlasticMicroVelocityGradient,
                                                                         currentPlasticMicroGradientVelocityGradient - delta,
                                                                         previousInversePlasticMicroDeformation,
                                                                         previousPlasticMicroGradient,
                                                                         previousPlasticMacroVelocityGradient,
                                                                         previousPlasticMicroVelocityGradient,
                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                         resultM, alpha );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testEvolvePlasticDeformation ){
    /*!
     * Evolve the plastic deformation.
     *
     */

    variableType Dt = 8.009359239014827;
    variableType alphaMacro = 0.48581069403548804;
    variableType alphaMicro = 0.827106668527372;
    variableType alphaMicroGrad = 0.608865458065455;

    variableVector currentPlasticMacroVelocityGradient = { -0.27014632, -0.35189361,  0.13880036,
                                                           -0.35584915, -0.48046671,  0.0912627 ,
                                                            0.47263309, -0.31133122, -0.2948745 };

    variableVector currentPlasticMicroVelocityGradient = { -0.01186714,  0.06090556,  0.40756696,
                                                            0.32909297,  0.28343008,  0.32993809,
                                                           -0.1899835 ,  0.41541764, -0.04761333 };

    variableVector currentPlasticMicroGradientVelocityGradient = { -0.44741758,  0.19632574,  0.20883732, -0.18740653,  0.30768804,
                                                                    0.46878663,  0.17691702, -0.07494477, -0.29871104, -0.46454855,
                                                                   -0.34339229, -0.24065171, -0.0737629 ,  0.48339896,  0.25954592,
                                                                    0.06579716, -0.43873995,  0.46482816,  0.20699556, -0.2417395 ,
                                                                   -0.46792573,  0.39616741, -0.22954849, -0.097399  ,  0.20759377,
                                                                   -0.35477122,  0.21945037 };

    variableVector previousPlasticDeformationGradient = { -0.45932572, -0.49691079, -0.13107169,
                                                           0.22431024, -0.61998429, -0.62903476,
                                                          -0.90544313,  0.41193656, -0.14749891 };

    variableVector previousPlasticMicroDeformation = { 0.28821738,  0.34407584,  0.06697721,
                                                       0.31626116, -0.62278299,  0.83699009,
                                                       0.16315786, -0.9434692 , -0.09991987 };

    variableVector previousPlasticMicroGradient = { -0.4908299 , -0.09107442,  0.11438472, -0.468681  , -0.09518528,
                                                     0.39306299,  0.30090607, -0.47472135, -0.39873768,  0.2071941 ,
                                                     0.3967767 , -0.40333359, -0.10263381, -0.39238343, -0.21260101,
                                                    -0.36133084, -0.32771357,  0.49025487,  0.01729508, -0.10338882,
                                                     0.15828675,  0.31044611, -0.30362281, -0.03972172, -0.377879  ,
                                                     0.11046526, -0.46997605 };

    variableVector previousPlasticMacroVelocityGradient = { -0.4362157 , -0.37728779,  0.40148895,
                                                            -0.15976616, -0.15839888,  0.26580961,
                                                             0.23215055, -0.02231055,  0.35815526 };

    variableVector previousPlasticMicroVelocityGradient = { 0.07222414, -0.40261712,  0.48556998,
                                                            0.45642476, -0.09409791, -0.49992287,
                                                            0.48252843,  0.48267731, -0.06878613 };

    variableVector previousPlasticMicroGradientVelocityGradient = { -0.15695492,  0.45802473, -0.38830808, -0.29241496, -0.34225716,
                                                                    -0.26127347, -0.41201504, -0.36558117, -0.24714561,  0.09340483,
                                                                     0.40887886, -0.1481247 ,  0.15128499, -0.04882876, -0.30445054,
                                                                     0.11557493, -0.20789811, -0.33913681,  0.3142761 ,  0.09871806,
                                                                     0.11316847, -0.45596559,  0.19030294, -0.33056333, -0.49391146,
                                                                     0.40129161, -0.06836289 };

    variableVector answerMacro = { -1.56924249,  1.75893631,  0.91035498,
                                    0.24241493, -0.43775652, -0.44317337,
                                   -2.6946146 ,  2.06645595,  0.86856673 };

    variableVector answerMicro = { 3.30247398,  2.51469186,  1.61673915,
                                   6.78559987, 11.60984671,  7.25516648,
                                   4.74326427,  4.32671877,  6.17698967 };

    variableVector answerMicroGrad = {  0.18184066,   0.06902067,  -4.65044761,  -4.6051252 ,
                                       -0.67925892,  -5.62411017,  -5.22449727,   0.2149805 ,
                                       -8.35852814, -14.60027145,  11.46879229,   1.1063327 ,
                                      -30.40158864,  24.16079446,  -3.60050806, -20.35119344,
                                       18.13574972,  -5.22922306,   4.00630219,  -2.39707255,
                                       -7.34817085,   3.12285942,  -6.24320576, -15.1020395 ,
                                       -5.11190824,   3.50558822,  -4.01621741 };

    variableVector resultMacro, resultMicro, resultMicroGrad;

    errorOut error = micromorphicElastoPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                             currentPlasticMicroVelocityGradient,
                                                                             currentPlasticMicroGradientVelocityGradient,
                                                                             previousPlasticDeformationGradient,
                                                                             previousPlasticMicroDeformation,
                                                                             previousPlasticMicroGradient,
                                                                             previousPlasticMacroVelocityGradient,
                                                                             previousPlasticMicroVelocityGradient,
                                                                             previousPlasticMicroGradientVelocityGradient,
                                                                             resultMacro, resultMicro, resultMicroGrad,
                                                                             alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );
    BOOST_CHECK( vectorTools::fuzzyEquals( resultMacro, answerMacro ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicro, answerMicro ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroGrad, answerMicroGrad ) );

    //Test the computation of the Jacobians
    variableVector resultMacroJ, resultMicroJ, resultMicroGradJ;
    variableMatrix dFdMacroL, dChidMicroL, dGradChidMacroL, dGradChidMicroL, dGradChidMicroGradL;

    error = micromorphicElastoPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                        currentPlasticMicroVelocityGradient,
                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                        previousPlasticDeformationGradient,
                                                                        previousPlasticMicroDeformation,
                                                                        previousPlasticMicroGradient,
                                                                        previousPlasticMacroVelocityGradient,
                                                                        previousPlasticMicroVelocityGradient,
                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                        resultMacroJ, resultMicroJ, resultMicroGradJ,
                                                                        dFdMacroL, dChidMicroL, dGradChidMacroL,
                                                                        dGradChidMicroL, dGradChidMicroGradL,
                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMacroJ, answerMacro ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroJ, answerMicro ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroGradJ, answerMicroGrad ) );

    //Test jacobians w.r.t. the current plastic macro velocity gradient
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < currentPlasticMacroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMacroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMacroVelocityGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        error = micromorphicElastoPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient + delta,
                                                                            currentPlasticMicroVelocityGradient,
                                                                            currentPlasticMicroGradientVelocityGradient,
                                                                            previousPlasticDeformationGradient,
                                                                            previousPlasticMicroDeformation,
                                                                            previousPlasticMicroGradient,
                                                                            previousPlasticMacroVelocityGradient,
                                                                            previousPlasticMicroVelocityGradient,
                                                                            previousPlasticMicroGradientVelocityGradient,
                                                                            resultMacroP, resultMicroP, resultMicroGradP,
                                                                            alphaMacro, alphaMicro, alphaMicroGrad );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient - delta,
                                                                            currentPlasticMicroVelocityGradient,
                                                                            currentPlasticMicroGradientVelocityGradient,
                                                                            previousPlasticDeformationGradient,
                                                                            previousPlasticMicroDeformation,
                                                                            previousPlasticMicroGradient,
                                                                            previousPlasticMacroVelocityGradient,
                                                                            previousPlasticMicroVelocityGradient,
                                                                            previousPlasticMicroGradientVelocityGradient,
                                                                            resultMacroM, resultMicroM, resultMicroGradM,
                                                                            alphaMacro, alphaMicro, alphaMicroGrad );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dFdMacroL[ j ][ i ] ) );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dGradChidMacroL[ j ][ i ] ) );
        }
    }

    for ( unsigned int i = 0; i < currentPlasticMicroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroVelocityGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        error = micromorphicElastoPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                            currentPlasticMicroVelocityGradient + delta,
                                                                            currentPlasticMicroGradientVelocityGradient,
                                                                            previousPlasticDeformationGradient,
                                                                            previousPlasticMicroDeformation,
                                                                            previousPlasticMicroGradient,
                                                                            previousPlasticMacroVelocityGradient,
                                                                            previousPlasticMicroVelocityGradient,
                                                                            previousPlasticMicroGradientVelocityGradient,
                                                                            resultMacroP, resultMicroP, resultMicroGradP,
                                                                            alphaMacro, alphaMicro, alphaMicroGrad );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                            currentPlasticMicroVelocityGradient - delta,
                                                                            currentPlasticMicroGradientVelocityGradient,
                                                                            previousPlasticDeformationGradient,
                                                                            previousPlasticMicroDeformation,
                                                                            previousPlasticMicroGradient,
                                                                            previousPlasticMacroVelocityGradient,
                                                                            previousPlasticMicroVelocityGradient,
                                                                            previousPlasticMicroGradientVelocityGradient,
                                                                            resultMacroM, resultMicroM, resultMicroGradM,
                                                                            alphaMacro, alphaMicro, alphaMicroGrad );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dChidMicroL[ j ][ i ] ) );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dGradChidMicroL[ j ][ i ] ) );
        }
    }

    for ( unsigned int i = 0; i < currentPlasticMicroGradientVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroGradientVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroGradientVelocityGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        error = micromorphicElastoPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                            currentPlasticMicroVelocityGradient,
                                                                            currentPlasticMicroGradientVelocityGradient + delta,
                                                                            previousPlasticDeformationGradient,
                                                                            previousPlasticMicroDeformation,
                                                                            previousPlasticMicroGradient,
                                                                            previousPlasticMacroVelocityGradient,
                                                                            previousPlasticMicroVelocityGradient,
                                                                            previousPlasticMicroGradientVelocityGradient,
                                                                            resultMacroP, resultMicroP, resultMicroGradP,
                                                                            alphaMacro, alphaMicro, alphaMicroGrad );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                            currentPlasticMicroVelocityGradient,
                                                                            currentPlasticMicroGradientVelocityGradient - delta,
                                                                            previousPlasticDeformationGradient,
                                                                            previousPlasticMicroDeformation,
                                                                            previousPlasticMicroGradient,
                                                                            previousPlasticMacroVelocityGradient,
                                                                            previousPlasticMicroVelocityGradient,
                                                                            previousPlasticMicroGradientVelocityGradient,
                                                                            resultMacroM, resultMicroM, resultMicroGradM,
                                                                            alphaMacro, alphaMicro, alphaMicroGrad );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dGradChidMicroGradL[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testEvolveStrainStateVariables){
    /*!
     * Test the evolution of the strain-like state variables.
     *
     */

    constantType Dt = 0.9173854839516489;

    parameterType alphaMacro = 0.5223948976356186;
    parameterType alphaMicro = 0.04303308197223432;
    parameterType alphaMicroGrad = 0.42138618879779943;

    variableType previousMacroISV = 0.7142072212214646;
    variableType previousMicroISV = 0.15443588642592365;
    variableVector previousMicroGradISV = { 0.05041408, 0.69624665, 0.8766614 };

    variableType currentMacroGamma = 0.17261697034297152;
    variableType currentMicroGamma = 0.02350474056209273;
    variableVector currentMicroGradientGamma = { 0.59406857, 0.79573586, 0.21213138 };

    variableType previousMacroGamma = 0.4782163679903608;
    variableType previousMicroGamma = 0.3546242470180624;
    variableVector previousMicroGradientGamma = { 0.6291708 , 0.48565385, 0.67132896 };

    variableType currentdGdMacroC = 0.2792660775467988;
    variableType currentdGdMicroC = 0.04579313341096025;
    variableMatrix currentdGdMicroGradC = { { 0.88556864, 0.08992741, 0.75316186 },
                                            { 0.76279627, 0.5635193 , 0.18529158 },
                                            { 0.05722408, 0.65275234, 0.97189144 } };

    variableType previousdGdMacroC = 0.6843499996792325;
    variableType previousdGdMicroC = 0.9662410574335287;
    variableMatrix previousdGdMicroGradC = { { 0.84374092, 0.21040392, 0.13887068 },
                                             { 0.34423717, 0.50801461, 0.28726825 },
                                             { 0.52590869, 0.36090934, 0.97602275 } };

    variableType answerMacroISV = 0.5362470356440043;
    variableType answerMicroISV = 0.13996373570538626;
    variableVector answerMicroGradISV = { -0.96380356,  0.11615298,  0.11045521 };

    variableType resultMacroISV, resultMicroISV;
    variableVector resultMicroGradISV;

    errorOut error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                               currentMicroGradientGamma, currentdGdMacroC,
                                                                               currentdGdMicroC, currentdGdMicroGradC,
                                                                               previousMacroISV, previousMicroISV,
                                                                               previousMicroGradISV, previousMacroGamma,
                                                                               previousMicroGamma, previousMicroGradientGamma,
                                                                               previousdGdMacroC, previousdGdMicroC,
                                                                               previousdGdMicroGradC, resultMacroISV,
                                                                               resultMicroISV, resultMicroGradISV,
                                                                               alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMacroISV, answerMacroISV ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroISV, answerMicroISV ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroGradISV, answerMicroGradISV ) );

    variableType resultMacroISVJ, resultMicroISVJ;
    variableVector resultMicroGradISVJ;

    variableType dMacroISVdMacroGamma, dMacroISVddGdMacroC;
    variableType dMicroISVdMicroGamma, dMicroISVddGdMicroC;
    variableMatrix dMicroGradientISVdMicroGradGamma, dMicroGradientISVddGdMicroGradGamma;

    //Test the Jacobians
    error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                      currentMicroGradientGamma, currentdGdMacroC,
                                                                      currentdGdMicroC, currentdGdMicroGradC,
                                                                      previousMacroISV, previousMicroISV,
                                                                      previousMicroGradISV, previousMacroGamma,
                                                                      previousMicroGamma, previousMicroGradientGamma,
                                                                      previousdGdMacroC, previousdGdMicroC,
                                                                      previousdGdMicroGradC, resultMacroISVJ,
                                                                      resultMicroISVJ, resultMicroGradISVJ,
                                                                      dMacroISVdMacroGamma, dMacroISVddGdMacroC,
                                                                      dMicroISVdMicroGamma, dMicroISVddGdMicroC,
                                                                      dMicroGradientISVdMicroGradGamma,
                                                                      dMicroGradientISVddGdMicroGradGamma,
                                                                      alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMacroISVJ, answerMacroISV ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroISVJ, answerMicroISV ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultMicroGradISVJ, answerMicroGradISV ) );

    //Test the Jacobians w.r.t. the current macro gamma
    constantType eps = 1e-6;

    constantType scalarDelta;
    scalarDelta = eps * fabs( currentMacroGamma ) + eps;

    variableType resultMacroISVP, resultMicroISVP;
    variableVector resultMicroGradISVP;
    variableType resultMacroISVM, resultMicroISVM;
    variableVector resultMicroGradISVM;

    error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma + scalarDelta, currentMicroGamma,
                                                                      currentMicroGradientGamma, currentdGdMacroC,
                                                                      currentdGdMicroC, currentdGdMicroGradC,
                                                                      previousMacroISV, previousMicroISV,
                                                                      previousMicroGradISV, previousMacroGamma,
                                                                      previousMicroGamma, previousMicroGradientGamma,
                                                                      previousdGdMacroC, previousdGdMicroC,
                                                                      previousdGdMicroGradC, resultMacroISVP,
                                                                      resultMicroISVP, resultMicroGradISVP,
                                                                      alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma - scalarDelta, currentMicroGamma,
                                                                      currentMicroGradientGamma, currentdGdMacroC,
                                                                      currentdGdMicroC, currentdGdMicroGradC,
                                                                      previousMacroISV, previousMicroISV,
                                                                      previousMicroGradISV, previousMacroGamma,
                                                                      previousMicroGamma, previousMicroGradientGamma,
                                                                      previousdGdMacroC, previousdGdMicroC,
                                                                      previousdGdMicroGradC, resultMacroISVM,
                                                                      resultMicroISVM, resultMicroGradISVM,
                                                                      alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    variableType gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradS, dMacroISVdMacroGamma ) );

    gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradS, 0. ) );

    variableVector gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * scalarDelta );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    //Test the Jacobians w.r.t. the current micro gamma
    scalarDelta = eps * fabs( currentMacroGamma ) + eps;

    error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma + scalarDelta,
                                                                      currentMicroGradientGamma, currentdGdMacroC,
                                                                      currentdGdMicroC, currentdGdMicroGradC,
                                                                      previousMacroISV, previousMicroISV,
                                                                      previousMicroGradISV, previousMacroGamma,
                                                                      previousMicroGamma, previousMicroGradientGamma,
                                                                      previousdGdMacroC, previousdGdMicroC,
                                                                      previousdGdMicroGradC, resultMacroISVP,
                                                                      resultMicroISVP, resultMicroGradISVP,
                                                                      alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma - scalarDelta,
                                                                      currentMicroGradientGamma, currentdGdMacroC,
                                                                      currentdGdMicroC, currentdGdMicroGradC,
                                                                      previousMacroISV, previousMicroISV,
                                                                      previousMicroGradISV, previousMacroGamma,
                                                                      previousMicroGamma, previousMicroGradientGamma,
                                                                      previousdGdMacroC, previousdGdMicroC,
                                                                      previousdGdMicroGradC, resultMacroISVM,
                                                                      resultMicroISVM, resultMicroGradISVM,
                                                                      alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradS, 0. ) );

    gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradS, dMicroISVdMicroGamma ) );

    gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * scalarDelta );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    //Test the Jacobians w.r.t. the current derivative of the macro plastic flow direction w.r.t. the macro cohesion
    scalarDelta = eps * fabs( currentdGdMacroC ) + eps;

    error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                      currentMicroGradientGamma, currentdGdMacroC + scalarDelta,
                                                                      currentdGdMicroC, currentdGdMicroGradC,
                                                                      previousMacroISV, previousMicroISV,
                                                                      previousMicroGradISV, previousMacroGamma,
                                                                      previousMicroGamma, previousMicroGradientGamma,
                                                                      previousdGdMacroC, previousdGdMicroC,
                                                                      previousdGdMicroGradC, resultMacroISVP,
                                                                      resultMicroISVP, resultMicroGradISVP,
                                                                      alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                      currentMicroGradientGamma, currentdGdMacroC - scalarDelta,
                                                                      currentdGdMicroC, currentdGdMicroGradC,
                                                                      previousMacroISV, previousMicroISV,
                                                                      previousMicroGradISV, previousMacroGamma,
                                                                      previousMicroGamma, previousMicroGradientGamma,
                                                                      previousdGdMacroC, previousdGdMicroC,
                                                                      previousdGdMicroGradC, resultMacroISVM,
                                                                      resultMicroISVM, resultMicroGradISVM,
                                                                      alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradS, dMacroISVddGdMacroC ) );

    gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradS, 0. ) );

    gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * scalarDelta );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    //Test the Jacobians w.r.t. the current derivative of the micro plastic flow direction w.r.t. the micro cohesion
    scalarDelta = eps * fabs( currentdGdMicroC ) + eps;

    error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                      currentMicroGradientGamma, currentdGdMacroC,
                                                                      currentdGdMicroC + scalarDelta, currentdGdMicroGradC,
                                                                      previousMacroISV, previousMicroISV,
                                                                      previousMicroGradISV, previousMacroGamma,
                                                                      previousMicroGamma, previousMicroGradientGamma,
                                                                      previousdGdMacroC, previousdGdMicroC,
                                                                      previousdGdMicroGradC, resultMacroISVP,
                                                                      resultMicroISVP, resultMicroGradISVP,
                                                                      alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                      currentMicroGradientGamma, currentdGdMacroC,
                                                                      currentdGdMicroC - scalarDelta, currentdGdMicroGradC,
                                                                      previousMacroISV, previousMicroISV,
                                                                      previousMicroGradISV, previousMacroGamma,
                                                                      previousMicroGamma, previousMicroGradientGamma,
                                                                      previousdGdMacroC, previousdGdMicroC,
                                                                      previousdGdMicroGradC, resultMacroISVM,
                                                                      resultMicroISVM, resultMicroGradISVM,
                                                                      alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_CHECK( !error );

    gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradS, 0. ) );

    gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * scalarDelta );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradS, dMicroISVddGdMicroC ) );

    gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * scalarDelta );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    //Test the Jacobians w.r.t. the micro gradient gamma
    for ( unsigned int i = 0; i < currentMicroGradientGamma.size(); i++ ){
        constantVector delta( currentMicroGradientGamma.size(), 0 );
        delta[ i ] = eps * fabs( currentMicroGradientGamma[ i ] ) + eps;

        error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                          currentMicroGradientGamma + delta, currentdGdMacroC,
                                                                          currentdGdMicroC, currentdGdMicroGradC,
                                                                          previousMacroISV, previousMicroISV,
                                                                          previousMicroGradISV, previousMacroGamma,
                                                                          previousMicroGamma, previousMicroGradientGamma,
                                                                          previousdGdMacroC, previousdGdMicroC,
                                                                          previousdGdMicroGradC, resultMacroISVP,
                                                                          resultMicroISVP, resultMicroGradISVP,
                                                                          alphaMacro, alphaMicro, alphaMicroGrad );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                          currentMicroGradientGamma - delta, currentdGdMacroC,
                                                                          currentdGdMicroC, currentdGdMicroGradC,
                                                                          previousMacroISV, previousMicroISV,
                                                                          previousMicroGradISV, previousMacroGamma,
                                                                          previousMicroGamma, previousMicroGradientGamma,
                                                                          previousdGdMacroC, previousdGdMicroC,
                                                                          previousdGdMicroGradC, resultMacroISVM,
                                                                          resultMicroISVM, resultMicroGradISVM,
                                                                          alphaMacro, alphaMicro, alphaMicroGrad );
    
        BOOST_CHECK( !error );
    
        gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * delta[ i ] );
    
        BOOST_CHECK( vectorTools::fuzzyEquals( gradS, 0. ) );
    
        gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * delta[ i ] );
    
        BOOST_CHECK( vectorTools::fuzzyEquals( gradS, 0. ) );
    
        gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * delta[ i ] );
    
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientISVdMicroGradGamma[ j ][ i ] ) );
        }
    }

    //Test the Jacobians w.r.t. the derivative of the micro gradient plastic potential w.r.t. the micro gradient cohesion
    for ( unsigned int i = 0; i < 9; i++ ){
        constantMatrix delta( 3, constantVector( 3, 0 ) );
        delta[ ( int )( i / 3) ][ i % 3 ] = eps * fabs( currentdGdMicroGradC[ ( int )( i / 3 ) ][ i % 3 ] ) + eps;

        error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                          currentMicroGradientGamma, currentdGdMacroC,
                                                                          currentdGdMicroC, currentdGdMicroGradC + delta,
                                                                          previousMacroISV, previousMicroISV,
                                                                          previousMicroGradISV, previousMacroGamma,
                                                                          previousMicroGamma, previousMicroGradientGamma,
                                                                          previousdGdMacroC, previousdGdMicroC,
                                                                          previousdGdMicroGradC, resultMacroISVP,
                                                                          resultMicroISVP, resultMicroGradISVP,
                                                                          alphaMacro, alphaMicro, alphaMicroGrad );
    
        BOOST_CHECK( !error );
    
        error = micromorphicElastoPlasticity::evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma,
                                                                          currentMicroGradientGamma, currentdGdMacroC,
                                                                          currentdGdMicroC, currentdGdMicroGradC - delta,
                                                                          previousMacroISV, previousMicroISV,
                                                                          previousMicroGradISV, previousMacroGamma,
                                                                          previousMicroGamma, previousMicroGradientGamma,
                                                                          previousdGdMacroC, previousdGdMicroC,
                                                                          previousdGdMicroGradC, resultMacroISVM,
                                                                          resultMicroISVM, resultMicroGradISVM,
                                                                          alphaMacro, alphaMicro, alphaMicroGrad );
    
        BOOST_CHECK( !error );
    
        gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );
    
        BOOST_CHECK( vectorTools::fuzzyEquals( gradS, 0. ) );
    
        gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );
    
        BOOST_CHECK( vectorTools::fuzzyEquals( gradS, 0. ) );
    
        gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );
    
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientISVddGdMicroGradGamma[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputeFlowDirections ){
    /*!
     * Test the computation of the flow directions
     *
     */

    variableVector PK2Stress = { 0.980355  ,  3.91541513,  1.67081083,
                                 4.54780598, -1.66521928,  4.96044281,
                                -0.12713388, -3.38729265, -2.25983197 };

    variableVector referenceMicroStress = { -2.01499646, -2.1610567 , -3.89718611,
                                            -2.1610567 , -0.43115989, -0.56605099,
                                            -3.89718611, -0.56605099,  4.22700636 };

    variableVector referenceHigherOrderStress = { -1.77580937, -1.71997457, -4.64296768,  3.20981958, -1.08382425,
                                                  -3.34229545,  4.85033094, -3.43682646,  4.09026212,  4.77299084,
                                                   3.93749094,  4.54620748, -2.79096599, -2.26200392,  4.61848971,
                                                  -4.82983291, -2.94604018,  0.91357337,  4.32181519, -4.30715354,
                                                   4.14245623,  4.62993794, -0.89968459, -0.31795132,  3.63048259,
                                                   2.65476864, -3.55307922 };

    variableType macroCohesion = 0.7183361521887071;
    variableType microCohesion = 0.6949323866523027;
    variableVector microGradientCohesion = { 0.713247  , 1.36791014, 1.5104695 };
    
    variableVector elasticRightCauchyGreen = { 4.59901102, 3.44800453, 0.77606303,
                                               3.44800453, 3.74960315, 0.52781964,
                                               0.77606303, 0.52781964, 1.93499034 };

    parameterVector macroFlowParameters = { 0.69418171, 0.48920844 };
    parameterVector microFlowParameters = { 0.92593374, 0.07788052 };
    parameterVector microGradientFlowParameters = { 0.75106525, 0.45320165 };

    variableVector answerMacroFlowDirection = { 3.38898459, 3.21365446, 0.74954484,
                                                3.24601877, 2.52116441, 0.6615649 ,
                                                0.65753016, 0.234347  , 1.18170509 };

    variableVector answerMicroFlowDirection = { 0.03163404, -0.5805495 , -0.3154144 ,
                                               -0.5805495 ,  0.22548302, -0.04324937,
                                               -0.3154144 , -0.04324937,  0.43080041 };

    variableVector answerMicroGradientFlowDirection = { 3.20270548,  0.        ,  0.        ,  3.07726673,  0.        ,
                                                        0.        ,  0.88460336,  0.        ,  0.        ,  3.15471725,
                                                        0.        ,  0.        ,  2.44225234,  0.        ,  0.        ,
                                                        0.16004057,  0.        ,  0.        ,  0.85841697,  0.        ,
                                                        0.        ,  0.62874426,  0.        ,  0.        ,  1.46940889,
                                                        0.        ,  0.        ,  0.        ,  0.76329042,  0.        ,
                                                        0.        ,  0.17122146,  0.        ,  0.        , -0.31051637,
                                                        0.        ,  0.        ,  0.73871243,  0.        ,  0.        ,
                                                        0.6214531 ,  0.        ,  0.        , -0.24342941,  0.        ,
                                                        0.        , -0.40887761,  0.        ,  0.        , -0.01215767,
                                                        0.        ,  0.        ,  0.79231477,  0.        ,  0.        ,
                                                        0.        ,  1.54713668,  0.        ,  0.        ,  1.17046377,
                                                        0.        ,  0.        ,  0.70180792,  0.        ,  0.        ,
                                                        1.88446996,  0.        ,  0.        ,  2.02250715,  0.        ,
                                                        0.        ,  0.30838528,  0.        ,  0.        ,  0.70653213,
                                                        0.        ,  0.        ,  0.19691721,  0.        ,  0.        ,
                                                        0.50658513 };

    variableType answerdGdMacroCohesion = -1.1365150471282142;
    variableType answerdGdMicroCohesion = -0.9616229151402682;
    variableMatrix answerdGdMicroGradientCohesion = { { -1.08210161, -0.        , -0.        },
                                                      { -0.        , -1.08210161, -0.        },
                                                      { -0.        , -0.        , -1.08210161 } };

    variableVector resultMacroFlowDirection, resultMicroFlowDirection, resultMicroGradientFlowDirection;

    variableType resultdGdMacroCohesion, resultdGdMicroCohesion;
    variableMatrix resultdGdMicroGradientCohesion;

    errorOut error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                          macroCohesion, microCohesion, microGradientCohesion,
                                                                          elasticRightCauchyGreen, macroFlowParameters,
                                                                          microFlowParameters, microGradientFlowParameters,
                                                                          resultMacroFlowDirection, resultMicroFlowDirection,
                                                                          resultMicroGradientFlowDirection, resultdGdMacroCohesion,
                                                                          resultdGdMicroCohesion, resultdGdMicroGradientCohesion );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMacroFlowDirection, resultMacroFlowDirection ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroFlowDirection, resultMicroFlowDirection ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroGradientFlowDirection, resultMicroGradientFlowDirection ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerdGdMacroCohesion, resultdGdMacroCohesion ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerdGdMicroCohesion, resultdGdMicroCohesion ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerdGdMicroGradientCohesion, resultdGdMicroGradientCohesion ) );

    //Tests of the Jacobian
    variableVector resultMacroFlowDirectionJ, resultMicroFlowDirectionJ, resultMicroGradientFlowDirectionJ;

    variableType resultdGdMacroCohesionJ, resultdGdMicroCohesionJ;
    variableMatrix resultdGdMicroGradientCohesionJ;

    variableMatrix dMacroFlowDirectiondPK2Stress, dMacroFlowDirectiondElasticRCG;
    variableMatrix dMicroFlowDirectiondSigma, dMicroFlowDirectiondElasticRCG;
    variableMatrix dMicroGradientFlowDirectiondM, dMicroGradientFlowDirectiondElasticRCG;

    error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                 macroCohesion, microCohesion, microGradientCohesion,
                                                                 elasticRightCauchyGreen, macroFlowParameters,
                                                                 microFlowParameters, microGradientFlowParameters,
                                                                 resultMacroFlowDirectionJ, resultMicroFlowDirectionJ,
                                                                 resultMicroGradientFlowDirectionJ, resultdGdMacroCohesionJ,
                                                                 resultdGdMicroCohesionJ, resultdGdMicroGradientCohesionJ,
                                                                 dMacroFlowDirectiondPK2Stress, dMacroFlowDirectiondElasticRCG,
                                                                 dMicroFlowDirectiondSigma, dMicroFlowDirectiondElasticRCG,
                                                                 dMicroGradientFlowDirectiondM,
                                                                 dMicroGradientFlowDirectiondElasticRCG );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMacroFlowDirection, resultMacroFlowDirectionJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroFlowDirection, resultMicroFlowDirectionJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerMicroGradientFlowDirection, resultMicroGradientFlowDirectionJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerdGdMacroCohesion, resultdGdMacroCohesionJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerdGdMicroCohesion, resultdGdMicroCohesionJ ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerdGdMicroGradientCohesion, resultdGdMicroGradientCohesionJ ) );

    //Test the Jacobians w.r.t. the PK2 stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < PK2Stress.size(); i++ ){
        constantVector delta( PK2Stress.size(), 0 );
        delta[ i ] = eps * fabs( PK2Stress[ i ] ) + eps;

        variableVector resultMacroFlowDirectionP, resultMicroFlowDirectionP, resultMicroGradientFlowDirectionP;
        variableVector resultMacroFlowDirectionM, resultMicroFlowDirectionM, resultMicroGradientFlowDirectionM;

        variableType resultdGdMacroCohesionP, resultdGdMicroCohesionP;
        variableType resultdGdMacroCohesionM, resultdGdMicroCohesionM;

        variableMatrix resultdGdMicroGradientCohesionP;
        variableMatrix resultdGdMicroGradientCohesionM;

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress + delta, referenceMicroStress, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionP, resultMicroFlowDirectionP,
                                                                     resultMicroGradientFlowDirectionP, resultdGdMacroCohesionP,
                                                                     resultdGdMicroCohesionP, resultdGdMicroGradientCohesionP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress - delta, referenceMicroStress, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMacroFlowDirectiondPK2Stress[ j ][ i ] ) );
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) );
            }
        }
    }

    //Test the jacobians w.r.t. the reference micro stress
    for ( unsigned int i = 0; i < referenceMicroStress.size(); i++ ){
        constantVector delta( referenceMicroStress.size(), 0 );
        delta[ i ] = eps * fabs( referenceMicroStress[ i ] ) + eps;

        variableVector resultMacroFlowDirectionP, resultMicroFlowDirectionP, resultMicroGradientFlowDirectionP;
        variableVector resultMacroFlowDirectionM, resultMicroFlowDirectionM, resultMicroGradientFlowDirectionM;

        variableType resultdGdMacroCohesionP, resultdGdMicroCohesionP;
        variableType resultdGdMacroCohesionM, resultdGdMicroCohesionM;

        variableMatrix resultdGdMicroGradientCohesionP;
        variableMatrix resultdGdMicroGradientCohesionM;

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress + delta, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionP, resultMicroFlowDirectionP,
                                                                     resultMicroGradientFlowDirectionP, resultdGdMacroCohesionP,
                                                                     resultdGdMicroCohesionP, resultdGdMicroGradientCohesionP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress - delta, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroFlowDirectiondSigma[ j ][ i ] ) );
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) );
            }
        }
    }

    //Test the jacobians w.r.t. the reference higher order stress
    for ( unsigned int i = 0; i < referenceHigherOrderStress.size(); i++ ){
        constantVector delta( referenceHigherOrderStress.size(), 0 );
        delta[ i ] = eps * fabs( referenceHigherOrderStress[ i ] ) + eps;

        variableVector resultMacroFlowDirectionP, resultMicroFlowDirectionP, resultMicroGradientFlowDirectionP;
        variableVector resultMacroFlowDirectionM, resultMicroFlowDirectionM, resultMicroGradientFlowDirectionM;

        variableType resultdGdMacroCohesionP, resultdGdMicroCohesionP;
        variableType resultdGdMacroCohesionM, resultdGdMicroCohesionM;

        variableMatrix resultdGdMicroGradientCohesionP;
        variableMatrix resultdGdMicroGradientCohesionM;

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress + delta,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionP, resultMicroFlowDirectionP,
                                                                     resultMicroGradientFlowDirectionP, resultdGdMacroCohesionP,
                                                                     resultdGdMicroCohesionP, resultdGdMicroGradientCohesionP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress - delta,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFlowDirectiondM[ j ][ i ] ) );
        }

        variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) );
            }
        }
    }

    //Test the jacobians w.r.t. the elastic right Cauchy-Green
    for ( unsigned int i = 0; i < elasticRightCauchyGreen.size(); i++ ){
        constantVector delta( elasticRightCauchyGreen.size(), 0 );
        delta[ i ] = eps * fabs( elasticRightCauchyGreen[ i ] ) + eps;

        variableVector resultMacroFlowDirectionP, resultMicroFlowDirectionP, resultMicroGradientFlowDirectionP;
        variableVector resultMacroFlowDirectionM, resultMicroFlowDirectionM, resultMicroGradientFlowDirectionM;

        variableType resultdGdMacroCohesionP, resultdGdMicroCohesionP;
        variableType resultdGdMacroCohesionM, resultdGdMicroCohesionM;

        variableMatrix resultdGdMicroGradientCohesionP;
        variableMatrix resultdGdMicroGradientCohesionM;

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen + delta, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionP, resultMicroFlowDirectionP,
                                                                     resultMicroGradientFlowDirectionP, resultdGdMacroCohesionP,
                                                                     resultdGdMicroCohesionP, resultdGdMicroGradientCohesionP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen - delta, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMacroFlowDirectiondElasticRCG[ j ][ i ] ) );
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroFlowDirectiondElasticRCG[ j ][ i ] ) );
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFlowDirectiondElasticRCG[ j ][ i ] ) );
        }

        variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) );
            }
        }
    }

    //Test the Jacobian w.r.t. the macro Cohesion
    constantType deltaScalar = eps * fabs( macroCohesion) + eps;

    variableVector resultMacroFlowDirectionP, resultMicroFlowDirectionP, resultMicroGradientFlowDirectionP;
    variableVector resultMacroFlowDirectionM, resultMicroFlowDirectionM, resultMicroGradientFlowDirectionM;

    variableType resultdGdMacroCohesionP, resultdGdMicroCohesionP;
    variableType resultdGdMacroCohesionM, resultdGdMicroCohesionM;

    variableMatrix resultdGdMicroGradientCohesionP;
    variableMatrix resultdGdMicroGradientCohesionM;

    error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                 macroCohesion + deltaScalar, microCohesion, microGradientCohesion,
                                                                 elasticRightCauchyGreen, macroFlowParameters,
                                                                 microFlowParameters, microGradientFlowParameters,
                                                                 resultMacroFlowDirectionP, resultMicroFlowDirectionP,
                                                                 resultMicroGradientFlowDirectionP, resultdGdMacroCohesionP,
                                                                 resultdGdMicroCohesionP, resultdGdMicroGradientCohesionP );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                 macroCohesion - deltaScalar, microCohesion, microGradientCohesion,
                                                                 elasticRightCauchyGreen, macroFlowParameters,
                                                                 microFlowParameters, microGradientFlowParameters,
                                                                 resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                 resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                 resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

    BOOST_CHECK( !error );

    variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * deltaScalar );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

    gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * deltaScalar );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

    variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradMat.size(); j++ ){
        for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) );
        }
    }

    //Test the Jacobian w.r.t. the micro Cohesion
    deltaScalar = eps * fabs( microCohesion) + eps;

    error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                 macroCohesion, microCohesion + deltaScalar, microGradientCohesion,
                                                                 elasticRightCauchyGreen, macroFlowParameters,
                                                                 microFlowParameters, microGradientFlowParameters,
                                                                 resultMacroFlowDirectionP, resultMicroFlowDirectionP,
                                                                 resultMicroGradientFlowDirectionP, resultdGdMacroCohesionP,
                                                                 resultdGdMicroCohesionP, resultdGdMicroGradientCohesionP );

    BOOST_CHECK( !error );

    error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                 macroCohesion, microCohesion - deltaScalar, microGradientCohesion,
                                                                 elasticRightCauchyGreen, macroFlowParameters,
                                                                 microFlowParameters, microGradientFlowParameters,
                                                                 resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                 resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                 resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

    BOOST_CHECK( !error );

    gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
    }

    gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * deltaScalar );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

    gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * deltaScalar );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

    gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradMat.size(); j++ ){
        for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) );
        }
    }

    //Test the jacobians w.r.t. the micro gradient cohesion
    for ( unsigned int i = 0; i < microGradientCohesion.size(); i++ ){
        constantVector delta( microGradientCohesion.size(), 0 );
        delta[ i ] = eps * fabs( microGradientCohesion[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion + delta,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionP, resultMicroFlowDirectionP,
                                                                     resultMicroGradientFlowDirectionP, resultdGdMacroCohesionP,
                                                                     resultdGdMicroCohesionP, resultdGdMicroGradientCohesionP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion - delta,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        BOOST_CHECK( !error );

        gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradScalar, 0. ) );

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) );
            }
        }
    }

}

BOOST_AUTO_TEST_CASE( testExtractMaterialParameters ){
    /*!
     * Test the extraction of the material parameters.
     *
     */

    std::vector< double > fparams = { 2, 0.53895133, 0.37172145,
                                      2, 0.37773052, 0.92739145,
                                      2, 0.53186824, 0.75454313,
                                      3, 0.95338442, 0.74042148, 0.09916127,
                                      3, 0.38093104, 0.49241325, 0.46187452,
                                      5, 0.82121039, 0.90566759, 0.50466975, 0.04830311, 0.85951495,
                                      3, 0.01166325, 0.05331896, 0.28081774,
                                      3, 0.32982199, 0.60161431, 0.33157768,
                                      5, 0.58881096, 0.11473813, 0.58001078, 0.83382529, 0.3260217,
                                      2, 1.7, 1.8,
                                      5, 2.8, .76, .15, 9.8, 5.4,
                                      11, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.,
                                      2, .76, 5.4,
                                      0.1, 0.2, 0.3, 1.1, 2.2 };

    parameterVector macroHardeningParameters;
    parameterVector microHardeningParameters;
    parameterVector microGradientHardeningParameters;
    parameterVector macroFlowParameters;
    parameterVector microFlowParameters;
    parameterVector microGradientFlowParameters;
    parameterVector macroYieldParameters;
    parameterVector microYieldParameters;
    parameterVector microGradientYieldParameters;
    parameterVector Amatrix;
    parameterVector Bmatrix;
    parameterVector Cmatrix;
    parameterVector Dmatrix;
    constantType alphaMacro;
    constantType alphaMicro;
    constantType alphaMicroGradient;
    constantType relativeTolerance;
    constantType absoluteTolerance;

    parameterVector answerMacroHardeningParameters = { 0.53895133, 0.37172145 };
    parameterVector answerMicroHardeningParameters = { 0.37773052, 0.92739145 };
    parameterVector answerMicroGradientHardeningParameters = { 0.53186824, 0.75454313 };
    parameterVector answerMacroFlowParameters = { 0.95338442, 0.74042148, 0.09916127 };
    parameterVector answerMicroFlowParameters = { 0.38093104, 0.49241325, 0.46187452 };
    parameterVector answerMicroGradientFlowParameters = { 0.82121039, 0.90566759, 0.50466975, 0.04830311, 0.85951495 };
    parameterVector answerMacroYieldParameters = { 0.01166325, 0.05331896, 0.28081774 };
    parameterVector answerMicroYieldParameters = { 0.32982199, 0.60161431, 0.33157768 };
    parameterVector answerMicroGradientYieldParameters = { 0.58881096, 0.11473813, 0.58001078, 0.83382529, 0.3260217 };


    parameterVector answerAmatrix;
    parameterVector answerBmatrix;
    parameterVector answerCmatrix;
    parameterVector answerDmatrix;

    errorOut error = micromorphicLinearElasticity::formIsotropicA( 1.7, 1.8, answerAmatrix );
    BOOST_CHECK( !error );

    error = micromorphicLinearElasticity::formIsotropicB( 2.8, 0.76, 0.15, 9.8, 5.4, answerBmatrix );
    BOOST_CHECK( !error );

    error = micromorphicLinearElasticity::formIsotropicC( { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.}, answerCmatrix );
    BOOST_CHECK( !error );

    error = micromorphicLinearElasticity::formIsotropicD( 0.76, 5.4, answerDmatrix );
    BOOST_CHECK( !error );

    constantType answerAlphaMacro = 0.1;
    constantType answerAlphaMicro = 0.2;
    constantType answerAlphaMicroGradient = 0.3;
    constantType answerRelativeTolerance = 1.1;
    constantType answerAbsoluteTolerance = 2.2;

    error = micromorphicElastoPlasticity::extractMaterialParameters( fparams,
                macroHardeningParameters, microHardeningParameters, microGradientHardeningParameters,
                macroFlowParameters, microFlowParameters, microGradientFlowParameters,
                macroYieldParameters, microYieldParameters, microGradientYieldParameters,
                Amatrix, Bmatrix, Cmatrix, Dmatrix, alphaMacro, alphaMicro, alphaMicroGradient,
                relativeTolerance, absoluteTolerance );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( macroHardeningParameters, answerMacroHardeningParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microHardeningParameters, answerMicroHardeningParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microGradientHardeningParameters, answerMicroGradientHardeningParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( macroFlowParameters, answerMacroFlowParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microFlowParameters, answerMicroFlowParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microGradientFlowParameters, answerMicroGradientFlowParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( macroYieldParameters, answerMacroYieldParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microYieldParameters, answerMicroYieldParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microGradientYieldParameters, answerMicroGradientYieldParameters ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( Amatrix, answerAmatrix ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( Bmatrix, answerBmatrix ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( Cmatrix, answerCmatrix ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( Dmatrix, answerDmatrix ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( alphaMacro, answerAlphaMacro ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( alphaMicro, answerAlphaMicro ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( alphaMicroGradient, answerAlphaMicroGradient ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( alphaMicro, answerAlphaMicro ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( relativeTolerance, answerRelativeTolerance ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( absoluteTolerance, answerAbsoluteTolerance ) );

}

BOOST_AUTO_TEST_CASE( testExtractStateVariables ){
    /*!
     * Test the extraction of the state variables.
     *
     */

    std::vector< double > SDVS( 55 );
    for ( unsigned int i = 0; i < SDVS.size(); i++ ){
        SDVS[ i ] = i + 1;
    }

    variableType   answerPreviousMacroStrainISV                    = 1;
    variableType   answerPreviousMicroStrainISV                    = 2;
    variableVector answerPreviousMicroGradientStrainISV            = { 3, 4, 5 };
    variableType   answerPreviousMacroGamma                        = 6;
    variableType   answerPreviousMicroGamma                        = 7;
    variableVector answerPreviousMicroGradientGamma                = { 8, 9, 10 };
    variableVector answerPreviousPlasticDeformationGradient        = { 12, 12, 13, 14, 16, 16, 17, 18, 20 };
    variableVector answerPreviousPlasticMicroDeformation           = { 21, 21, 22, 23, 25, 25, 26, 27, 29 };
    variableVector answerPreviousPlasticGradientMicroDeformation   = { 29, 30, 31, 32, 33, 34, 35, 36, 37,
                                                                       38, 39, 40, 41, 42, 43, 44, 45, 46,
                                                                       47, 48, 49, 50, 51, 52, 53, 54, 55 };

    variableType   previousMacroStrainISV;
    variableType   previousMicroStrainISV;
    variableVector previousMicroGradientStrainISV;
    variableType   previousMacroGamma;
    variableType   previousMicroGamma;
    variableVector previousMicroGradientGamma;
    variableVector previousPlasticDeformationGradient;
    variableVector previousPlasticMicroDeformation;
    variableVector previousPlasticGradientMicroDeformation;

    errorOut error = micromorphicElastoPlasticity::extractStateVariables( SDVS,
                         previousMacroStrainISV, previousMicroStrainISV,
                         previousMicroGradientStrainISV,
                         previousMacroGamma, previousMicroGamma,
                         previousMicroGradientGamma,
                         previousPlasticDeformationGradient,
                         previousPlasticMicroDeformation,
                         previousPlasticGradientMicroDeformation );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerPreviousMacroStrainISV, previousMacroStrainISV ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerPreviousMicroStrainISV, previousMicroStrainISV ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerPreviousMicroGradientStrainISV, previousMicroGradientStrainISV ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerPreviousMacroGamma, previousMacroGamma ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerPreviousMicroGamma, previousMicroGamma ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerPreviousMicroGradientGamma, previousMicroGradientGamma ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerPreviousPlasticDeformationGradient, previousPlasticDeformationGradient ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerPreviousPlasticMicroDeformation, previousPlasticMicroDeformation ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( answerPreviousPlasticGradientMicroDeformation, previousPlasticGradientMicroDeformation ) );

}

BOOST_AUTO_TEST_CASE( testCout_redirect ){
    /*!
     * Test the utility function which redirects cout to a string buffer.
     *
     */

    std::stringbuf buffer;
    micromorphicElastoPlasticity::cout_redirect rd( &buffer );
    
    std::string answer = "hello world\n";

    std::cout << answer;

    BOOST_CHECK( answer.compare( buffer.str() ) == 0 );
}

BOOST_AUTO_TEST_CASE( testCerr_redirect ){
    /*!
     * Test the utility function which redirects cerr to a string buffer.
     *
     */

    std::stringbuf buffer;
    micromorphicElastoPlasticity::cerr_redirect rd( &buffer );
    
    std::string answer = "hello world\n";

    std::cerr << answer;

    BOOST_CHECK( answer.compare( buffer.str() ) == 0 );
}

BOOST_AUTO_TEST_CASE( testAssembleFundamentalDeformationMeasures ){
    /*!
     * Assemble the fundamental deformation measures from the degrees of freedom.
     *
     */

    double grad_u[ 3 ][ 3 ] = { { 1, 2, 3 },
                                { 4, 5, 6 },
                                { 7, 8, 9 } };

    double phi[ 9 ] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    double grad_phi[ 9 ][ 3 ] = { {  1,  2,  3 },
                                  {  4,  5,  6 },
                                  {  7,  8,  9 },
                                  { 10, 11, 12 },
                                  { 13, 14, 15 },
                                  { 16, 17, 18 },
                                  { 19, 20, 21 },
                                  { 22, 23, 24 },
                                  { 25, 26, 27 } };

    variableVector answerDeformationGradient = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector answerMicroDeformation = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector answerGradientMicroDeformation = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                     10, 11, 12, 13, 14, 15, 16, 17, 18,
                                                     19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector resultF, resultChi, resultGradChi;

    errorOut error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, grad_phi,
                                                                                           resultF, resultChi, resultGradChi );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultF, answerDeformationGradient ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultChi, answerMicroDeformation ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultGradChi, answerGradientMicroDeformation ) );

    //Test the Jacobians
    variableVector resultFJ, resultChiJ, resultGradChiJ;
    variableMatrix dFdGradU, dChidPhi, dGradChidGradPhi;

    error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, grad_phi,
                                                                                  resultFJ, resultChiJ, resultGradChiJ,
                                                                                  dFdGradU, dChidPhi, dGradChidGradPhi );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultFJ, answerDeformationGradient ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultChiJ, answerMicroDeformation ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultGradChiJ, answerGradientMicroDeformation ) );

    //Test the jacobians w.r.t. the gradient of the displacement
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){
        constantMatrix delta( 3, constantVector( 3, 0 ) );
        unsigned int ii, ij;
        ii = ( int )( i / 3 );
        ij = i % 3;
        delta[ ii ][ ij ] = eps * fabs( grad_u[ ii ][ ij ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 3 ][ 3 ] =
        { 
                { grad_u[ 0 ][ 0 ] + delta[ 0 ][ 0 ], grad_u[ 0 ][ 1 ] + delta[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] + delta[ 0 ][ 2 ] },
                { grad_u[ 1 ][ 0 ] + delta[ 1 ][ 0 ], grad_u[ 1 ][ 1 ] + delta[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] + delta[ 1 ][ 2 ] },
                { grad_u[ 2 ][ 0 ] + delta[ 2 ][ 0 ], grad_u[ 2 ][ 1 ] + delta[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] + delta[ 2 ][ 2 ] }
        };

        double negative_perturb[ 3 ][ 3 ] =
        { 
                { grad_u[ 0 ][ 0 ] - delta[ 0 ][ 0 ], grad_u[ 0 ][ 1 ] - delta[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] - delta[ 0 ][ 2 ] },
                { grad_u[ 1 ][ 0 ] - delta[ 1 ][ 0 ], grad_u[ 1 ][ 1 ] - delta[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] - delta[ 1 ][ 2 ] },
                { grad_u[ 2 ][ 0 ] - delta[ 2 ][ 0 ], grad_u[ 2 ][ 1 ] - delta[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] - delta[ 2 ][ 2 ] }
        };

        error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( positive_perturb, phi, grad_phi,
                                                                                      FP, chiP, gradChiP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( negative_perturb, phi, grad_phi,
                                                                                      FM, chiM, gradChiM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dFdGradU[ j ][ i ] ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    for ( unsigned int i = 0; i < 9; i++ ){
        constantVector delta( 9, 0 );

        delta[ i ] = eps * fabs( phi[ i ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 9 ] = { phi[ 0 ] + delta[ 0 ], phi[ 1 ] + delta[ 1 ], phi[ 2 ] + delta[ 2 ],
                                         phi[ 3 ] + delta[ 3 ], phi[ 4 ] + delta[ 4 ], phi[ 5 ] + delta[ 5 ],
                                         phi[ 6 ] + delta[ 6 ], phi[ 7 ] + delta[ 7 ], phi[ 8 ] + delta[ 8 ] };

        double negative_perturb[ 9 ] = { phi[ 0 ] - delta[ 0 ], phi[ 1 ] - delta[ 1 ], phi[ 2 ] - delta[ 2 ],
                                         phi[ 3 ] - delta[ 3 ], phi[ 4 ] - delta[ 4 ], phi[ 5 ] - delta[ 5 ],
                                         phi[ 6 ] - delta[ 6 ], phi[ 7 ] - delta[ 7 ], phi[ 8 ] - delta[ 8 ] };

        error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, positive_perturb, grad_phi,
                                                                                      FP, chiP, gradChiP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, negative_perturb, grad_phi,
                                                                                      FM, chiM, gradChiM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dChidPhi[ j ][ i ] ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    for ( unsigned int i = 0; i < 27; i++ ){
        constantMatrix delta( 9, constantVector( 3, 0 ) );
        unsigned int ii, ij;
        ii = ( int )( i / 3 );
        ij = i % 3;
        delta[ ii ][ ij ] = eps * fabs( grad_phi[ ii ][ ij ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 9 ][ 3 ] =
        { 
                { grad_phi[ 0 ][ 0 ] + delta[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ] + delta[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] + delta[ 0 ][ 2 ] },
                { grad_phi[ 1 ][ 0 ] + delta[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ] + delta[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] + delta[ 1 ][ 2 ] },
                { grad_phi[ 2 ][ 0 ] + delta[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ] + delta[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] + delta[ 2 ][ 2 ] },
                { grad_phi[ 3 ][ 0 ] + delta[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ] + delta[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] + delta[ 3 ][ 2 ] },
                { grad_phi[ 4 ][ 0 ] + delta[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ] + delta[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] + delta[ 4 ][ 2 ] },
                { grad_phi[ 5 ][ 0 ] + delta[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ] + delta[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] + delta[ 5 ][ 2 ] },
                { grad_phi[ 6 ][ 0 ] + delta[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ] + delta[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] + delta[ 6 ][ 2 ] },
                { grad_phi[ 7 ][ 0 ] + delta[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ] + delta[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] + delta[ 7 ][ 2 ] },
                { grad_phi[ 8 ][ 0 ] + delta[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ] + delta[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] + delta[ 8 ][ 2 ] }
        };

        double negative_perturb[ 9 ][ 3 ] =
        { 
                { grad_phi[ 0 ][ 0 ] - delta[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ] - delta[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] - delta[ 0 ][ 2 ] },
                { grad_phi[ 1 ][ 0 ] - delta[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ] - delta[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] - delta[ 1 ][ 2 ] },
                { grad_phi[ 2 ][ 0 ] - delta[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ] - delta[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] - delta[ 2 ][ 2 ] },
                { grad_phi[ 3 ][ 0 ] - delta[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ] - delta[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] - delta[ 3 ][ 2 ] },
                { grad_phi[ 4 ][ 0 ] - delta[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ] - delta[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] - delta[ 4 ][ 2 ] },
                { grad_phi[ 5 ][ 0 ] - delta[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ] - delta[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] - delta[ 5 ][ 2 ] },
                { grad_phi[ 6 ][ 0 ] - delta[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ] - delta[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] - delta[ 6 ][ 2 ] },
                { grad_phi[ 7 ][ 0 ] - delta[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ] - delta[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] - delta[ 7 ][ 2 ] },
                { grad_phi[ 8 ][ 0 ] - delta[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ] - delta[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] - delta[ 8 ][ 2 ] }
        };


        error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, positive_perturb,
                                                                                      FP, chiP, gradChiP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, negative_perturb,
                                                                                      FM, chiM, gradChiM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dGradChidGradPhi[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testEvaluateYieldFunctions ){
    /*!
     * Test the evaluation of the yield functions
     *
     */

    variableVector PK2Stress = { 0.65888653,  0.28604558, -0.44426912,
                                -0.89889651,  0.56484547, -0.49370633,
                                -0.72307938, -0.43117725, -0.90089686 };

    variableVector SigmaStress = { -0.7631354 , -0.94148632,  0.31898031,
                                    0.89228848,  0.84753794,  0.99573296,
                                    0.0069115 ,  0.28727492,  0.82575735 };

    variableVector M = { -0.56449022, -0.15407578, -0.06806782, -0.48688421, -0.76199393,
                         -0.0703529 ,  0.36092558, -0.30212673,  0.81992342, -0.13140498,
                          0.25451202, -0.11514735, -0.75631506,  0.49592347, -0.31198427,
                          0.43412034, -0.03285418,  0.49030702,  0.76019333,  0.95293589,
                          0.23643603, -0.1391037 , -0.6446606 ,  0.17408786, -0.41085493,
                         -0.60240843,  0.70821091 };

    variableVector elasticC = { 1.13538358, 0.11776884, 0.13734646,
                                0.11776884, 1.00826486, 0.0787335 ,
                                0.13734646, 0.0787335 , 1.02987928 };

    variableType   macroCohesion = 0.8860253435105722;
    variableType   microCohesion = 0.18795616072213417;
    variableVector microGradientCohesion = { 0.23149834, 0.43664147, 0.30363429 };

    variableType macroPhi = 0.6598820660198553;
    variableType microPhi = 0.140325385292413;
    variableType microGradientPhi = 0.3194920699960573;

    variableType macroBeta = 0.1366298928020886;
    variableType microBeta = 0.7392551172586164;
    variableType microGradientBeta = 0.5515550713283653;

    variableVector macroYieldParameters = { macroPhi, macroBeta };
    variableVector microYieldParameters = { microPhi, microBeta };
    variableVector microGradientYieldParameters = { microGradientPhi, microGradientBeta };

    constantType answerMacroYield = 0.8067218309832007;
    constantType answerMicroYield = 1.9183045092031357;
    constantVector answerMicroGradientYield = { 0.45351001, 0.94557435, 0.91947682 }; 

    variableVector answer = { answerMacroYield, answerMicroYield, answerMicroGradientYield[ 0 ],
                              answerMicroGradientYield[ 1 ], answerMicroGradientYield[ 2 ] };

    variableVector result;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG;
    errorOut error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                           microGradientCohesion, elasticC, macroYieldParameters,
                                                                           microYieldParameters, microGradientYieldParameters,
                                                                           result, DEBUG );
#else
    errorOut error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                           microGradientCohesion, elasticC, macroYieldParameters,
                                                                           microYieldParameters, microGradientYieldParameters,
                                                                           result );
#endif

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    //Test the jacobians
    variableVector resultJ;

    variableVector dMacroFdPK2;
    variableType   dMacroFdMacroC;
    variableVector dMacroFdElasticRCG;
    variableVector dMicroFdSigma;
    variableType   dMicroFdMicroC;
    variableVector dMicroFdElasticRCG;
    variableMatrix dMicroGradientFdM;
    variableMatrix dMicroGradientFdMicroGradientC;
    variableMatrix dMicroGradientFdElasticRCG;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG_JACOBIAN;
    error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                  microGradientCohesion, elasticC, macroYieldParameters,
                                                                  microYieldParameters, microGradientYieldParameters,
                                                                  resultJ,
                                                                  dMacroFdPK2, dMacroFdMacroC, dMacroFdElasticRCG,
                                                                  dMicroFdSigma, dMicroFdMicroC, dMicroFdElasticRCG,
                                                                  dMicroGradientFdM, dMicroGradientFdMicroGradientC,
                                                                  dMicroGradientFdElasticRCG,
                                                                  DEBUG_JACOBIAN );
#else
    error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                  microGradientCohesion, elasticC, macroYieldParameters,
                                                                  microYieldParameters, microGradientYieldParameters,
                                                                  resultJ,
                                                                  dMacroFdPK2, dMacroFdMacroC, dMacroFdElasticRCG,
                                                                  dMicroFdSigma, dMicroFdMicroC, dMicroFdElasticRCG,
                                                                  dMicroGradientFdM, dMicroGradientFdMicroGradientC,
                                                                  dMicroGradientFdElasticRCG );
#endif

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( resultJ, answer ) );

    //Test Jacobians w.r.t. the PK2 stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < PK2Stress.size(); i++ ){
        constantVector delta( PK2Stress.size(), 0 );
        delta[ i ] = eps * fabs( PK2Stress[ i ] ) + eps;

        variableVector resultP;
        variableVector resultM;

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector > DEBUGP, DEBUGM;
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress + delta, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP, DEBUGP );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress + delta, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP );
#endif

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress - delta, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM, DEBUGM );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress - delta, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM );
#endif

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 0 ], dMacroFdPK2[ i ] ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 1 ], 0. ) );

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    //Test Jacobians w.r.t. the reference micro stress
    for ( unsigned int i = 0; i < SigmaStress.size(); i++ ){
        constantVector delta( SigmaStress.size(), 0 );
        delta[ i ] = eps * fabs( SigmaStress[ i ] ) + eps;

        variableVector resultP;
        variableVector resultM;

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector > DEBUGP, DEBUGM;
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress + delta, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP, DEBUGP );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress + delta, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP );
#endif

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress - delta, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM, DEBUGM );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress - delta, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM );
#endif

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 0 ], 0. ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 1 ], dMicroFdSigma[ i ] ) );

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    //Test Jacobians w.r.t. the higher order stress
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[ i ] = eps * fabs( M[ i ] ) + eps;

        variableVector resultP;
        variableVector resultM;

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector > DEBUGP, DEBUGM;
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M + delta, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP, DEBUGP );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M + delta, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP );
#endif

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M - delta, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM, DEBUGM );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M - delta, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM );
#endif

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 0 ], 0. ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 1 ], 0. ) );

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFdM[ j - 2 ][ i ] ) );
        }
    }

    //Test jacobians w.r.t. the macro cohesion
    for ( unsigned int i = 0; i < 1; i++ ){
        constantType delta;
        delta = eps * fabs( macroCohesion ) + eps;

        variableVector resultP;
        variableVector resultM;

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector > DEBUGP, DEBUGM;
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion + delta, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP, DEBUGP );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion + delta, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP );
#endif

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion - delta, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM, DEBUGM );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion - delta, microCohesion,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM );
#endif

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 0 ], dMacroFdMacroC ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 1 ], 0. ) );

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    //Test jacobians w.r.t. the micro cohesion
    for ( unsigned int i = 0; i < 1; i++ ){
        constantType delta;
        delta = eps * fabs( microCohesion ) + eps;

        variableVector resultP;
        variableVector resultM;

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector > DEBUGP, DEBUGM;
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion + delta,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP, DEBUGP );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion + delta,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP );
#endif

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion - delta,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM, DEBUGM );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion - delta,
                                                                      microGradientCohesion, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM );
#endif

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 0 ], 0. ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 1 ], dMicroFdMicroC ) );

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    //Test jacobians w.r.t. the micro gradient cohesion
    for ( unsigned int i = 0; i < microGradientCohesion.size(); i++ ){
        constantVector delta( microGradientCohesion.size(), 0 );
        delta[ i ] = eps * fabs( microGradientCohesion[ i ] ) + eps;

        variableVector resultP;
        variableVector resultM;

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector > DEBUGP, DEBUGM;
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion + delta, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP, DEBUGP );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion  + delta, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP );
#endif

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion - delta, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM, DEBUGM );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion - delta, elasticC, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM );
#endif

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 0 ], 0. ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 1 ], 0. ) );

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFdMicroGradientC[ j - 2 ][ i ] ) );
        }
    }

    for ( unsigned int i = 0; i < elasticC.size(); i++ ){
        constantVector delta( elasticC.size(), 0 );
        delta[ i ] = eps * fabs( elasticC[ i ] ) + eps;

        variableVector resultP;
        variableVector resultM;

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector > DEBUGP, DEBUGM;
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC + delta, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP, DEBUGP );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC + delta, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultP );
#endif

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC - delta, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM, DEBUGM );
#else
        error = micromorphicElastoPlasticity::evaluateYieldFunctions( PK2Stress, SigmaStress, M, macroCohesion, microCohesion,
                                                                      microGradientCohesion, elasticC - delta, macroYieldParameters,
                                                                      microYieldParameters, microGradientYieldParameters,
                                                                      resultM );
#endif

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 0 ], dMacroFdElasticRCG[ i ] ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ 1 ], dMicroFdElasticRCG[ i ] ) );

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFdElasticRCG[ j - 2 ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputeCohesion ){
    /*!
     * Test the computation of the cohesion
     *
     */

    variableType   macroStrainISV = 1.0;
    variableType   microStrainISV = 2.0;
    variableVector microGradientStrainISV = { 3.0, 4.0, 5.0 };

    parameterVector macroHardeningParameters = { 0.1, 0.2 };
    parameterVector microHardeningParameters = { 0.3, 0.4 };
    parameterVector microGradientHardeningParameters = { 0.5, 0.6 };

    variableType   answerMacroC = 0.1 + 0.2 * 1.0;
    variableType   answerMicroC = 0.3 + 0.4 * 2.0;
    variableVector answerMicroGradientC = { 0.5 + 0.6 * 3.0, 0.5 + 0.6 * 4.0, 0.5 + 0.6 * 5.0 };

    variableType   macroC, microC;
    variableVector microGradientC;

    errorOut error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV, microStrainISV, microGradientStrainISV,
                                                                    macroHardeningParameters, microHardeningParameters,
                                                                    microGradientHardeningParameters, macroC, microC,
                                                                    microGradientC );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( macroC, answerMacroC ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microC, answerMicroC ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microGradientC, answerMicroGradientC ) );

    //Test the Jacobians
    variableType   macroCJ, microCJ;
    variableVector microGradientCJ;

    variableType dMacroCdMacroStrainISV, dMicroCdMicroStrainISV;
    variableMatrix dMicroGradientCdMicroStrainISV;

    error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV, microStrainISV, microGradientStrainISV,
                                                           macroHardeningParameters, microHardeningParameters,
                                                           microGradientHardeningParameters, macroCJ, microCJ,
                                                           microGradientCJ, dMacroCdMacroStrainISV,
                                                           dMicroCdMicroStrainISV, dMicroGradientCdMicroStrainISV );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( macroCJ, answerMacroC ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microCJ, answerMicroC ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( microGradientCJ, answerMicroGradientC ) );

    //Test Jacobian w.r.t. the macro strain ISV
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 1; i++ ){
        constantType delta = eps * fabs( macroStrainISV ) + eps;

        variableType macroCP, microCP;
        variableType macroCM, microCM;

        variableVector microGradientCP;
        variableVector microGradientCM;

        error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV + delta, microStrainISV, microGradientStrainISV,
                                                               macroHardeningParameters, microHardeningParameters,
                                                               microGradientHardeningParameters, macroCP, microCP,
                                                               microGradientCP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV - delta, microStrainISV, microGradientStrainISV,
                                                               macroHardeningParameters, microHardeningParameters,
                                                               microGradientHardeningParameters, macroCM, microCM,
                                                               microGradientCM );

        BOOST_CHECK( !error );

        variableType grad = ( macroCP - macroCM ) / ( 2 * delta );

        BOOST_CHECK( vectorTools::fuzzyEquals( grad, dMacroCdMacroStrainISV ) );

        grad = ( microCP - microCM ) / ( 2 * delta );

        BOOST_CHECK( vectorTools::fuzzyEquals( grad, 0. ) );

        variableVector gradCol = ( microGradientCP - microGradientCM ) / ( 2 * delta );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    //Test Jacobian w.r.t. the micro strain ISV
    for ( unsigned int i = 0; i < 1; i++ ){
        constantType delta = eps * fabs( microStrainISV ) + eps;

        variableType macroCP, microCP;
        variableType macroCM, microCM;

        variableVector microGradientCP;
        variableVector microGradientCM;

        error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV, microStrainISV + delta, microGradientStrainISV,
                                                               macroHardeningParameters, microHardeningParameters,
                                                               microGradientHardeningParameters, macroCP, microCP,
                                                               microGradientCP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV, microStrainISV - delta, microGradientStrainISV,
                                                               macroHardeningParameters, microHardeningParameters,
                                                               microGradientHardeningParameters, macroCM, microCM,
                                                               microGradientCM );

        BOOST_CHECK( !error );

        variableType grad = ( macroCP - macroCM ) / ( 2 * delta );

        BOOST_CHECK( vectorTools::fuzzyEquals( grad, 0. ) );

        grad = ( microCP - microCM ) / ( 2 * delta );

        BOOST_CHECK( vectorTools::fuzzyEquals( grad, dMicroCdMicroStrainISV ) );

        variableVector gradCol = ( microGradientCP - microGradientCM ) / ( 2 * delta );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    //Test Jacobian w.r.t. the micro gradient strain ISV
    for ( unsigned int i = 0; i < microGradientStrainISV.size(); i++ ){
        constantVector delta( microGradientStrainISV.size(), 0 );
        delta[ i ] = eps * fabs( microGradientStrainISV[ i ] ) + eps;

        variableType macroCP, microCP;
        variableType macroCM, microCM;

        variableVector microGradientCP;
        variableVector microGradientCM;

        error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV, microStrainISV, microGradientStrainISV + delta,
                                                               macroHardeningParameters, microHardeningParameters,
                                                               microGradientHardeningParameters, macroCP, microCP,
                                                               microGradientCP );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV, microStrainISV, microGradientStrainISV - delta,
                                                               macroHardeningParameters, microHardeningParameters,
                                                               microGradientHardeningParameters, macroCM, microCM,
                                                               microGradientCM );

        BOOST_CHECK( !error );

        variableType grad = ( macroCP - macroCM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( grad, 0. ) );

        grad = ( microCP - microCM ) / ( 2 * delta[ i ] );

        BOOST_CHECK( vectorTools::fuzzyEquals( grad, 0. ) );

        variableVector gradCol = ( microGradientCP - microGradientCM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientCdMicroStrainISV[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testEvaluate_model){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time
    std::vector< double > time = { 10., 2.5 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 2.4e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.4e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e0, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double current_grad_u[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
                                        {0.100, 0.001, 0.000 },
                                        {0.000, 0.000, 0.000 } };

    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
                                         {0, 0, 0},
                                         {0, 0, 0} };
    //Initialize the micro displacement
//    double current_phi[ 9 ] = { 0.84729289,  0.40617104,  0.59534561,  
//                                0.44195587,  0.34121966, -0.79098944, 
//                               -0.43965428,  0.88466225,  0.1684519 };
//
//    double previous_phi[ 9 ] = { -0.99935855, -0.21425717,  0.0668254 ,
//                                 -0.11111872, -0.07416114, -1.01048108,
//                                  0.1804018 , -1.01116291,  0.03248007 };

    double current_phi[ 9 ] = { 0.100, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000 };

    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
                                          { -0.18572739,  0.06847269,  0.22931628 },
                                          { -0.01829735, -0.48731265, -0.25277529 },
                                          {  0.26626212,  0.4844646 , -0.31965177 },
                                          {  0.49197846,  0.19051656, -0.0365349  },
                                          { -0.06607774, -0.33526875, -0.15803078 },
                                          {  0.09738707, -0.49482218, -0.39584868 },
                                          { -0.45599864,  0.08585038, -0.09432794 },
                                          {  0.23055539,  0.07564162,  0.24051469 } };

//    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
//                                           {  0.00294417,  0.34480654, -0.34450988 },
//                                           {  0.21056511, -0.28113967, -0.45726839 },
//                                           { -0.26431286, -0.09985721,  0.47322301 },
//                                           { -0.18156887, -0.32226199, -0.37295847 },
//                                           {  0.15062371,  0.09439471,  0.09167948 },
//                                           { -0.46869859,  0.018301  ,  0.45013866 },
//                                           { -0.15455446,  0.40552715, -0.4216042  },
//                                           { -0.38930237,  0.10974753, -0.31188239 } };

//    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0} };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };
                                           

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > current_PK2( 9, 0 );

    std::vector< double > current_SIGMA( 9, 0 );

    std::vector< double > current_M( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

#ifdef DEBUG_MODE
    solverTools::homotopyMap homotopyDEBUG;
#endif

    solverTools::floatVector PK2Answer = { 172.484,   15.3785,   -0.917177,
                                            13.4848, 142.823,    -0.0214307,
                                            -1.7635,   1.77719, 141.069 };

    solverTools::floatVector SigmaAnswer = { 176.916,   15.8646,   -2.83731,
                                              15.8646, 144.538,     1.85836,
                                              -2.83731,  1.85836, 142.013 };

    solverTools::floatVector MAnswer = { 0.598283, -0.512218,  0.620664,    3.22636,   1.16682,
                                         1.20593,   0.562825, -2.52317,     1.62616,  -2.61391,
                                        -0.60994,  -1.02147,   0.668187,    0.49348,  -0.23916,
                                        -2.77419,   0.760483,  1.71784,    -0.499389,  2.62828,
                                        -0.761044,  1.23369,  -0.00778206, -2.25643,  -0.729551,
                                         0.743204,  0.910521 };

    solverTools::floatVector SDVSAnswer = { -1.79592e-24,  0.0243222,    0.0822384,    0.0430345,   0.0435752,
                                            -8.96006e-25,  0.00852191,   0.0465339,    0.0243507,   0.0246566,
                                             0.00742998,   0.00500421,  -0.000296486,  0.00498757, -0.00260492,
                                             0.000284355, -0.000367318,  0.000222511, -0.0015603,   0.00863313,
                                             0.00537105,  -0.000347686,  0.00643802,  -0.00298667,  0.000297105,
                                            -0.000422398,  0.000293946, -0.00181986,   0.0385522,  -0.0244021,
                                            -0.00912035,   0.0105171,   -0.000328615,  0.0222843,   0.0185626,
                                            -0.0118234,   -0.00785555,   0.0451085,    0.031607,    0.0212748,
                                             0.0116981,    0.0161821,    0.00126031,   0.0147688,   0.00691805,
                                            -0.0241431,    0.00942608,  -0.0350366,    0.0221571,  -0.0249697,
                                             0.00935849,   0.0214931,    0.0169609,    0.0177352,   0.0203838 };

    std::vector< double > SDVS = SDVSDefault;

    int errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                                  current_grad_u,  current_phi,  current_grad_phi,
                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                  SDVS,
                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                  current_PK2, current_SIGMA, current_M,
                                                                  ADD_TERMS,
                                                                  output_message
#ifdef DEBUG_MODE
                                                                  , homotopyDEBUG
#endif
                                                                  );

    BOOST_CHECK( errorCode == 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVSAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2Answer, current_PK2, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SigmaAnswer, current_SIGMA, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( MAnswer, current_M, 1e-5 ) );

    //Test the Jacobians
    std::vector< std::vector< double > > DPK2Dgrad_u, DPK2Dphi, DPK2Dgrad_phi, DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                         DMDgrad_u, DMDphi, DMDgrad_phi;

    std::vector< std::vector< std::vector< double > > > ADD_JACOBIANS;

    SDVS = SDVSDefault;

#ifdef DEBUG_MODE
    homotopyDEBUG.clear();
#endif

    errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                              current_grad_u,  current_phi,  current_grad_phi,
                                                              previous_grad_u, previous_phi, previous_grad_phi,
                                                              SDVS,
                                                              current_ADD_DOF,  current_ADD_grad_DOF,
                                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                                              current_PK2, current_SIGMA, current_M,
                                                              DPK2Dgrad_u, DPK2Dphi, DPK2Dgrad_phi,
                                                              DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                                              DMDgrad_u, DMDphi, DMDgrad_phi,
                                                              ADD_TERMS, ADD_JACOBIANS,
                                                              output_message
#ifdef DEBUG_MODE
                                                              , homotopyDEBUG
#endif
                                                              );
#ifdef DEBUG_MODE
    solverTools::debugMap DEBUG = homotopyDEBUG[ "converged_values" ][ "converged_values" ];
    BOOST_CHECK( DEBUG.size() != 0 );
#endif

    BOOST_CHECK( errorCode == 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVSAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2Answer, current_PK2, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SigmaAnswer, current_SIGMA, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( MAnswer, current_M, 1e-5 ) );

    //Test the jacobians w.r.t. the gradient of the macro displacement
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){
        constantMatrix delta( 3, constantVector( 3, 0 ) );
        delta[ i / 3 ][ i % 3 ] = eps * fabs( current_grad_u[ i / 3 ][ i % 3 ] ) + eps;

        constantMatrix gradU =
        { 
            { current_grad_u[ 0 ][ 0 ], current_grad_u[ 0 ][ 1 ], current_grad_u[ 0 ][ 2 ] },
            { current_grad_u[ 1 ][ 0 ], current_grad_u[ 1 ][ 1 ], current_grad_u[ 1 ][ 2 ] },
            { current_grad_u[ 2 ][ 0 ], current_grad_u[ 2 ][ 1 ], current_grad_u[ 2 ][ 2 ] }
        };

        constantMatrix gradU_P = gradU + delta;
        constantMatrix gradU_M = gradU - delta;

        double current_grad_u_P[ 3 ][ 3 ] =
        {
            { gradU_P[ 0 ][ 0 ], gradU_P[ 0 ][ 1 ], gradU_P[ 0 ][ 2 ] },
            { gradU_P[ 1 ][ 0 ], gradU_P[ 1 ][ 1 ], gradU_P[ 1 ][ 2 ] },
            { gradU_P[ 2 ][ 0 ], gradU_P[ 2 ][ 1 ], gradU_P[ 2 ][ 2 ] }
        };

        double current_grad_u_M[ 3 ][ 3 ] =
        {
            { gradU_M[ 0 ][ 0 ], gradU_M[ 0 ][ 1 ], gradU_M[ 0 ][ 2 ] },
            { gradU_M[ 1 ][ 0 ], gradU_M[ 1 ][ 1 ], gradU_M[ 1 ][ 2 ] },
            { gradU_M[ 2 ][ 0 ], gradU_M[ 2 ][ 1 ], gradU_M[ 2 ][ 2 ] }
        };

        solverTools::floatVector PK2_P, SIGMA_P, M_P;
        solverTools::floatVector PK2_M, SIGMA_M, M_M;

        solverTools::floatVector SDVS_P = SDVSDefault;
        solverTools::floatVector SDVS_M = SDVSDefault;

#ifdef DEBUG_MODE
        solverTools::homotopyMap DEBUG_P, DEBUG_M;
#endif

        errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                                  current_grad_u_P,  current_phi,  current_grad_phi,
                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                  SDVS_P,
                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                  PK2_P, SIGMA_P, M_P,
                                                                  ADD_TERMS,
                                                                  output_message
#ifdef DEBUG_MODE
                                                                  , DEBUG_P
#endif
                                                                  );

        BOOST_CHECK( errorCode == 0 );

        errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                                  current_grad_u_M,  current_phi,  current_grad_phi,
                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                  SDVS_M,
                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                  PK2_M, SIGMA_M, M_M,
                                                                  ADD_TERMS,
                                                                  output_message
#ifdef DEBUG_MODE
                                                                  , DEBUG_M
#endif
                                                                  );

        BOOST_CHECK( errorCode == 0 );

#ifdef DEBUG_MODE
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P[ "converged_values" ][ "converged_values" ].begin();
                   it_P != DEBUG_P[ "converged_values" ][ "converged_values" ].end(); it_P++ ){
            
            auto it_M = DEBUG_M[ "converged_values" ][ "converged_values" ].find( it_P->first );
            BOOST_CHECK( it_M != DEBUG_M[ "converged_values" ][ "converged_values" ].end() );

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i / 3 ][ i % 3 ] ) );

        }

        BOOST_CHECK( numericGradients.size() != 0 );

        //Check the total Jacobians of the plastic deformation measures
        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticDeformationGradient" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticDeformationGradient" ][ j ],
                                            DEBUG[ "totaldPlasticDeformationGradientdDeformationGradient" ][ 9 * j + i ],
                                            1e-6, 1e-6 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticMicroDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticMicroDeformationdDeformationGradient" ][ 9 * j + i ],
                                            1e-6, 1e-6 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticGradientMicroDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticGradientMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticGradientMicroDeformationdDeformationGradient" ][ 9 * j + i ],
                                            1e-6, 1e-6 ) );
        }

        //Check the total Jacobians of the intermediate stresses
        for ( unsigned int j = 0; j < numericGradients[ "intermediatePK2Stress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "intermediatePK2Stress" ][ j ],
                                            DEBUG[ "totaldPK2StressdDeformationGradient" ][ 9 * j + i ],
                                            1e-5, 1e-5 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceMicroStress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceMicroStress" ][ j ],
                                            DEBUG[ "totaldReferenceMicroStressdDeformationGradient" ][ 9 * j + i ],
                                            1e-5, 1e-5 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceHigherOrderStress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceHigherOrderStress" ][ j ],
                                            DEBUG[ "totaldReferenceHigherOrderStressdDeformationGradient" ][ 9 * j + i ],
                                            1e-5, 1e-5 ) );
        }
#endif

        solverTools::floatVector gradCol = ( PK2_P - PK2_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], DPK2Dgrad_u[ j ][ i ], 1e-4, 1e-5 ) );
        }

        gradCol = ( SIGMA_P - SIGMA_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], DSIGMADgrad_u[ j ][ i ], 1e-4, 1e-5 ) );
        }

        gradCol = ( M_P - M_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], DMDgrad_u[ j ][ i ], 1e-4, 1e-5 ) );
        }
    }

    //Test the jacobians w.r.t. the micro displacement
    for ( unsigned int i = 0; i < 9; i++ ){
        constantVector delta( 9, 0 );
        delta[ i ] = eps * fabs( current_phi[ i ] ) + eps;

        constantVector phi =
        {
            current_phi[ 0 ], current_phi[ 1 ], current_phi[ 2 ],
            current_phi[ 3 ], current_phi[ 4 ], current_phi[ 5 ],
            current_phi[ 6 ], current_phi[ 7 ], current_phi[ 8 ]
        };

        constantVector phi_P = phi + delta;
        constantVector phi_M = phi - delta;

        double current_phi_P[ 9 ] =
        {
            phi_P[ 0 ], phi_P[ 1 ], phi_P[ 2 ],
            phi_P[ 3 ], phi_P[ 4 ], phi_P[ 5 ],
            phi_P[ 6 ], phi_P[ 7 ], phi_P[ 8 ]
        };

        double current_phi_M[ 9 ] =
        {
            phi_M[ 0 ], phi_M[ 1 ], phi_M[ 2 ],
            phi_M[ 3 ], phi_M[ 4 ], phi_M[ 5 ],
            phi_M[ 6 ], phi_M[ 7 ], phi_M[ 8 ]
        };

        solverTools::floatVector PK2_P, SIGMA_P, M_P;
        solverTools::floatVector PK2_M, SIGMA_M, M_M;

        solverTools::floatVector SDVS_P = SDVSDefault;
        solverTools::floatVector SDVS_M = SDVSDefault;

#ifdef DEBUG_MODE
        solverTools::homotopyMap DEBUG_P, DEBUG_M;
#endif

        errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                                  current_grad_u,  current_phi_P,  current_grad_phi,
                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                  SDVS_P,
                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                  PK2_P, SIGMA_P, M_P,
                                                                  ADD_TERMS,
                                                                  output_message
#ifdef DEBUG_MODE
                                                                  , DEBUG_P
#endif
                                                                  );

        BOOST_CHECK( errorCode == 0 );

        errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                                  current_grad_u,  current_phi_M,  current_grad_phi,
                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                  SDVS_M,
                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                  PK2_M, SIGMA_M, M_M,
                                                                  ADD_TERMS,
                                                                  output_message
#ifdef DEBUG_MODE
                                                                  , DEBUG_M
#endif
                                                                  );

        BOOST_CHECK( errorCode == 0 );

#ifdef DEBUG_MODE
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P[ "converged_values" ][ "converged_values" ].begin();
                   it_P != DEBUG_P[ "converged_values" ][ "converged_values" ].end(); it_P++ ){
            
            auto it_M = DEBUG_M[ "converged_values" ][ "converged_values" ].find( it_P->first );
            BOOST_CHECK( it_M != DEBUG_M[ "converged_values" ][ "converged_values" ].end() );

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );

        }

        BOOST_CHECK( numericGradients.size() != 0 );

        //Check the total Jacobians of the plastic deformation measures
        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticDeformationGradient" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticDeformationGradient" ][ j ],
                                            DEBUG[ "totaldPlasticDeformationGradientdMicroDeformation" ][ 9 * j + i ],
                                            1e-6, 1e-6 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticMicroDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticMicroDeformationdMicroDeformation" ][ 9 * j + i ],
                                            1e-6, 1e-6 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticGradientMicroDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticGradientMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticGradientMicroDeformationdMicroDeformation" ][ 9 * j + i ],
                                            1e-6, 1e-6 ) );
        }

        //Check the total Jacobians of the intermediate stresses
        for ( unsigned int j = 0; j < numericGradients[ "intermediatePK2Stress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "intermediatePK2Stress" ][ j ],
                                            DEBUG[ "totaldPK2StressdMicroDeformation" ][ 9 * j + i ],
                                            1e-5, 1e-5 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceMicroStress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceMicroStress" ][ j ],
                                            DEBUG[ "totaldReferenceMicroStressdMicroDeformation" ][ 9 * j + i ],
                                            1e-5, 1e-5 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceHigherOrderStress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceHigherOrderStress" ][ j ],
                                            DEBUG[ "totaldReferenceHigherOrderStressdMicroDeformation" ][ 9 * j + i ],
                                            1e-5, 1e-5 ) );
        }

#endif

        solverTools::floatVector gradCol = ( PK2_P - PK2_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], DPK2Dphi[ j ][ i ], 1e-4 ) );
        }

        gradCol = ( SIGMA_P - SIGMA_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], DSIGMADphi[ j ][ i ], 1e-4 ) );
        }

        gradCol = ( M_P - M_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], DMDphi[ j ][ i ], 1e-4 ) );
        }
    }

    //Test the jacobians w.r.t. the gradient of the micro displacement
    for ( unsigned int i = 0; i < 27; i++ ){
        constantMatrix delta( 9, constantVector( 3, 0 ) );
        delta[ i / 3 ][ i % 3 ] = eps * fabs( current_grad_phi[ i / 3 ][ i % 3 ] ) + eps;

        constantMatrix gradPhi =
        { 
            { current_grad_phi[ 0 ][ 0 ], current_grad_phi[ 0 ][ 1 ], current_grad_phi[ 0 ][ 2 ] },
            { current_grad_phi[ 1 ][ 0 ], current_grad_phi[ 1 ][ 1 ], current_grad_phi[ 1 ][ 2 ] },
            { current_grad_phi[ 2 ][ 0 ], current_grad_phi[ 2 ][ 1 ], current_grad_phi[ 2 ][ 2 ] },
            { current_grad_phi[ 3 ][ 0 ], current_grad_phi[ 3 ][ 1 ], current_grad_phi[ 3 ][ 2 ] },
            { current_grad_phi[ 4 ][ 0 ], current_grad_phi[ 4 ][ 1 ], current_grad_phi[ 4 ][ 2 ] },
            { current_grad_phi[ 5 ][ 0 ], current_grad_phi[ 5 ][ 1 ], current_grad_phi[ 5 ][ 2 ] },
            { current_grad_phi[ 6 ][ 0 ], current_grad_phi[ 6 ][ 1 ], current_grad_phi[ 6 ][ 2 ] },
            { current_grad_phi[ 7 ][ 0 ], current_grad_phi[ 7 ][ 1 ], current_grad_phi[ 7 ][ 2 ] },
            { current_grad_phi[ 8 ][ 0 ], current_grad_phi[ 8 ][ 1 ], current_grad_phi[ 8 ][ 2 ] }
        };

        constantMatrix gradPhi_P = gradPhi + delta;
        constantMatrix gradPhi_M = gradPhi - delta;

        double current_grad_phi_P[ 9 ][ 3 ] =
        {
            { gradPhi_P[ 0 ][ 0 ], gradPhi_P[ 0 ][ 1 ], gradPhi_P[ 0 ][ 2 ] },
            { gradPhi_P[ 1 ][ 0 ], gradPhi_P[ 1 ][ 1 ], gradPhi_P[ 1 ][ 2 ] },
            { gradPhi_P[ 2 ][ 0 ], gradPhi_P[ 2 ][ 1 ], gradPhi_P[ 2 ][ 2 ] },
            { gradPhi_P[ 3 ][ 0 ], gradPhi_P[ 3 ][ 1 ], gradPhi_P[ 3 ][ 2 ] },
            { gradPhi_P[ 4 ][ 0 ], gradPhi_P[ 4 ][ 1 ], gradPhi_P[ 4 ][ 2 ] },
            { gradPhi_P[ 5 ][ 0 ], gradPhi_P[ 5 ][ 1 ], gradPhi_P[ 5 ][ 2 ] },
            { gradPhi_P[ 6 ][ 0 ], gradPhi_P[ 6 ][ 1 ], gradPhi_P[ 6 ][ 2 ] },
            { gradPhi_P[ 7 ][ 0 ], gradPhi_P[ 7 ][ 1 ], gradPhi_P[ 7 ][ 2 ] },
            { gradPhi_P[ 8 ][ 0 ], gradPhi_P[ 8 ][ 1 ], gradPhi_P[ 8 ][ 2 ] }
        };

        double current_grad_phi_M[ 9 ][ 3 ] =
        {
            { gradPhi_M[ 0 ][ 0 ], gradPhi_M[ 0 ][ 1 ], gradPhi_M[ 0 ][ 2 ] },
            { gradPhi_M[ 1 ][ 0 ], gradPhi_M[ 1 ][ 1 ], gradPhi_M[ 1 ][ 2 ] },
            { gradPhi_M[ 2 ][ 0 ], gradPhi_M[ 2 ][ 1 ], gradPhi_M[ 2 ][ 2 ] },
            { gradPhi_M[ 3 ][ 0 ], gradPhi_M[ 3 ][ 1 ], gradPhi_M[ 3 ][ 2 ] },
            { gradPhi_M[ 4 ][ 0 ], gradPhi_M[ 4 ][ 1 ], gradPhi_M[ 4 ][ 2 ] },
            { gradPhi_M[ 5 ][ 0 ], gradPhi_M[ 5 ][ 1 ], gradPhi_M[ 5 ][ 2 ] },
            { gradPhi_M[ 6 ][ 0 ], gradPhi_M[ 6 ][ 1 ], gradPhi_M[ 6 ][ 2 ] },
            { gradPhi_M[ 7 ][ 0 ], gradPhi_M[ 7 ][ 1 ], gradPhi_M[ 7 ][ 2 ] },
            { gradPhi_M[ 8 ][ 0 ], gradPhi_M[ 8 ][ 1 ], gradPhi_M[ 8 ][ 2 ] }
        };

        solverTools::floatVector PK2_P, SIGMA_P, M_P;
        solverTools::floatVector PK2_M, SIGMA_M, M_M;

        solverTools::floatVector SDVS_P = SDVSDefault;
        solverTools::floatVector SDVS_M = SDVSDefault;

#ifdef DEBUG_MODE
        solverTools::homotopyMap DEBUG_P, DEBUG_M;
#endif

        errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                                  current_grad_u,  current_phi,  current_grad_phi_P,
                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                  SDVS_P,
                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                  PK2_P, SIGMA_P, M_P,
                                                                  ADD_TERMS,
                                                                  output_message
#ifdef DEBUG_MODE
                                                                  , DEBUG_P
#endif
                                                                  );

        BOOST_CHECK( errorCode == 0 );

        errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                                  current_grad_u,  current_phi,  current_grad_phi_M,
                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                  SDVS_M,
                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                  PK2_M, SIGMA_M, M_M,
                                                                  ADD_TERMS,
                                                                  output_message
#ifdef DEBUG_MODE
                                                                  , DEBUG_M
#endif
                                                                  );

        BOOST_CHECK( errorCode == 0 );

#ifdef DEBUG_MODE
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P[ "converged_values" ][ "converged_values" ].begin();
                   it_P != DEBUG_P[ "converged_values" ][ "converged_values" ].end(); it_P++ ){
            
            auto it_M = DEBUG_M[ "converged_values" ][ "converged_values" ].find( it_P->first );
            BOOST_CHECK( it_M != DEBUG_M[ "converged_values" ][ "converged_values" ].end() );

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i / 3 ][ i % 3 ] ) );

        }

        BOOST_CHECK( numericGradients.size() != 0 );

        //Check the total Jacobians of the plastic deformation measures
        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticDeformationGradient" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticDeformationGradient" ][ j ],
                                            DEBUG[ "totaldPlasticDeformationGradientdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-6, 1e-8 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticMicroDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticMicroDeformationdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-6, 1e-8 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticGradientMicroDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticGradientMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticGradientMicroDeformationdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-6, 1e-8 ) );
        }

        //Check the total Jacobians of the intermediate stresses
        for ( unsigned int j = 0; j < numericGradients[ "intermediatePK2Stress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "intermediatePK2Stress" ][ j ],
                                            DEBUG[ "totaldPK2StressdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-5, 1e-5 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceMicroStress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceMicroStress" ][ j ],
                                            DEBUG[ "totaldReferenceMicroStressdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-5, 1e-5 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceHigherOrderStress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceHigherOrderStress" ][ j ],
                                            DEBUG[ "totaldReferenceHigherOrderStressdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-5, 1e-5 ) );
        }
#endif

        solverTools::floatVector gradCol = ( PK2_P - PK2_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], DPK2Dgrad_phi[ j ][ i ] ) );
        }

//        std::cout << "numeric DPK2Dgrad_u:\n"; vectorTools::print( gradCol );

        gradCol = ( SIGMA_P - SIGMA_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], DSIGMADgrad_phi[ j ][ i ], 1e-5 ) );
        }

//        std::cout << "numeric DSIGMADphi:\n"; vectorTools::print( gradCol );

        gradCol = ( M_P - M_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], DMDgrad_phi[ j ][ i ], 1e-5 ) );
        }

//        std::cout << "numeric DMDgrad_phi:\n"; vectorTools::print( gradCol );

    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticDeformationResidual ){
    /*!
     * Test the computation of the plastic deformation residual.
     *
     */

    std::vector< double > fparams = { 2, 1.0e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.5e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e2, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    parameterVector macroHardeningParameters;
    parameterVector microHardeningParameters;
    parameterVector microGradientHardeningParameters;

    parameterVector macroFlowParameters;
    parameterVector microFlowParameters;
    parameterVector microGradientFlowParameters;

    parameterVector macroYieldParameters;
    parameterVector microYieldParameters;
    parameterVector microGradientYieldParameters;

    parameterVector Amatrix, Bmatrix, Cmatrix, Dmatrix;
    parameterType alphaMacro, alphaMicro, alphaMicroGradient;
    parameterType relativeTolerance, absoluteTolerance;

    errorOut error = micromorphicElastoPlasticity::extractMaterialParameters( fparams,
                                                                              macroHardeningParameters, microHardeningParameters,
                                                                              microGradientHardeningParameters,
                                                                              macroFlowParameters, microFlowParameters,
                                                                              microGradientFlowParameters,
                                                                              macroYieldParameters, microYieldParameters,
                                                                              microGradientYieldParameters,
                                                                              Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                              alphaMacro, alphaMicro, alphaMicroGradient,
                                                                              relativeTolerance, absoluteTolerance );

    BOOST_CHECK( !error );

    constantType   Dt = 2.5;
    variableType   currentMacroGamma = 0.01;
    variableType   currentMicroGamma = 0.02;
    variableVector currentMicroGradientGamma = {0.011, 0.021, 0.031};
    variableType   previousMacroGamma = 0.;
    variableType   previousMicroGamma = 0.;
    variableVector previousMicroGradientGamma( 3, 0 );
    variableType   previousMacroStrainISV = 0.;
    variableType   previousMicroStrainISV = 0.;
    variableVector previousMicroGradientStrainISV( 3., 0. );
    variableVector currentDeformationGradient = { 0.04482969,  0.88562312, -0.38710144,
                                                 -0.93722716,  0.19666568,  0.41677155,
                                                  0.46929057,  0.33779672,  0.81392228 };
    variableVector currentMicroDeformation = {  0.51930689,  0.27954023, -0.85955731,
                                                0.09469279, -0.99381243, -0.23218079,
                                               -0.82281393,  0.09643296, -0.54637704 };
    variableVector currentGradientMicroDeformation = { 0.04176306, -0.0151958 , -0.00090558, -0.01844751,  0.04512391,
                                                       0.02174263, -0.00508239, -0.01827377,  0.00541031, -0.01330239,
                                                      -0.02479987,  0.02914825,  0.00168841,  0.00230506,  0.00994845,
                                                       0.00413116,  0.04555686, -0.00431862, -0.0138286 , -0.04412473,
                                                       0.02016718, -0.03868735,  0.03842166, -0.0009337 ,  0.02977617,
                                                       0.02310445,  0.02827616 };
    variableVector previousPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticGradientMicroDeformation = variableVector( 27, 0 );
    variableVector previousPlasticMacroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroGradientVelocityGradient = variableVector( 27, 0 );

    variableVector currentElasticDeformationGradient = currentDeformationGradient;
    variableVector currentElasticMicroDeformation = currentMicroDeformation;
    variableVector currentElasticGradientMicroDeformation = currentGradientMicroDeformation;
    variableVector currentPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticGradientMicroDeformation = variableVector( 27, 0 );

    variableType   previousdMacroGdMacroCohesion = 0;
    variableType   previousdMicroGdMicroCohesion = 0;
    variableMatrix previousdMicroGradientGdMicroGradientCohesion( 3, variableVector( 3, 0 ) );

    variableType currentMacroStrainISV = 0;
    variableType currentMicroStrainISV = 0;
    variableVector currentMicroGradientStrainISV( 3, 0 );

    BOOST_CHECK( !error );

    solverTools::floatMatrix floatArgsDefault =
        {
            { Dt },
            currentDeformationGradient,
            currentMicroDeformation,
            currentGradientMicroDeformation,
            { currentMacroGamma },
            { currentMicroGamma },
            currentMicroGradientGamma,
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
//            macroHardeningParameters,
//            microHardeningParameters,
//            microGradientHardeningParameters,
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

    variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

    solverTools::floatMatrix floatOutsDefault = 
        {
            currentPK2Stress,
            currentReferenceMicroStress,
            currentReferenceHigherOrderStress,
            { currentMacroStrainISV },
            { currentMicroStrainISV },
            currentMicroGradientStrainISV,
            {},
            {},
            {},
            {}
        };

    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intArgs = { { 0 } };
    solverTools::intMatrix intOutsDefault = { };
    
    solverTools::intMatrix intOuts = intOutsDefault;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG;
#endif

    solverTools::floatVector x( 45, 0 );
    for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
        x[ i + 0 ] = currentPlasticDeformationGradient[ i ];
        x[ i + 9 ] = currentPlasticMicroDeformation[ i ];
    }
    for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
        x[ i + 18 ] = currentPlasticGradientMicroDeformation[ i ];
    }
//    x[ 45 ] = currentMacroGamma;
//    x[ 46 ] = currentMicroGamma;
//    for ( unsigned int i = 0; i < currentMicroGradientGamma.size(); i++ ){
//        x[ 47 + i ] = currentMicroGradientGamma[ i ];
//    }

    solverTools::floatVector residual;
    solverTools::floatMatrix jacobian;
    solverTools::floatVector residualAnswer =
    {
       -0.01815695,  0.01501355, -0.00834171,  0.0143829 ,  0.00306498,
       -0.02539779, -0.0152216 , -0.01990943, -0.00694793,  0.00268528,
        0.00234368, -0.0228185 ,  0.01557574,  0.00628433, -0.00204614,
        0.00363844, -0.02058245,  0.00657132, -0.00085036, -0.01079597,
        0.00694313, -0.00685743,  0.01366363,  0.00825854,  0.01981985,
        0.0137994 ,  0.01596549, -0.00819879, -0.00149723,  0.00232586,
       -0.00973346, -0.01658729,  0.00511054, -0.01596974,  0.02714507,
        0.01888991,  0.00869238,  0.00112477,  0.01152114, -0.01415901,
        0.02438177,  0.00265287,  0.0312893 ,  0.01743839,  0.00304193
    }; 

    error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                             floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                             , DEBUG
#endif
                                                                             );

#ifdef DEBUG_MODE
    //Check if the floatOuts are what is expected.
    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 0 ], DEBUG[ "currentPK2Stress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 1 ], DEBUG[ "currentReferenceMicroStress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 2 ], DEBUG[ "currentReferenceHigherOrderStress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 3 ], DEBUG[ "currentMacroStrainISV" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 4 ], DEBUG[ "currentMicroStrainISV" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 5 ], DEBUG[ "currentMicroGradientStrainISV" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 6 ], DEBUG[ "currentMacroCohesion" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 7 ], DEBUG[ "currentMicroCohesion" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 8 ], DEBUG[ "currentMicroGradientCohesion" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 9 ], DEBUG[ "currentElasticRightCauchyGreen" ] ) );
#endif

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( residualAnswer, residual ) );

    //Test the plastic deformation jacobians
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < x.size(); i++ ){
        constantVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector residual_P, residual_M;
        solverTools::floatMatrix jacobian_P, jacobian_M;

        solverTools::floatMatrix fA, fO;
        solverTools::intMatrix iA, iO;

        fA = floatArgsDefault;
        iA = intArgs;

#ifdef DEBUG_MODE
        solverTools::debugMap DEBUG_P, DEBUG_M;
#endif

        fO = floatOutsDefault;
        iO = intOutsDefault;

        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x + delta, fA, iA, residual_P, jacobian_P,
                                                                                 fO, iO
#ifdef DEBUG_MODE
                                                                                 , DEBUG_P
#endif
                                                                               );

        BOOST_CHECK( !error );

        fO = floatOutsDefault;
        iO = intOutsDefault;

        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x - delta, fA, iA, residual_M, jacobian_M,
                                                                                 fO, iO
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE

        //Debug each of the sub-jacobians if required. This can be very slow so it isn't done all the time

        //Assemble the numeric Jacobians
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P.begin(); it_P != DEBUG_P.end(); it_P++ ){
            
            auto it_M = DEBUG_M.find( it_P->first );
            BOOST_CHECK( it_M != DEBUG_M.end() );

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );
        }

        if ( i < 9 ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                DEBUG[ "dElasticDeformationGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticRightCauchyGreendPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) );
            }

            /*!=================
            | Strain-Like ISVS |
            ==================*/
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!================
            | Cohesion values |
            =================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-8, 1e-8 ) );
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-8 ) );
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

//            /*!====================
//            | The yield equations |
//            =====================*/
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMacroYielddPlasticDeformationGradient" ][ i ],
//                                            1e-5, 1e-8 ) );
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMicroYielddPlasticDeformationGradient" ][ i ],
//                                            1e-5, 1e-8 ) );
//
//            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroGradientYielddPlasticDeformationGradient" ][ 9 * j + i ],
//                                                1e-5, 1e-8 ) );
//            }
        }

        if ( ( i >= 9 ) && ( i < 18 ) ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticMicroRightCauchyGreendPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) );
            }

            /*!=================
            | Strain-Like ISVS |
            ==================*/
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!================
            | Cohesion values |
            =================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-8, 1e-8 ) );
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-8 ) );
            }

//            /*!====================
//            | The yield equations |
//            =====================*/
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMacroYielddPlasticMicroDeformation" ][ i - 9 ],
//                                            1e-5, 1e-8 ) );
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMicroYielddPlasticMicroDeformation" ][ i - 9 ],
//                                            1e-5, 1e-8 ) );
//
//            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroGradientYielddPlasticMicroDeformation" ][ 9 * j + i - 9 ],
//                                                1e-5, 1e-7 ) );
//            }
        }

        if ( ( i >= 18 ) && ( i < 45 ) ){
            
            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) );
            }

            /*!=================
            | Strain-Like ISVS |
            ==================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!================
            | Cohesion values |
            =================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) );
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-7 ) );
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-8 ) );
            }

//            /*!====================
//            | The yield equations |
//            =====================*/
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMacroYielddPlasticGradientMicroDeformation" ][ i - 18 ],
//                                            1e-5, 1e-8 ) );
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMicroYielddPlasticGradientMicroDeformation" ][ i - 18 ],
//                                            1e-5, 1e-8 ) );
//
//            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroGradientYielddPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
//                                                1e-5, 1e-8 ) );
//            }
        }
//
//        //Jacobians w.r.t. the plastic multipliers
//        if ( ( i >= 45 ) && ( i < 50 ) ){
//
//            /*==============================
//            | Elastic Deformation Measures |
//            ==============================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            /*!=====================================
//            | Elastic Derived Deformation Measures |
//            ======================================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            /*!================
//            | Stress Measures |
//            =================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
//                                                0.,
//                                                1e-6, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
//                                                0.,
//                                                1e-6, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
//                                                0.,
//                                                1e-6, 1e-9 ) );
//            }
//
//            /*!=================
//            | Strain-Like ISVS |
//            ==================*/
//
//            if ( i == 45 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
//                                                DEBUG[ "dCurrentMacroStrainISVdMacroGamma" ],
//                                                1e-9, 1e-9 ) );
//
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
//                                                variableVector( 3, 0 ),
//                                                1e-9, 1e-9 ) );
//            }
//
//            if ( i == 46 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
//                                                DEBUG[ "dCurrentMicroStrainISVdMicroGamma" ],
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
//                                                variableVector( 3, 0 ),
//                                                1e-9, 1e-9 ) );
//            }
//
//            if ( i >= 47 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                for ( unsigned int j = 0; j < 3; j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ][ j ],
//                                                    DEBUG[ "dCurrentMicroGradientStrainISVdMicroGradientGamma" ][ 3 * j + i - 47 ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            /*!================
//            | Cohesion values |
//            =================*/
//
//            if ( i == 45 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
//                                                DEBUG[ "dCurrentMacroCohesiondMacroGamma" ],
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
//                                                variableVector( 3, 0 ),
//                                                1e-9, 1e-9 ) );
//            }
//
//            if ( i == 46 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
//                                                DEBUG[ "dCurrentMicroCohesiondMicroGamma" ],
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
//                                                variableVector( 3, 0 ),
//                                                1e-9, 1e-9 ) );
//            }
//
//            if ( i >= 47 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                for ( unsigned int j = 0; j < 3; j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ][ j ],
//                                                    DEBUG[ "dCurrentMicroGradientCohesiondMicroGradientGamma" ][ 3 * j + i - 47 ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            /*=================
//            | Flow Directions |
//            =================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            /*!===========================
//            | Plastic velocity gradients |
//            ============================*/
//
//            if ( i == 45 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMacroVelocityGradientdMacroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            if ( i == 46 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMacroVelocityGradientdMicroGamma" ][ j ],
//                                                    1e-9, 1e-8 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMicroVelocityGradientdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            if ( i >= 47 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGradientGamma" ][ 3 * j + i - 47 ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            /*!======================================
//            | Expected plastic deformation measures |
//            =======================================*/
//
//            if ( i == 45 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticDeformationGradientdMacroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMacroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            else if ( i == 46 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticDeformationGradientdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticMicroDeformationdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            else if ( i >= 47 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGradientGamma" ][ 3 * j + i - 47 ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            /*!====================
//            | The yield equations |
//            =====================*/
//
//            if ( i == 45 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
//                                                DEBUG[ "dMacroYielddMacroGamma" ],
//                                                1e-5, 1e-8 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
//                                                variableVector( 1, 0 ),
//                                                1e-5, 1e-8 ) );
//    
//                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                    0.,
//                                                    1e-5, 1e-8 ) );
//                }
//            }
//            else if ( i == 46 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
//                                                variableVector( 1, 0. ),
//                                                1e-5, 1e-8 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
//                                                DEBUG[ "dMicroYielddMicroGamma" ],
//                                                1e-5, 1e-8 ) );
//    
//                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                    0.,
//                                                    1e-5, 1e-8 ) );
//                }
//            }
//            else if ( i >= 47 ){
//
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
//                                                variableVector( 1, 0. ),
//                                                1e-5, 1e-8 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
//                                                variableVector( 1, 0. ),
//                                                1e-5, 1e-8 ) );
//    
//                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                    DEBUG[ "dMicroGradientYielddMicroGradientGamma" ][ 3 * ( i - 47 ) + j ],
//                                                    1e-5, 1e-8 ) );
//                }
//            } 
//        }
//
//        if ( i >= 50 ){
//            std::cout << "ERROR in the test!\n";
//            assert( 1 == 0 );
//        }

#endif

        //Jacobian test

        solverTools::floatVector gradCol = ( residual_P - residual_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ], 1e-5, 1e-7 ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticDeformationResidual2 ){
    /*!
     * Second test the computation of the plastic deformation residual.
     *
     */

    std::vector< double > fparams = { 2, 1.0e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.5e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e2, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    parameterVector macroHardeningParameters;
    parameterVector microHardeningParameters;
    parameterVector microGradientHardeningParameters;

    parameterVector macroFlowParameters;
    parameterVector microFlowParameters;
    parameterVector microGradientFlowParameters;

    parameterVector macroYieldParameters;
    parameterVector microYieldParameters;
    parameterVector microGradientYieldParameters;

    parameterVector Amatrix, Bmatrix, Cmatrix, Dmatrix;
    parameterType alphaMacro, alphaMicro, alphaMicroGradient;
    parameterType relativeTolerance, absoluteTolerance;

    errorOut error = micromorphicElastoPlasticity::extractMaterialParameters( fparams,
                                                                              macroHardeningParameters, microHardeningParameters,
                                                                              microGradientHardeningParameters,
                                                                              macroFlowParameters, microFlowParameters,
                                                                              microGradientFlowParameters,
                                                                              macroYieldParameters, microYieldParameters,
                                                                              microGradientYieldParameters,
                                                                              Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                              alphaMacro, alphaMicro, alphaMicroGradient,
                                                                              relativeTolerance, absoluteTolerance );

    BOOST_CHECK( !error );

    constantType   Dt = 2.5;
    variableType   currentMacroGamma = 0.01;
    variableType   currentMicroGamma = 0.02;
    variableVector currentMicroGradientGamma = {0.011, 0.021, 0.031};
    variableType   previousMacroGamma = 0.;
    variableType   previousMicroGamma = 0.;
    variableVector previousMicroGradientGamma( 3, 0 );
    variableType   currentMacroStrainISV = 0.;
    variableType   currentMicroStrainISV = 0.;
    variableVector currentMicroGradientStrainISV( 3., 0. );
    variableType   previousMacroStrainISV = 0.;
    variableType   previousMicroStrainISV = 0.;
    variableVector previousMicroGradientStrainISV( 3., 0. );
    variableVector currentDeformationGradient = { 0.04482969,  0.88562312, -0.38710144,
                                                 -0.93722716,  0.19666568,  0.41677155,
                                                  0.46929057,  0.33779672,  0.81392228 };
    variableVector currentMicroDeformation = {  0.51930689,  0.27954023, -0.85955731,
                                                0.09469279, -0.99381243, -0.23218079,
                                               -0.82281393,  0.09643296, -0.54637704 };
    variableVector currentGradientMicroDeformation = { 0.04176306, -0.0151958 , -0.00090558, -0.01844751,  0.04512391,
                                                       0.02174263, -0.00508239, -0.01827377,  0.00541031, -0.01330239,
                                                      -0.02479987,  0.02914825,  0.00168841,  0.00230506,  0.00994845,
                                                       0.00413116,  0.04555686, -0.00431862, -0.0138286 , -0.04412473,
                                                       0.02016718, -0.03868735,  0.03842166, -0.0009337 ,  0.02977617,
                                                       0.02310445,  0.02827616 };
    variableVector previousPlasticDeformationGradient = { -0.02449435, -0.57317707,  0.81909948,
                                                          -0.95484326, -0.21673221, -0.17800334,
                                                           0.28314298, -0.78904695, -0.55302402 };
    variableVector previousPlasticMicroDeformation = { 0.49143847, -0.61330346, -0.61490607,
                                                       0.17446661,  0.7585011 , -0.62811593,
                                                       0.85044752,  0.20128837,  0.48228538 };
    variableVector previousPlasticGradientMicroDeformation = {  1.17975391e-03, -4.91057046e-03, -1.97589125e-03,  4.28967411e-03,
                                                                7.47454477e-05, -3.05112255e-03, -1.37754405e-03,  1.84975854e-03,
                                                               -2.79115642e-03, -3.79286235e-03,  9.99338104e-04, -4.74822369e-03,
                                                               -4.45479780e-03,  1.71442686e-03,  4.86890092e-03,  9.74233152e-04,
                                                                3.42092901e-03, -9.45725029e-04,  3.52072403e-03, -1.58992133e-03,
                                                               -9.98608169e-04,  4.06025431e-03,  1.69580352e-03, -8.54127175e-04,
                                                               -2.81532605e-03, -4.68166468e-03, -1.28948806e-03 };
    variableVector previousPlasticMacroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroGradientVelocityGradient = variableVector( 27, 0 );

    variableVector currentElasticDeformationGradient = currentDeformationGradient;
    variableVector currentElasticMicroDeformation = currentMicroDeformation;
    variableVector currentElasticGradientMicroDeformation = currentGradientMicroDeformation;
    variableVector currentPlasticDeformationGradient = { -0.02449435, -0.57317707,  0.81909948,
                                                         -0.95484326, -0.21673221, -0.17800334,
                                                          0.28314298, -0.78904695, -0.55302402 };
    variableVector currentPlasticMicroDeformation = { 0.49143847, -0.61330346, -0.61490607,
                                                      0.17446661,  0.7585011 , -0.62811593,
                                                      0.85044752,  0.20128837,  0.48228538 };
    variableVector currentPlasticGradientMicroDeformation = {  1.17975391e-03, -4.91057046e-03, -1.97589125e-03,  4.28967411e-03,
                                                               7.47454477e-05, -3.05112255e-03, -1.37754405e-03,  1.84975854e-03,
                                                              -2.79115642e-03, -3.79286235e-03,  9.99338104e-04, -4.74822369e-03,
                                                              -4.45479780e-03,  1.71442686e-03,  4.86890092e-03,  9.74233152e-04,
                                                               3.42092901e-03, -9.45725029e-04,  3.52072403e-03, -1.58992133e-03,
                                                              -9.98608169e-04,  4.06025431e-03,  1.69580352e-03, -8.54127175e-04,
                                                              -2.81532605e-03, -4.68166468e-03, -1.28948806e-03 };

    variableType   previousdMacroGdMacroCohesion = 0;
    variableType   previousdMicroGdMicroCohesion = 0;
    variableMatrix previousdMicroGradientGdMicroGradientCohesion( 3, variableVector( 3, 0 ) );

    BOOST_CHECK( !error );

    solverTools::floatMatrix floatArgsDefault =
        {
            { Dt },
            currentDeformationGradient,
            currentMicroDeformation,
            currentGradientMicroDeformation,
            { currentMacroGamma },
            { currentMicroGamma },
            currentMicroGradientGamma,
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
//            macroHardeningParameters,
//            microHardeningParameters,
//            microGradientHardeningParameters,
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

    variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

    solverTools::floatMatrix floatOutsDefault = 
        {
            currentPK2Stress,
            currentReferenceMicroStress,
            currentReferenceHigherOrderStress,
            { currentMacroStrainISV },
            { currentMicroStrainISV },
            currentMicroGradientStrainISV,
            {}, {}, {},
            {},
            {}, {}, {}, {}, {}, {}, {}
        };

    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intArgs = { { 1 } };
    solverTools::intMatrix intOutsDefault = { };
    
    solverTools::intMatrix intOuts = intOutsDefault;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG;
#endif

    solverTools::floatVector x( 45, 0 );
    for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
        x[ i + 0 ] = currentPlasticDeformationGradient[ i ];
        x[ i + 9 ] = currentPlasticMicroDeformation[ i ];
    }
    for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
        x[ i + 18 ] = currentPlasticGradientMicroDeformation[ i ];
    }

    solverTools::floatVector residual;
    solverTools::floatMatrix jacobian;
    solverTools::floatVector residualAnswer =
        {
             0.00240362,  0.0248923,   -0.000583893,  0.0265404,   -0.000308798,
             0.000169618, 0.00397561,  -0.00324285,  -0.0273723,   -0.00330406,
             0.00550047, -0.0170485,   -0.0129665,    0.0095253,    0.00738336,
             0.0135288,   0.0172144,    0.00998783,   0.0147325,   -0.00133292,
             0.0252213,  -0.000932423, -0.00783062,   0.0069921,    0.0118855,
             0.0166555,  -0.00558663,  -0.00134957,   0.037993,     0.000269495,
             0.00305745,  0.00555468,   0.024898,    -0.00886859,   0.0112834,
            -0.0120454,   0.0339332,    0.00802589,   0.000475772,  0.00990966,
             0.00691379, -0.00383809,  -0.000962487,  0.000727374, -0.00245113
        }; 

    error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                             floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                             , DEBUG
#endif
                                                                             );

    BOOST_CHECK( !error );

#ifdef DEBUG_MODE
    //Check if the floatOuts are what is expected.
    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 0 ], DEBUG[ "currentPK2Stress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 1 ], DEBUG[ "currentReferenceMicroStress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 2 ], DEBUG[ "currentReferenceHigherOrderStress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 3 ], DEBUG[ "currentMacroStrainISV" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 4 ], DEBUG[ "currentMicroStrainISV" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 5 ], DEBUG[ "currentMicroGradientStrainISV" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 6 ], DEBUG[ "currentMacroCohesion" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 7 ], DEBUG[ "currentMicroCohesion" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 8 ], DEBUG[ "currentMicroGradientCohesion" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 9 ], DEBUG[ "currentElasticRightCauchyGreen" ] ) );
#endif

    BOOST_CHECK( vectorTools::fuzzyEquals( residualAnswer, residual, 1e-5 ) );

    //Test the plastic deformation jacobians
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < x.size(); i++ ){
        constantVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector residual_P, residual_M;
        solverTools::floatMatrix jacobian_P, jacobian_M;

        solverTools::floatMatrix fA, fO_P, fO_M;
        solverTools::intMatrix iA, iO_P, iO_M;

        fA = floatArgsDefault;
        iA = intArgs;

#ifdef DEBUG_MODE
        solverTools::debugMap DEBUG_P, DEBUG_M;
#endif

        fO_P = floatOutsDefault;
        fO_M = floatOutsDefault;

        iO_P = intOutsDefault;
        iO_M = intOutsDefault;

        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x + delta, fA, iA, residual_P, jacobian_P,
                                                                                 fO_P, iO_P
#ifdef DEBUG_MODE
                                                                                 , DEBUG_P
#endif
                                                                               );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x - delta, fA, iA, residual_M, jacobian_M,
                                                                                 fO_M, iO_M
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE

        //Debug each of the sub-jacobians if required. This can be very slow so it isn't done all the time

        //Assemble the numeric Jacobians
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P.begin(); it_P != DEBUG_P.end(); it_P++ ){
            
            auto it_M = DEBUG_M.find( it_P->first );
            BOOST_CHECK( it_M != DEBUG_M.end() );

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );
        }

        if ( i < 9 ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                DEBUG[ "dElasticDeformationGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticRightCauchyGreendPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) );
            }

            /*!=================
            | Strain-Like ISVS |
            ==================*/
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!================
            | Cohesion values |
            =================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-8, 1e-8 ) );
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-8 ) );
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

//            /*!====================
//            | The yield equations |
//            =====================*/
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMacroYielddPlasticDeformationGradient" ][ i ],
//                                            1e-5, 1e-8 ) );
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMicroYielddPlasticDeformationGradient" ][ i ],
//                                            1e-5, 1e-8 ) );
//
//            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroGradientYielddPlasticDeformationGradient" ][ 9 * j + i ],
//                                                1e-5, 1e-8 ) );
//            }
        }

        if ( ( i >= 9 ) && ( i < 18 ) ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticMicroRightCauchyGreendPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) );
            }

            /*!=================
            | Strain-Like ISVS |
            ==================*/
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!================
            | Cohesion values |
            =================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-8, 1e-8 ) );
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-8 ) );
            }

//            /*!====================
//            | The yield equations |
//            =====================*/
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMacroYielddPlasticMicroDeformation" ][ i - 9 ],
//                                            1e-5, 1e-8 ) );
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMicroYielddPlasticMicroDeformation" ][ i - 9 ],
//                                            1e-5, 1e-8 ) );
//
//            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroGradientYielddPlasticMicroDeformation" ][ 9 * j + i - 9 ],
//                                                1e-5, 1e-7 ) );
//            }
        }

        if ( ( i >= 18 ) && ( i < 45 ) ){
            
            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-7 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-7 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-7 ) );
            }

            /*!=================
            | Strain-Like ISVS |
            ==================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!================
            | Cohesion values |
            =================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) );
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-7 ) );
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-8 ) );
            }

//            /*!====================
//            | The yield equations |
//            =====================*/
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMacroYielddPlasticGradientMicroDeformation" ][ i - 18 ],
//                                            1e-5, 1e-7 ) );
//
//            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
//                                            DEBUG[ "dMicroYielddPlasticGradientMicroDeformation" ][ i - 18 ],
//                                            1e-5, 1e-7 ) );
//
//            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroGradientYielddPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
//                                                1e-5, 1e-7 ) );
//            }
        }

//        //Jacobians w.r.t. the plastic multipliers
//        if ( ( i >= 45 ) && ( i < 50 ) ){
//
//            /*==============================
//            | Elastic Deformation Measures |
//            ==============================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            /*!=====================================
//            | Elastic Derived Deformation Measures |
//            ======================================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            /*!================
//            | Stress Measures |
//            =================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
//                                                0.,
//                                                1e-6, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
//                                                0.,
//                                                1e-6, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
//                                                0.,
//                                                1e-6, 1e-9 ) );
//            }
//
//            /*!=================
//            | Strain-Like ISVS |
//            ==================*/
//
//            if ( i == 45 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
//                                                DEBUG[ "dCurrentMacroStrainISVdMacroGamma" ],
//                                                1e-9, 1e-9 ) );
//
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
//                                                variableVector( 3, 0 ),
//                                                1e-9, 1e-9 ) );
//            }
//
//            if ( i == 46 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
//                                                DEBUG[ "dCurrentMicroStrainISVdMicroGamma" ],
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
//                                                variableVector( 3, 0 ),
//                                                1e-9, 1e-9 ) );
//            }
//
//            if ( i >= 47 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                for ( unsigned int j = 0; j < 3; j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ][ j ],
//                                                    DEBUG[ "dCurrentMicroGradientStrainISVdMicroGradientGamma" ][ 3 * j + i - 47 ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            /*!================
//            | Cohesion values |
//            =================*/
//
//            if ( i == 45 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
//                                                DEBUG[ "dCurrentMacroCohesiondMacroGamma" ],
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
//                                                variableVector( 3, 0 ),
//                                                1e-9, 1e-9 ) );
//            }
//
//            if ( i == 46 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
//                                                DEBUG[ "dCurrentMicroCohesiondMicroGamma" ],
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
//                                                variableVector( 3, 0 ),
//                                                1e-9, 1e-9 ) );
//            }
//
//            if ( i >= 47 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
//                                                variableVector( 1, 0 ),
//                                                1e-9, 1e-9 ) );
//    
//                for ( unsigned int j = 0; j < 3; j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ][ j ],
//                                                    DEBUG[ "dCurrentMicroGradientCohesiondMicroGradientGamma" ][ 3 * j + i - 47 ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            /*=================
//            | Flow Directions |
//            =================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
//                                                0.,
//                                                1e-9, 1e-9 ) );
//            }
//
//            /*!===========================
//            | Plastic velocity gradients |
//            ============================*/
//
//            if ( i == 45 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMacroVelocityGradientdMacroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            if ( i == 46 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMacroVelocityGradientdMicroGamma" ][ j ],
//                                                    1e-9, 1e-8 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMicroVelocityGradientdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            if ( i >= 47 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
//                                                    DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGradientGamma" ][ 3 * j + i - 47 ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            /*!======================================
//            | Expected plastic deformation measures |
//            =======================================*/
//
//            if ( i == 45 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticDeformationGradientdMacroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMacroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            else if ( i == 46 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticDeformationGradientdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticMicroDeformationdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGamma" ][ j ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            else if ( i >= 47 ){
//
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
//                                                    0.,
//                                                    1e-9, 1e-9 ) );
//                }
//    
//                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
//                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGradientGamma" ][ 3 * j + i - 47 ],
//                                                    1e-9, 1e-9 ) );
//                }
//            }
//
//            /*!====================
//            | The yield equations |
//            =====================*/
//
//            if ( i == 45 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
//                                                DEBUG[ "dMacroYielddMacroGamma" ],
//                                                1e-5, 1e-8 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
//                                                variableVector( 1, 0 ),
//                                                1e-5, 1e-8 ) );
//    
//                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                    0.,
//                                                    1e-5, 1e-8 ) );
//                }
//            }
//            else if ( i == 46 ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
//                                                variableVector( 1, 0. ),
//                                                1e-5, 1e-8 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
//                                                DEBUG[ "dMicroYielddMicroGamma" ],
//                                                1e-5, 1e-8 ) );
//    
//                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                    0.,
//                                                    1e-5, 1e-8 ) );
//                }
//            }
//            else if ( i >= 47 ){
//
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
//                                                variableVector( 1, 0. ),
//                                                1e-5, 1e-8 ) );
//    
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
//                                                variableVector( 1, 0. ),
//                                                1e-5, 1e-8 ) );
//    
//                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                    BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                    DEBUG[ "dMicroGradientYielddMicroGradientGamma" ][ 3 * ( i - 47 ) + j ],
//                                                    1e-5, 1e-8 ) );
//                }
//            } 
//        }
//
//        if ( i >= 50 ){
//            std::cout << "ERROR in the test!\n";
//            assert( 1 == 0 );
//        }

#endif

        //Jacobian test

        solverTools::floatVector gradCol = ( residual_P - residual_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ], 1e-5, 1e-7 ) );
        }

        //Test the partial derivative of the stress measures w.r.t. the plastic fundamental deformation measures
        gradCol = vectorTools::appendVectors( { fO_P[ 0 ] - fO_M[ 0 ], fO_P[ 1 ] - fO_M[ 1 ], fO_P[ 2 ] - fO_M[ 2 ] } ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 11 ][ 45 * j + i ], 1e-5 ) );
        }

        //Test the partial derivative of the elastic Right Cauchy Green deformation tensor w.r.t. the plastic deformation gradient
        if ( i < 9 ){
            gradCol = ( fO_P[ 9 ] - fO_M[ 9 ] ) / ( 2 * delta[ i ] );
            for ( unsigned int j = 0; j < gradCol.size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 10 ][ 9 * j + i ] ) );
            }
        }
    }

    //Test the jacobians w.r.t. the deformation measures
    //construct the jacobian matrix w.r.t. the fundamental deformation measures
    solverTools::floatVector baseDeformationVector = vectorTools::appendVectors( { currentDeformationGradient,
                                                                                   currentMicroDeformation,
                                                                                   currentGradientMicroDeformation } );
    solverTools::floatMatrix deformationJacobian( baseDeformationVector.size(),
                                                  solverTools::floatVector( baseDeformationVector.size(), 0 ) );

    for ( unsigned int i = 0; i < baseDeformationVector.size(); i++ ){
        constantVector delta( baseDeformationVector.size(), 0 );
        delta[ i ] = eps * fabs( baseDeformationVector[ i ] ) + eps;

        solverTools::floatVector perturbedVector_P = baseDeformationVector + delta;
        solverTools::floatVector perturbedVector_M = baseDeformationVector - delta;

        //Extract the perturbed deformation vectors
        solverTools::floatVector cF_P( perturbedVector_P.begin(), perturbedVector_P.begin() + 9 );
        solverTools::floatVector cF_M( perturbedVector_M.begin(), perturbedVector_M.begin() + 9 );

        //Extract the perturbed micro deformation vectors
        solverTools::floatVector cChi_P( perturbedVector_P.begin() + 9, perturbedVector_P.begin() + 18 );
        solverTools::floatVector cChi_M( perturbedVector_M.begin() + 9, perturbedVector_M.begin() + 18 );

        //Extract the perturbed gradient micro deformation
        solverTools::floatVector cGradChi_P( perturbedVector_P.begin() + 18, perturbedVector_P.begin() + 45 );
        solverTools::floatVector cGradChi_M( perturbedVector_M.begin() + 18, perturbedVector_M.begin() + 45 );

        //Assemble the floatArgs matrices
        solverTools::floatMatrix floatArgs_P = floatArgsDefault;
        solverTools::floatMatrix floatArgs_M = floatArgsDefault;

        floatArgs_P[ 1 ] = cF_P;
        floatArgs_P[ 2 ] = cChi_P;
        floatArgs_P[ 3 ] = cGradChi_P;

        floatArgs_M[ 1 ] = cF_M;
        floatArgs_M[ 2 ] = cChi_M;
        floatArgs_M[ 3 ] = cGradChi_M;

        //Evaluate the residual
        solverTools::floatMatrix fO_P, fO_M;
        fO_P = floatOutsDefault;
        fO_M = floatOutsDefault;

        solverTools::intMatrix iO_P, iO_M;
        iO_P = intOutsDefault;
        iO_M = intOutsDefault;

        solverTools::floatVector residual_P, residual_M;
        solverTools::floatMatrix _J;

#ifdef DEBUG_MODE
        solverTools::debugMap DEBUG_P, DEBUG_M;
#endif

        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x, floatArgs_P, intArgs, residual_P, _J,
                                                                                 fO_P, iO_P
#ifdef DEBUG_MODE
                                                                                 , DEBUG_P
#endif
                                                                               );

        BOOST_CHECK( !error );
        
        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x, floatArgs_M, intArgs, residual_M, _J,
                                                                                 fO_M, iO_M
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        BOOST_CHECK( !error );

        //Debug each of the sub-jacobians if required. This can be very slow so it isn't done all the time

#ifdef DEBUG_MODE

        //Assemble the numeric Jacobians w.r.t. the deformation parameters
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P.begin(); it_P != DEBUG_P.end(); it_P++ ){
            
            auto it_M = DEBUG_M.find( it_P->first );
            BOOST_CHECK( it_M != DEBUG_M.end() );

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );
        }

        //Debug the Jacobians w.r.t. the current deformation gradient
        if ( i < 9 ) {

            /*==========================
            | The Elastic Deformations |
            ==========================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                DEBUG[ "dElasticDeformationGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            /*!=========================================
            | The Elastic Derived Deformation Measures |
            ==========================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticRightCauchyGreendDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-9 ) );
            }

            /*!=================
            | Strain-Like ISVS |
            ==================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!==============
            | The Cohesions |
            ===============*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            solverTools::floatVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!====================
            | The Flow Directions |
            =====================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondDeformationGradient" ][ 9 * j + i ],
                                                1e-8, 1e-8 ) );
            }

            /*!===========================
            | Plastic Velocity Gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroGradientVelocityGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            /*!==============================
            | Plastic Deformation Gradients |
            ===============================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) );
            }

//            /*!====================
//            | The Yield Equations |
//            =====================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "macroYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ j ],
//                                                DEBUG[ "dMacroYielddDeformationGradient" ][ 9 * j + i ],
//                                                1e-6, 1e-8 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "microYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroYielddDeformationGradient" ][ 9 * j + i ],
//                                                1e-6, 1e-8 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroGradientYielddDeformationGradient" ][ 9 * j + i ],
//                                                1e-6, 1e-8 ) );
//            }
        }

        if ( ( i >= 9 ) && ( i < 18 ) ){

            /*==========================
            | The Elastic Deformations |
            ==========================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticMicroDeformationdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!=========================================
            | The Elastic Derived Deformation Measures |
            ==========================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticMicroRightCauchyGreendMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-9 ) );
            }

            /*!=================
            | Strain-Like ISVS |
            ==================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!==============
            | The Cohesions |
            ===============*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            solverTools::floatVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!====================
            | The Flow Directions |
            =====================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!===========================
            | Plastic Velocity Gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            /*!==============================
            | Plastic Deformation Gradients |
            ===============================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) );
            }

//            /*!====================
//            | The Yield Equations |
//            =====================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "macroYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ j ],
//                                                DEBUG[ "dMacroYielddMicroDeformation" ][ 9 * j + i - 9 ],
//                                                1e-6, 1e-8 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "microYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroYielddMicroDeformation" ][ 9 * j + i - 9 ],
//                                                1e-6, 1e-8 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroGradientYielddMicroDeformation" ][ 9 * j + i - 9 ],
//                                                1e-6, 1e-7 ) );
//            }
        }

        if ( ( i >= 18 ) && ( i < 45 ) ){

            /*==========================
            | The Elastic Deformations |
            ==========================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            /*!=========================================
            | The Elastic Derived Deformation Measures |
            ==========================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) );
            }

            /*!=================
            | Strain-Like ISVS |
            ==================*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!==============
            | The Cohesions |
            ===============*/

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            solverTools::floatVector( 3, 0 ),
                                            1e-9, 1e-9 ) );

            /*!====================
            | The Flow Directions |
            =====================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-8, 1e-7 ) );
            }

            /*!===========================
            | Plastic Velocity Gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroGradientVelocityGradientdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            /*!==============================
            | Plastic Deformation Gradients |
            ===============================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) );
            }

//            /*!====================
//            | The Yield Equations |
//            =====================*/
//
//            for ( unsigned int j = 0; j < numericGradients[ "macroYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ j ],
//                                                DEBUG[ "dMacroYielddGradientMicroDeformation" ][ 27 * j + i - 18 ],
//                                                1e-6, 1e-7 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "microYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroYielddGradientMicroDeformation" ][ 27 * j + i - 18 ],
//                                                1e-6, 1e-7 ) );
//            }
//
//            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
//                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
//                                                DEBUG[ "dMicroGradientYielddGradientMicroDeformation" ][ 27 * j + i - 18 ],
//                                                1e-6, 1e-7 ) );
//            }
        }

        if ( i >= 45 ){
            assert( 1 == 0 );
        }

#endif

        //Test the partial derivative of the residual w.r.t. the fundamental deformation measures 
        solverTools::floatVector gradCol = ( residual_P - residual_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 13 ][ 45 * j + i ] ) );
        }

        //Test the partial derivative of the stress measures w.r.t. the fundamental deformation measures
        gradCol = vectorTools::appendVectors( { fO_P[ 0 ] - fO_M[ 0 ], fO_P[ 1 ] - fO_M[ 1 ], fO_P[ 2 ] - fO_M[ 2 ] } ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 12 ][ 45 * j + i ], 1e-5 ) );
        }

        //Test the partial derivative of the Elastic Right Cauchy Green deformation tensor w.r.t. the deformation gradient
        if ( i < 9 ) {
            gradCol = ( fO_P[ 9 ] - fO_M[ 9 ] ) / ( 2 * delta[ i ] );
            for ( unsigned int j = 0; j < gradCol.size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 14 ][ 9 * j + i ] ) );
            }
        }
    }

    //Test the Jacobians w.r.t. the plastic multipliers
    solverTools::floatVector baseGammaVector =
        {
            currentMacroGamma,
            currentMicroGamma,
            currentMicroGradientGamma[ 0 ],
            currentMicroGradientGamma[ 1 ],
            currentMicroGradientGamma[ 2 ]
        };

    for ( unsigned int i = 0; i < baseGammaVector.size(); i++ ){
        constantVector delta( baseGammaVector.size(), 0 );
        delta[ i ] = eps * fabs( baseGammaVector[ i ] ) + eps;

        solverTools::floatVector perturbedVector_P = baseGammaVector + delta;
        solverTools::floatVector perturbedVector_M = baseGammaVector - delta;

        //Extract the perturbed deformation vectors
        solverTools::floatVector cMacG_P( perturbedVector_P.begin(), perturbedVector_P.begin() + 1 );
        solverTools::floatVector cMacG_M( perturbedVector_M.begin(), perturbedVector_M.begin() + 1 );

        //Extract the perturbed micro deformation vectors
        solverTools::floatVector cMicG_P( perturbedVector_P.begin() + 1, perturbedVector_P.begin() + 2 );
        solverTools::floatVector cMicG_M( perturbedVector_M.begin() + 1, perturbedVector_M.begin() + 2 );

        //Extract the perturbed gradient micro deformation
        solverTools::floatVector cMicGradG_P( perturbedVector_P.begin() + 2, perturbedVector_P.begin() + 5 );
        solverTools::floatVector cMicGradG_M( perturbedVector_M.begin() + 2, perturbedVector_M.begin() + 5 );

        //Assemble the floatArgs matrices
        solverTools::floatMatrix floatArgs_P = floatArgsDefault;
        solverTools::floatMatrix floatArgs_M = floatArgsDefault;

        floatArgs_P[ 4 ] = cMacG_P;
        floatArgs_P[ 5 ] = cMicG_P;
        floatArgs_P[ 6 ] = cMicGradG_P;

        floatArgs_M[ 4 ] = cMacG_M;
        floatArgs_M[ 5 ] = cMicG_M;
        floatArgs_M[ 6 ] = cMicGradG_M;

        //Evaluate the residual
        solverTools::floatMatrix fO_P, fO_M;
        fO_P = floatOutsDefault;
        fO_M = floatOutsDefault;

        solverTools::intMatrix iO_P, iO_M;
        iO_P = intOutsDefault;
        iO_M = intOutsDefault;

        solverTools::floatVector residual_P, residual_M;
        solverTools::floatMatrix _J;

#ifdef DEBUG_MODE
        solverTools::debugMap DEBUG_P, DEBUG_M;
#endif

        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x, floatArgs_P, intArgs, residual_P, _J,
                                                                                 fO_P, iO_P
#ifdef DEBUG_MODE
                                                                                 , DEBUG_P
#endif
                                                                               );

        BOOST_CHECK( !error );
        
        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x, floatArgs_M, intArgs, residual_M, _J,
                                                                                 fO_M, iO_M
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        BOOST_CHECK( !error );

        //Debug each of the sub-jacobians if required. This can be very slow so it isn't done all the time

#ifdef DEBUG_MODE

        //Assemble the numeric Jacobians w.r.t. the deformation parameters
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P.begin(); it_P != DEBUG_P.end(); it_P++ ){
            
            auto it_M = DEBUG_M.find( it_P->first );
            BOOST_CHECK( it_M != DEBUG_M.end() );

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );
        }

        //Jacobians w.r.t. the plastic multipliers

        /*==============================
        | Elastic Deformation Measures |
        ==============================*/
  
        for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }
  
        for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }
  
        for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }

        /*!=====================================
        | Elastic Derived Deformation Measures |
        ======================================*/
  
        for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }
  
        for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }
  
        for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }
  
        for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }

        /*!================
        | Stress Measures |
        =================*/

        for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                            0.,
                                            1e-6, 1e-9 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                            0.,
                                            1e-6, 1e-9 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                            0.,
                                            1e-6, 1e-9 ) );
        }

        /*!=================
        | Strain-Like ISVS |
        ==================*/

        if ( i == 0 ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            DEBUG[ "dCurrentMacroStrainISVdMacroGamma" ],
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );
        }

        if ( i == 1 ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            DEBUG[ "dCurrentMicroStrainISVdMicroGamma" ],
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );
        }

        if ( i >= 2 ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            for ( unsigned int j = 0; j < 3; j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientStrainISV" ][ j ],
                                                DEBUG[ "dCurrentMicroGradientStrainISVdMicroGradientGamma" ][ 3 * j + i - 2 ],
                                                1e-9, 1e-9 ) );
            }
        }

        /*!================
        | Cohesion values |
        =================*/

        if ( i == 0 ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            DEBUG[ "dCurrentMacroCohesiondMacroGamma" ],
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) );
        }

        if ( i == 1 ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                                DEBUG[ "dCurrentMicroCohesiondMicroGamma" ],
                                                1e-9, 1e-9 ) );
    
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                                variableVector( 3, 0 ),
                                                1e-9, 1e-9 ) );
            }

        if ( i >= 2 ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) );

            for ( unsigned int j = 0; j < 3; j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ][ j ],
                                                DEBUG[ "dCurrentMicroGradientCohesiondMicroGradientGamma" ][ 3 * j + i - 2 ],
                                                1e-9, 1e-9 ) );
            }
        }

        /*=================
        | Flow Directions |
        =================*/

        for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                            0.,
                                            1e-9, 1e-9 ) );
        }

        /*!===========================
        | Plastic velocity gradients |
        ============================*/

        if ( i == 0 ){

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdMacroGamma" ][ j ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }
        }

        if ( i == 1 ){

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdMicroGamma" ][ j ],
                                                1e-9, 1e-8 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdMicroGamma" ][ j ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGamma" ][ j ],
                                                1e-9, 1e-9 ) );
            }
        }

        if ( i >= 2 ){

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGradientGamma" ][ 3 * j + i - 2 ],
                                                1e-9, 1e-9 ) );
            }
        }

        /*!======================================
        | Expected plastic deformation measures |
        =======================================*/

        if ( i == 0 ){

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdMacroGamma" ][ j ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdMacroGamma" ][ j ],
                                                1e-9, 1e-9 ) );
            }
        }

        else if ( i == 1 ){

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdMicroGamma" ][ j ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdMicroGamma" ][ j ],
                                                1e-9, 1e-9 ) );
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGamma" ][ j ],
                                                1e-9, 1e-9 ) );
            }
        }

        else if ( i >= 2 ){

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }
 
            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) );
            }
 
            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGradientGamma" ][ 3 * j + i - 2 ],
                                                1e-9, 1e-9 ) );
            }
        }
    
        if ( i >= 5 ){
            std::cout << "ERROR in the test!\n";
            assert( 1 == 0 );
        }
#endif

        solverTools::floatVector gradCol = ( residual_P - residual_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 15 ][ 5 * j + i ] ) );
        }

        //Check the Jacobians of the cohesions w.r.t. the Gammas
        gradCol = vectorTools::appendVectors( 
            {
                ( fO_P[ 6 ] - fO_M[ 6 ] ) / ( 2 * delta[ i ] ),
                ( fO_P[ 7 ] - fO_M[ 7 ] ) / ( 2 * delta[ i ] ),
                ( fO_P[ 8 ] - fO_M[ 8 ] ) / ( 2 * delta[ i ] )
            } );

        for ( unsigned int j = 0; j < 5; j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 16 ][ 5 * j + i ] ) );
        }

    }

}

//BOOST_AUTO_TEST_CASE( testMaterialLibraryInterface ){
//    /*!
//     * Test the interface to the linear elastic model
//     * via the material library.
//     *
//     */
//
//    //Initialize the model
//    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
//    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
//    auto material = factory.GetMaterial(_model_name);
//
//    //Set up the inputs
//    //Initialize the time
//    std::vector< double > time = { 10., 2.5 };
//
//    //Initialize the material parameters
//    std::vector< double > fparams = { 2, 2.4e2, 1.5e1,             //Macro hardening parameters
//                                      2, 1.4e2, 2.0e1,             //Micro hardening parameters
//                                      2, 2.0e0, 2.7e1,             //Micro gradient hardening parameters
//                                      2, 0.56, 0.2,                //Macro flow parameters
//                                      2, 0.15,-0.2,                //Micro flow parameters
//                                      2, 0.82, 0.1,                //Micro gradient flow parameters
//                                      2, 0.70, 0.3,                //Macro yield parameters
//                                      2, 0.40,-0.3,                //Micro yield parameters
//                                      2, 0.52, 0.4,                //Micro gradient yield parameters
//                                      2, 696.47, 65.84,            //A stiffness tensor parameters
//                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
//                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
//                                      2, -51.92, 5.13,             //D stiffness tensor parameters
//                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
//                                    };
//
//    //Initialize the gradient of the macro displacement
////    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
////                                        { -0.56632666, -0.21399259,  0.16367238 },
////                                        { -0.29129789, -0.22367825, -2.0632945  } };
////
////    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
////                                         {  0.31303067, -1.23910631, -0.93837662 },
////                                         { -0.32571524, -0.95306342, -0.93025257 } };
//
//    double current_grad_u[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
//                                        {0.100, 0.001, 0.000 },
//                                        {0.000, 0.000, 0.000 } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
//                                         {0, 0, 0},
//                                         {0, 0, 0} };
//    //Initialize the micro displacement
////    double current_phi[ 9 ] = { 0.84729289,  0.40617104,  0.59534561,  
////                                0.44195587,  0.34121966, -0.79098944, 
////                               -0.43965428,  0.88466225,  0.1684519 };
////
////    double previous_phi[ 9 ] = { -0.99935855, -0.21425717,  0.0668254 ,
////                                 -0.11111872, -0.07416114, -1.01048108,
////                                  0.1804018 , -1.01116291,  0.03248007 };
//
//    double current_phi[ 9 ] = { 0.100, 0.000, 0.000,
//                                0.000, 0.000, 0.000,
//                                0.000, 0.000, 0.000 };
//
//    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//
//    //Initialize the gradient of the micro displacement
//    double current_grad_phi[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
//                                          { -0.18572739,  0.06847269,  0.22931628 },
//                                          { -0.01829735, -0.48731265, -0.25277529 },
//                                          {  0.26626212,  0.4844646 , -0.31965177 },
//                                          {  0.49197846,  0.19051656, -0.0365349  },
//                                          { -0.06607774, -0.33526875, -0.15803078 },
//                                          {  0.09738707, -0.49482218, -0.39584868 },
//                                          { -0.45599864,  0.08585038, -0.09432794 },
//                                          {  0.23055539,  0.07564162,  0.24051469 } };
//
////    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
////                                           {  0.00294417,  0.34480654, -0.34450988 },
////                                           {  0.21056511, -0.28113967, -0.45726839 },
////                                           { -0.26431286, -0.09985721,  0.47322301 },
////                                           { -0.18156887, -0.32226199, -0.37295847 },
////                                           {  0.15062371,  0.09439471,  0.09167948 },
////                                           { -0.46869859,  0.018301  ,  0.45013866 },
////                                           { -0.15455446,  0.40552715, -0.4216042  },
////                                           { -0.38930237,  0.10974753, -0.31188239 } };
//
////    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0} };
//
//    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0} };
//                                           
//
//    //Initialize the state variable vector
//    std::vector< double > SDVSDefault( 55, 0 );
//
//    //Initialize the additional degree of freedom vectors
//    std::vector< double > current_ADD_DOF;
//    std::vector< std::vector< double > > current_ADD_grad_DOF;
//
//    std::vector< double > previous_ADD_DOF;
//    std::vector< std::vector< double > > previous_ADD_grad_DOF;
//
//    //Initialize the stress measures
//    std::vector< double > current_PK2( 9, 0 );
//
//    std::vector< double > current_SIGMA( 9, 0 );
//
//    std::vector< double > current_M( 27, 0 );
//
//    //Initialize the additional terms vector
//    std::vector< std::vector< double > > ADD_TERMS;
//
//    //Initialize the output message string
//    std::string output_message;
//
//#ifdef DEBUG_MODE
//    solverTools::homotopyMap DEBUG;
//#endif
//
//    solverTools::floatVector PK2_answer = { 172.484,   15.3785,   -0.917177,
//                                             13.4848, 142.823,    -0.0214307,
//                                             -1.7635,   1.77719, 141.069 };
//
//    solverTools::floatVector SIGMA_answer = { 176.916,   15.8646,   -2.83731,
//                                               15.8646, 144.538,     1.85836,
//                                               -2.83731,  1.85836, 142.013 };
//
//    solverTools::floatVector M_answer = { 0.598283, -0.512218,  0.620664,    3.22636,   1.16682,
//                                          1.20593,   0.562825, -2.52317,     1.62616,  -2.61391,
//                                         -0.60994,  -1.02147,   0.668187,    0.49348,  -0.23916,
//                                         -2.77419,   0.760483,  1.71784,    -0.499389,  2.62828,
//                                         -0.761044,  1.23369,  -0.00778206, -2.25643,  -0.729551,
//                                          0.743204,  0.910521 };
//
//    solverTools::floatVector SDVS_answer = { -1.79592e-24,  0.0243222,    0.0822384,    0.0430345,   0.0435752,
//                                             -8.96006e-25,  0.00852191,   0.0465339,    0.0243507,   0.0246566,
//                                              0.00742998,   0.00500421,  -0.000296486,  0.00498757, -0.00260492,
//                                              0.000284355, -0.000367318,  0.000222511, -0.0015603,   0.00863313,
//                                              0.00537105,  -0.000347686,  0.00643802,  -0.00298667,  0.000297105,
//                                             -0.000422398,  0.000293946, -0.00181986,   0.0385522,  -0.0244021,
//                                             -0.00912035,   0.0105171,   -0.000328615,  0.0222843,   0.0185626,
//                                             -0.0118234,   -0.00785555,   0.0451085,    0.031607,    0.0212748,
//                                              0.0116981,    0.0161821,    0.00126031,   0.0147688,   0.00691805,
//                                             -0.0241431,    0.00942608,  -0.0350366,    0.0221571,  -0.0249697,
//                                              0.00935849,   0.0214931,    0.0169609,    0.0177352,   0.0203838 };
//
//    std::vector< double > SDVS = SDVSDefault;
//
//    std::vector< double > PK2_result, SIGMA_result, M_result;
//
//    //Evaluate the model
//    int errorCode = material->evaluate_model( time, fparams,
//                                              current_grad_u, current_phi, current_grad_phi,
//                                              previous_grad_u, previous_phi, previous_grad_phi,
//                                              SDVS,
//                                              current_ADD_DOF, current_ADD_grad_DOF,
//                                              previous_ADD_DOF, previous_ADD_grad_DOF,
//                                              PK2_result, SIGMA_result, M_result,
//                                              ADD_TERMS,
//                                              output_message
//#ifdef DEBUG_MODE
//                                              , DEBUG
//#endif
//                                            );
//
//    BOOST_CHECK( errorCode == 0 );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVS_answer ) );
//
//    //Check the Jacobian using the previously tested jacobian
//    std::vector< std::vector< double > > DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer,
//                                         DSIGMADgrad_u_answer, DSIGMADphi_answer, DSIGMADgrad_phi_answer,
//                                         DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer;
//
//    std::vector< std::vector< std::vector< double > > > ADD_JACOBIANS;
//
//    SDVS = SDVSDefault;
//
//#ifdef DEBUG_MODE
//    DEBUG.clear();
//#endif
//
//    errorCode = micromorphicElastoPlasticity::evaluate_model(
//                                time, fparams,
//                                current_grad_u, current_phi, current_grad_phi,
//                                previous_grad_u, previous_phi, previous_grad_phi,
//                                SDVS,
//                                current_ADD_DOF, current_ADD_grad_DOF,
//                                previous_ADD_DOF, previous_ADD_grad_DOF,
//                                PK2_result, SIGMA_result, M_result,
//                                DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer,
//                                DSIGMADgrad_u_answer, DSIGMADphi_answer, DSIGMADgrad_phi_answer,
//                                DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer,
//                                ADD_TERMS, ADD_JACOBIANS,
//                                output_message
//#ifdef DEBUG_MODE
//                                , DEBUG
//#endif
//                              );
//
//    BOOST_CHECK( errorCode <= 0 );
//
//    PK2_result.clear();
//    SIGMA_result.clear();
//    M_result.clear();
//
//    SDVS = SDVSDefault;
//
//    std::vector< std::vector< double > > DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
//                                         DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
//                                         DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result;
//
//    errorCode = material->evaluate_model( time, fparams,
//                                          current_grad_u, current_phi, current_grad_phi,
//                                          previous_grad_u, previous_phi, previous_grad_phi,
//                                          SDVS,
//                                          current_ADD_DOF, current_ADD_grad_DOF,
//                                          previous_ADD_DOF, previous_ADD_grad_DOF,
//                                          PK2_result, SIGMA_result, M_result,
//                                          DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
//                                          DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
//                                          DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result,
//                                          ADD_TERMS, ADD_JACOBIANS,
//                                          output_message
//#ifdef DEBUG_MODE
//                                          , DEBUG
//#endif
//                                        );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS_answer, SDVS ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dgrad_u_result, DPK2Dgrad_u_answer ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dphi_result, DPK2Dphi_answer ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADgrad_u_result, DSIGMADgrad_u_answer ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADphi_result, DSIGMADphi_answer ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DMDgrad_u_result, DMDgrad_u_answer ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DMDphi_result, DMDphi_answer ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DMDgrad_phi_result, DMDgrad_phi_answer ) );
//
//#ifdef DEBUG_MODE
//    DEBUG.clear();
//#endif
//
//    //Test the computed numeric Jacobian values
//    SDVS = SDVSDefault;
//    errorCode = material->evaluate_model_numeric_gradients( time, fparams,
//                                                            current_grad_u, current_phi, current_grad_phi,
//                                                            previous_grad_u, previous_phi, previous_grad_phi,
//                                                            SDVS,
//                                                            current_ADD_DOF, current_ADD_grad_DOF,
//                                                            previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                            PK2_result, SIGMA_result, M_result,
//                                                            DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
//                                                            DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
//                                                            DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result,
//                                                            ADD_TERMS, ADD_JACOBIANS,
//                                                            output_message,
//#ifdef DEBUG_MODE
//                                                            DEBUG,
//#endif
//                                                            1e-6 );
//
//    BOOST_CHECK( errorCode <= 0 );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVS_answer ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dgrad_u_result, DPK2Dgrad_u_answer, 1e-4 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dphi_result, DPK2Dphi_answer, 1e-4 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer, 1e-4, 1e-5 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADgrad_u_result, DSIGMADgrad_u_answer, 1e-4 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADphi_result, DSIGMADphi_answer, 1e-4 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer, 1e-4 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DMDgrad_u_result, DMDgrad_u_answer, 1e-4 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DMDphi_result, DMDphi_answer, 1e-4 ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( DMDgrad_phi_result, DMDgrad_phi_answer, 1e-4 ) );
//
//}
//
//BOOST_AUTO_TEST_CASE( testMaterialLibraryInterface2 ){
//    /*!
//     * Test the interface to the linear elastic model
//     * via the material library.
//     *
//     * NOTE: This function mostly exists to perform debugging
//     *       on the implementation of the function into a 
//     *       larger solver code.
//     *
//     */
//
//    //Initialize the model
//    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
//    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
//    auto material = factory.GetMaterial(_model_name);
//
//    //Set up the inputs
//    //Initialize the time
//    std::vector< double > time = { 0.045, 0.01 };
//
//    //Initialize the material parameters
//    std::vector< double > fparams = { 2, 170, 15, 2, 140, 20, 2, 2, 27, 2, 0.56, 0.2, 2, 0.15, 0.3, 2, 0.82, 0.1, 2, 0.42, 0.3, 2, 0.05, 0.2, 2, 0.52, 0.4, 2, 29480, 25480, 5, 1000, 400, -1500, -1400, -3000, 11, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0, 2, 400, -3000, 0.5, 0.5, 0.5, 1e-09, 1e-09 };
//
//    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] =
//    {
//        { -0.00124343, -6.55319e-14, 3.99657e-13},
//        { 0, 0.0045, 0},
//        { -1.75135e-13, -1.35481e-13, -0.00124343 },
//    };
//
//    double previous_grad_u[ 3 ][ 3 ] =
//    {
//        { -0.00123858, -1.22379e-17, 5.04154e-18},
//        { 0, 0.004, 0},
//        { -1.47723e-18, 4.44523e-18, -0.00123858 },
//    };
//
//    //Initialize the micro displacement
//    double current_phi[ 9 ] = { -0.00153489, -3.04626e-13, 5.16537e-13, 1.58771e-13, 0.00303407, 4.29828e-14, -4.38368e-13, -1.80694e-13, -0.00153489 };
//
//    double previous_phi[ 9 ] = { -0.00164749, -2.63663e-17, 1.35603e-17, 8.65138e-19, 0.00325613, -2.13082e-20, -1.17433e-17, 2.24626e-18, -0.00164749 };
//
//    //Initialize the gradient of the micro displacement
//    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0} };
//
//    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0} };
//                                           
//
//    //Initialize the state variable vector
//    std::vector< double > SDVSDefault( 55, 0 );
//
//    //Initialize the additional degree of freedom vectors
//    std::vector< double > current_ADD_DOF;
//    std::vector< std::vector< double > > current_ADD_grad_DOF;
//
//    std::vector< double > previous_ADD_DOF;
//    std::vector< std::vector< double > > previous_ADD_grad_DOF;
//
//    //Initialize the stress measures
//    std::vector< double > current_PK2( 9, 0 );
//
//    std::vector< double > current_SIGMA( 9, 0 );
//
//    std::vector< double > current_M( 27, 0 );
//
//    //Initialize the additional terms vector
//    std::vector< std::vector< double > > ADD_TERMS;
//
//    //Initialize the output message string
//    std::string output_message;
//
//#ifdef DEBUG_MODE
//    solverTools::homotopyMap DEBUG;
//#endif
//
//    std::vector< double > SDVS = SDVSDefault;
//
//    std::vector< double > PK2_result, SIGMA_result, M_result;
//
//    //Evaluate the model
//    int errorCode = material->evaluate_model( time, fparams,
//                                              current_grad_u, current_phi, current_grad_phi,
//                                              previous_grad_u, previous_phi, previous_grad_phi,
//                                              SDVS,
//                                              current_ADD_DOF, current_ADD_grad_DOF,
//                                              previous_ADD_DOF, previous_ADD_grad_DOF,
//                                              PK2_result, SIGMA_result, M_result,
//                                              ADD_TERMS,
//                                              output_message
//#ifdef DEBUG_MODE
//                                              , DEBUG
//#endif
//                                            );
//
//    BOOST_CHECK( errorCode <= 0 );
//
//    std::cout << "SDVS:\n"; vectorTools::print( SDVS );
//    std::cout << "PK2:\n"; vectorTools::print( PK2_result );
//    std::cout << "SIGMA:\n"; vectorTools::print( SIGMA_result );
//    std::cout << "M:\n"; vectorTools::print( M_result );
//
//#ifdef DEBUG_MODE
//    for ( auto step = DEBUG.begin(); step != DEBUG.end(); step++ ){
//        std::cout << step->first << "\n";
//        for ( auto iteration = step->second.begin(); iteration != step->second.end(); iteration++ ){
//            std::cout << "    " << iteration->first << "\n";
//            for ( auto value = iteration->second.begin(); value != iteration->second.end(); value++ ){
//                if ( value->second.size() <= 27 ) {
//                    std::cout << "        " << value->first << "\n";
//                    std::cout << "            "; vectorTools::print( value->second );
//                }
//            }
//        }
//    }
//#endif
//
//    return 0;
//}

BOOST_AUTO_TEST_CASE( testEvaluate_model_continuation){
    /*!
     * Test the evaluation of the constitutive model when being
     * continued from a previous plastic deformation
     *
     */

    //Initialize the time
    std::vector< double > time = { 10., 2.5 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 2.4e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.4e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e0, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double current_grad_u0[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
                                         {0.100, 0.001, 0.000 },
                                         {0.000, 0.000, 0.000 } };

    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
                                         {0, 0, 0},
                                         {0, 0, 0} };
    //Initialize the micro displacement
//    double current_phi[ 9 ] = { 0.84729289,  0.40617104,  0.59534561,  
//                                0.44195587,  0.34121966, -0.79098944, 
//                               -0.43965428,  0.88466225,  0.1684519 };
//
//    double previous_phi[ 9 ] = { -0.99935855, -0.21425717,  0.0668254 ,
//                                 -0.11111872, -0.07416114, -1.01048108,
//                                  0.1804018 , -1.01116291,  0.03248007 };

    double current_phi0[ 9 ] = { 0.100, 0.000, 0.000,
                                 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000 };

    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi0[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
                                           { -0.18572739,  0.06847269,  0.22931628 },
                                           { -0.01829735, -0.48731265, -0.25277529 },
                                           {  0.26626212,  0.4844646 , -0.31965177 },
                                           {  0.49197846,  0.19051656, -0.0365349  },
                                           { -0.06607774, -0.33526875, -0.15803078 },
                                           {  0.09738707, -0.49482218, -0.39584868 },
                                           { -0.45599864,  0.08585038, -0.09432794 },
                                           {  0.23055539,  0.07564162,  0.24051469 } };

//    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
//                                           {  0.00294417,  0.34480654, -0.34450988 },
//                                           {  0.21056511, -0.28113967, -0.45726839 },
//                                           { -0.26431286, -0.09985721,  0.47322301 },
//                                           { -0.18156887, -0.32226199, -0.37295847 },
//                                           {  0.15062371,  0.09439471,  0.09167948 },
//                                           { -0.46869859,  0.018301  ,  0.45013866 },
//                                           { -0.15455446,  0.40552715, -0.4216042  },
//                                           { -0.38930237,  0.10974753, -0.31188239 } };

//    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0} };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };
                                           

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > current_PK20( 9, 0 );

    std::vector< double > current_SIGMA0( 9, 0 );

    std::vector< double > current_M0( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

#ifdef DEBUG_MODE
    solverTools::homotopyMap homotopyDEBUG0;
#endif
    solverTools::floatVector PK2Answer0 = { 172.484,  15.3785,    -0.917177,
                                             13.4848, 142.823,    -0.0214307,
                                             -1.7635,   1.77719, 141.069 };

    solverTools::floatVector SigmaAnswer0 = { 176.916,   15.8646,   -2.83731,
                                               15.8646, 144.538,     1.85836,
                                               -2.83731,  1.85836, 142.013 };

    solverTools::floatVector MAnswer0 = { 0.598283, -0.512218,  0.620664,    3.22636,  1.16682,
                                          1.20593,   0.562825, -2.52317,     1.62616, -2.61391,
                                         -0.60994,  -1.02147,   0.668187,    0.49348, -0.23916,
                                         -2.77419,   0.760483,  1.71784,    -0.499389, 2.62828,
                                         -0.761044,  1.23369,  -0.00778206, -2.25643, -0.729551,
                                          0.743204,  0.910521 };

    solverTools::floatVector SDVSAnswer0 = { -1.79592e-24,  0.0243222,    0.0822384,    0.0430345,   0.0435752,
                                             -8.96006e-25,  0.00852191,   0.0465339,    0.0243507,   0.0246566,
                                              0.00742998,   0.00500421,  -0.000296486,  0.00498757, -0.00260492,
                                              0.000284355, -0.000367318,  0.000222511, -0.0015603,   0.00863313,
                                              0.00537105,  -0.000347686,  0.00643802,  -0.00298667,  0.000297105,
                                             -0.000422398,  0.000293946, -0.00181986,   0.0385522,  -0.0244021,
                                             -0.00912035,   0.0105171,   -0.000328615,  0.0222843,   0.0185626,
                                             -0.0118234,   -0.00785555,   0.0451085,    0.031607,    0.0212748,
                                              0.0116981,    0.0161821,    0.00126031,   0.0147688,   0.00691805,
                                             -0.0241431,    0.00942608,  -0.0350366,    0.0221571,  -0.0249697,
                                              0.00935849,   0.0214931,    0.0169609,    0.0177352,   0.0203838 };

    std::vector< double > SDVS = SDVSDefault;

    int errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                                  current_grad_u0,  current_phi0,  current_grad_phi0,
                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                  SDVS,
                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                  current_PK20, current_SIGMA0, current_M0,
                                                                  ADD_TERMS,
                                                                  output_message
#ifdef DEBUG_MODE
                                                                  , homotopyDEBUG0
#endif
                                                                  );

    BOOST_CHECK( errorCode == 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVSAnswer0 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2Answer0, current_PK20, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SigmaAnswer0, current_SIGMA0, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( MAnswer0, current_M0, 1e-5 ) );

#ifdef DEBUG_MODE
    solverTools::homotopyMap homotopyDEBUG1;
#endif

    double current_grad_u1[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u1[ i ][ j ] = current_grad_u0[ i ][ j ];
        }
    }

    double current_phi1[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi1[ i ] = current_phi0[ i ];
    }

    double current_grad_phi1[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi1[ i ][ j ] = current_grad_phi0[ i ][ j ];
        }
    }

    //Compute the new stress
    std::vector< double > current_PK21( 9, 0 );

    std::vector< double > current_SIGMA1( 9, 0 );

    std::vector< double > current_M1( 27, 0 );

    std::vector< double  > SDVSAnswer1 = { -2.9932e-24,   0.0347461,    0.126521,     0.0662069,   0.0670387,
                                            0,            0,            0,            0,           0,
                                            0.0123833,    0.00834036,  -0.000494143,  0.00831262, -0.00434154,
                                            0.000473925, -0.000612196,  0.000370852, -0.0026005,   0.012333,
                                            0.00767293,  -0.000496694,  0.00919718,  -0.00426667,  0.000424436,
                                           -0.000603425,  0.000419922, -0.0025998,    0.0593111,  -0.0375416,
                                           -0.0140313,    0.0161802,   -0.000505561,  0.0342836,   0.0285578,
                                           -0.0181899,   -0.0120855,    0.0693977,    0.0486262,   0.0327304,
                                            0.0179971,    0.0248956,    0.00193894,   0.0227212,   0.0106431,
                                           -0.0371433,    0.0145017,   -0.0539024,    0.0340878,  -0.038415,
                                            0.0143977,    0.0330663,    0.0260937,    0.0272849,   0.0313596 };

    std::vector< double > PK2Answer1 = { 1.66486805e+02,  1.34771535e+01, -6.65420292e-01,
                                         1.19267084e+01,  1.40317255e+02,  1.45612760e-03,
                                        -1.42979700e+00,  1.58068588e+00,  1.38540666e+02 };

    std::vector< double > SIGMAAnswer1 = { 170.72445498,  14.01774493,  -2.35073512,
                                            14.01774493, 142.08713184,   1.75286437,
                                            -2.35073512,   1.75286437, 139.65940169 };

    std::vector< double > MAnswer1 = { 0.41224802, -0.49821007,  0.57679082,  2.85720915,  1.13385061,
                                       1.09114466,  0.49605673, -2.3437262 ,  1.54092504, -2.39874562,
                                      -0.70509639, -0.98570884,  0.6169238 ,  0.36560779, -0.26418717,
                                      -2.5593553 ,  0.73214846,  1.54474417, -0.47291518,  2.49428362,
                                      -0.79126967,  1.05309185, -0.03402865, -2.10298011, -0.7864686 ,
                                       0.63759613,  0.82637607 };

    errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                              current_grad_u1, current_phi1, current_grad_phi1,
                                                              current_grad_u0, current_phi0, current_grad_phi0,
                                                              SDVS,
                                                              current_ADD_DOF,  current_ADD_grad_DOF,
                                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                                              current_PK21, current_SIGMA1, current_M1,
                                                              ADD_TERMS,
                                                              output_message
#ifdef DEBUG_MODE
                                                              , homotopyDEBUG1
#endif
                                                              );
    std::vector< double > SDVS1 = SDVS;

    //Make sure that the state variables are monotonically increasing
    solverTools::floatVector ISVsDelta = solverTools::floatVector( SDVS.begin(), SDVS.begin() + 5 )
                                       - solverTools::floatVector( SDVSAnswer0.begin(), SDVSAnswer0.begin() + 5 );

    for ( unsigned int i = 0; i < ISVsDelta.size(); i++ ){
        BOOST_CHECK( ISVsDelta[ i ] >= -1e-9 );
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVSAnswer1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_PK21, PK2Answer1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_SIGMA1, SIGMAAnswer1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_M1, MAnswer1 ) );

#ifdef DEBUG_MODE

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPK2Stress" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "intermediatePK2Stress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousReferenceMicroStress" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "intermediateReferenceMicroStress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousReferenceHigherOrderStress" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "intermediateReferenceHigherOrderStress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPlasticDeformationGradient" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "convergedPlasticDeformationGradient" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPlasticMicroDeformation" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "convergedPlasticMicroDeformation" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPlasticGradientMicroDeformation" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "convergedPlasticGradientMicroDeformation" ] ) );

#endif
    
    //Compute the new stress
    std::vector< double > current_PK22( 9, 0 );

    std::vector< double > current_SIGMA2( 9, 0 );

    std::vector< double > current_M2( 27, 0 );

#ifdef DEBUG_MODE
    solverTools::homotopyMap homotopyDEBUG2;
#endif

    errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                              current_grad_u1, current_phi1, current_grad_phi1,
                                                              current_grad_u1, current_phi1, current_grad_phi1,
                                                              SDVS,
                                                              current_ADD_DOF,  current_ADD_grad_DOF,
                                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                                              current_PK22, current_SIGMA2, current_M2,
                                                              ADD_TERMS,
                                                              output_message
#ifdef DEBUG_MODE
                                                              , homotopyDEBUG2
#endif
                                                              );

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVS1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_PK21, current_PK22 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_SIGMA1, current_SIGMA2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_M1, current_M2 ) );

}

BOOST_AUTO_TEST_CASE( testEvaluate_model_continuation2){
    /*!
     * Test the evaluation of the constitutive model when being
     * continued from a previous plastic deformation
     *
     */

    //Initialize the time
    std::vector< double > time = { 10., 2.5 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 2.4e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.4e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e0, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double current_grad_u0[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
                                         {0.100, 0.001, 0.000 },
                                         {0.000, 0.000, 0.000 } };

    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
                                         {0, 0, 0},
                                         {0, 0, 0} };
    //Initialize the micro displacement
//    double current_phi[ 9 ] = { 0.84729289,  0.40617104,  0.59534561,  
//                                0.44195587,  0.34121966, -0.79098944, 
//                               -0.43965428,  0.88466225,  0.1684519 };
//
//    double previous_phi[ 9 ] = { -0.99935855, -0.21425717,  0.0668254 ,
//                                 -0.11111872, -0.07416114, -1.01048108,
//                                  0.1804018 , -1.01116291,  0.03248007 };

    double current_phi0[ 9 ] = { 0.100, 0.000, 0.000,
                                 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000 };

    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi0[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
                                           { -0.18572739,  0.06847269,  0.22931628 },
                                           { -0.01829735, -0.48731265, -0.25277529 },
                                           {  0.26626212,  0.4844646 , -0.31965177 },
                                           {  0.49197846,  0.19051656, -0.0365349  },
                                           { -0.06607774, -0.33526875, -0.15803078 },
                                           {  0.09738707, -0.49482218, -0.39584868 },
                                           { -0.45599864,  0.08585038, -0.09432794 },
                                           {  0.23055539,  0.07564162,  0.24051469 } };

//    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
//                                           {  0.00294417,  0.34480654, -0.34450988 },
//                                           {  0.21056511, -0.28113967, -0.45726839 },
//                                           { -0.26431286, -0.09985721,  0.47322301 },
//                                           { -0.18156887, -0.32226199, -0.37295847 },
//                                           {  0.15062371,  0.09439471,  0.09167948 },
//                                           { -0.46869859,  0.018301  ,  0.45013866 },
//                                           { -0.15455446,  0.40552715, -0.4216042  },
//                                           { -0.38930237,  0.10974753, -0.31188239 } };

//    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0} };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };
                                           

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > current_PK20( 9, 0 );

    std::vector< double > current_SIGMA0( 9, 0 );

    std::vector< double > current_M0( 27, 0 );

    std::vector< std::vector< double > > DPK2Dgrad_u,   DPK2Dphi,   DPK2Dgrad_phi,
                                         DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                         DMDgrad_u,     DMDphi,     DMDgrad_phi;

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;
    std::vector< std::vector< std::vector< double > > > ADD_JACOBIANS;

    //Initialize the output message string
    std::string output_message;

#ifdef DEBUG_MODE
    solverTools::homotopyMap homotopyDEBUG0;
#endif

    solverTools::floatVector PK2Answer0 = { 172.484,  15.3785,    -0.917177,
                                             13.4848, 142.823,    -0.0214307,
                                             -1.7635,   1.77719, 141.069 };

    solverTools::floatVector SigmaAnswer0 = { 176.916,   15.8646,   -2.83731,
                                               15.8646, 144.538,     1.85836,
                                               -2.83731,  1.85836, 142.013 };

    solverTools::floatVector MAnswer0 = { 0.598283, -0.512218,  0.620664,    3.22636,  1.16682,
                                          1.20593,   0.562825, -2.52317,     1.62616, -2.61391,
                                         -0.60994,  -1.02147,   0.668187,    0.49348, -0.23916,
                                         -2.77419,   0.760483,  1.71784,    -0.499389, 2.62828,
                                         -0.761044,  1.23369,  -0.00778206, -2.25643, -0.729551,
                                          0.743204,  0.910521 };

    solverTools::floatVector SDVSAnswer0 = { -1.79592e-24,  0.0243222,    0.0822384,    0.0430345,   0.0435752,
                                             -8.96006e-25,  0.00852191,   0.0465339,    0.0243507,   0.0246566,
                                              0.00742998,   0.00500421,  -0.000296486,  0.00498757, -0.00260492,
                                              0.000284355, -0.000367318,  0.000222511, -0.0015603,   0.00863313,
                                              0.00537105,  -0.000347686,  0.00643802,  -0.00298667,  0.000297105,
                                             -0.000422398,  0.000293946, -0.00181986,   0.0385522,  -0.0244021,
                                             -0.00912035,   0.0105171,   -0.000328615,  0.0222843,   0.0185626,
                                             -0.0118234,   -0.00785555,   0.0451085,    0.031607,    0.0212748,
                                              0.0116981,    0.0161821,    0.00126031,   0.0147688,   0.00691805,
                                             -0.0241431,    0.00942608,  -0.0350366,    0.0221571,  -0.0249697,
                                              0.00935849,   0.0214931,    0.0169609,    0.0177352,   0.0203838 };

    std::vector< double > SDVS = SDVSDefault;

    int errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                                  current_grad_u0,  current_phi0,  current_grad_phi0,
                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                  SDVS,
                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                  current_PK20, current_SIGMA0, current_M0,
                                                                  DPK2Dgrad_u,   DPK2Dphi,   DPK2Dgrad_phi,
                                                                  DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                                                  DMDgrad_u,     DMDphi,     DMDgrad_phi,
                                                                  ADD_TERMS, ADD_JACOBIANS,
                                                                  output_message
#ifdef DEBUG_MODE
                                                                  , homotopyDEBUG0
#endif
                                                                  );

    BOOST_CHECK( errorCode == 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVSAnswer0 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2Answer0, current_PK20, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SigmaAnswer0, current_SIGMA0, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( MAnswer0, current_M0, 1e-5 ) );

#ifdef DEBUG_MODE
    solverTools::homotopyMap homotopyDEBUG1;
#endif

    double current_grad_u1[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u1[ i ][ j ] = current_grad_u0[ i ][ j ];
        }
    }

    double current_phi1[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi1[ i ] = current_phi0[ i ];
    }

    double current_grad_phi1[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi1[ i ][ j ] = current_grad_phi0[ i ][ j ];
        }
    }

    //Compute the new stress
    std::vector< double > current_PK21( 9, 0 );

    std::vector< double > current_SIGMA1( 9, 0 );

    std::vector< double > current_M1( 27, 0 );

    std::vector< double  > SDVSAnswer1 = { -2.9932e-24,   0.0347461,    0.126521,     0.0662069,   0.0670387,
                                            0,            0,            0,            0,           0,
                                            0.0123833,    0.00834036,  -0.000494143,  0.00831262, -0.00434154,
                                            0.000473925, -0.000612196,  0.000370852, -0.0026005,   0.012333,
                                            0.00767293,  -0.000496694,  0.00919718,  -0.00426667,  0.000424436,
                                           -0.000603425,  0.000419922, -0.0025998,    0.0593111,  -0.0375416,
                                           -0.0140313,    0.0161802,   -0.000505561,  0.0342836,   0.0285578,
                                           -0.0181899,   -0.0120855,    0.0693977,    0.0486262,   0.0327304,
                                            0.0179971,    0.0248956,    0.00193894,   0.0227212,   0.0106431,
                                           -0.0371433,    0.0145017,   -0.0539024,    0.0340878,  -0.038415,
                                            0.0143977,    0.0330663,    0.0260937,    0.0272849,   0.0313596 };

    std::vector< double > PK2Answer1 = { 1.66486805e+02,  1.34771535e+01, -6.65420292e-01,
                                         1.19267084e+01,  1.40317255e+02,  1.45612760e-03,
                                        -1.42979700e+00,  1.58068588e+00,  1.38540666e+02 };

    std::vector< double > SIGMAAnswer1 = { 170.72445498,  14.01774493,  -2.35073512,
                                            14.01774493, 142.08713184,   1.75286437,
                                            -2.35073512,   1.75286437, 139.65940169 };

    std::vector< double > MAnswer1 = { 0.41224802, -0.49821007,  0.57679082,  2.85720915,  1.13385061,
                                       1.09114466,  0.49605673, -2.3437262 ,  1.54092504, -2.39874562,
                                      -0.70509639, -0.98570884,  0.6169238 ,  0.36560779, -0.26418717,
                                      -2.5593553 ,  0.73214846,  1.54474417, -0.47291518,  2.49428362,
                                      -0.79126967,  1.05309185, -0.03402865, -2.10298011, -0.7864686 ,
                                       0.63759613,  0.82637607 };

    std::vector< std::vector< double > > DPK2Dgrad_u1,   DPK2Dphi1,   DPK2Dgrad_phi1,
                                         DSIGMADgrad_u1, DSIGMADphi1, DSIGMADgrad_phi1,
                                         DMDgrad_u1,     DMDphi1,     DMDgrad_phi1;

    errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                              current_grad_u1, current_phi1, current_grad_phi1,
                                                              current_grad_u0, current_phi0, current_grad_phi0,
                                                              SDVS,
                                                              current_ADD_DOF,  current_ADD_grad_DOF,
                                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                                              current_PK21, current_SIGMA1, current_M1,
                                                              DPK2Dgrad_u1,   DPK2Dphi1,   DPK2Dgrad_phi1,
                                                              DSIGMADgrad_u1, DSIGMADphi1, DSIGMADgrad_phi1,
                                                              DMDgrad_u1,     DMDphi1,     DMDgrad_phi1,
                                                              ADD_TERMS, ADD_JACOBIANS,
                                                              output_message
#ifdef DEBUG_MODE
                                                              , homotopyDEBUG1
#endif
                                                              );
    std::vector< double > SDVS1 = SDVS;

    //Make sure that the state variables are monotonically increasing
    solverTools::floatVector ISVsDelta = solverTools::floatVector( SDVS.begin(), SDVS.begin() + 5 )
                                       - solverTools::floatVector( SDVSAnswer0.begin(), SDVSAnswer0.begin() + 5 );

    for ( unsigned int i = 0; i < ISVsDelta.size(); i++ ){
        BOOST_CHECK( ISVsDelta[ i ] >= -1e-9 );
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVSAnswer1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_PK21, PK2Answer1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_SIGMA1, SIGMAAnswer1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_M1, MAnswer1 ) );

#ifdef DEBUG_MODE

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPK2Stress" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "intermediatePK2Stress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousReferenceMicroStress" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "intermediateReferenceMicroStress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousReferenceHigherOrderStress" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "intermediateReferenceHigherOrderStress" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPlasticDeformationGradient" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "convergedPlasticDeformationGradient" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPlasticMicroDeformation" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "convergedPlasticMicroDeformation" ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPlasticGradientMicroDeformation" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "convergedPlasticGradientMicroDeformation" ] ) );

#endif
    
    //Compute the new stress
    std::vector< double > current_PK22( 9, 0 );

    std::vector< double > current_SIGMA2( 9, 0 );

    std::vector< double > current_M2( 27, 0 );

    std::vector< std::vector< double > > DPK2Dgrad_u2,   DPK2Dphi2,   DPK2Dgrad_phi2,
                                         DSIGMADgrad_u2, DSIGMADphi2, DSIGMADgrad_phi2,
                                         DMDgrad_u2,     DMDphi2,     DMDgrad_phi2;

#ifdef DEBUG_MODE
    solverTools::homotopyMap homotopyDEBUG2;
#endif

    errorCode = micromorphicElastoPlasticity::evaluate_model( time, fparams,
                                                              current_grad_u1, current_phi1, current_grad_phi1,
                                                              current_grad_u1, current_phi1, current_grad_phi1,
                                                              SDVS,
                                                              current_ADD_DOF,  current_ADD_grad_DOF,
                                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                                              current_PK22, current_SIGMA2, current_M2,
                                                              DPK2Dgrad_u2,   DPK2Dphi2,   DPK2Dgrad_phi2,
                                                              DSIGMADgrad_u2, DSIGMADphi2, DSIGMADgrad_phi2,
                                                              DMDgrad_u2,     DMDphi2,     DMDgrad_phi2,
                                                              ADD_TERMS, ADD_JACOBIANS,
                                                              output_message
#ifdef DEBUG_MODE
                                                              , homotopyDEBUG2
#endif
                                                              );

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVS1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_PK21, current_PK22 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_SIGMA1, current_SIGMA2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( current_M1, current_M2 ) );

}

//BOOST_AUTO_TEST_CASE( testEvaluate_model_history ){
//    /*!
//     * Test the material model undergoing a time history.
//     *
//     */
//
//#ifdef DEBUG_MODE
//    std::ofstream output_file;
//    output_file.open( "output_file.txt" );
//#endif
//
//    //Initialize the model
//    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
//    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
//    auto material = factory.GetMaterial(_model_name);
//
//    std::vector< std::vector< double > > grad_u_0 = { { 0, 0, 0 },
//                                                      { 0, 0, 0 },
//                                                      { 0, 0, 0 } };
//
//    std::vector< std::vector< double > > grad_u_f = { { 0.5, 0, 0 },
//                                                      { 0.0, 0, 0 },
//                                                      { 0.0, 0, 0 } };
//
//    std::vector< double > phi_0 = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//    std::vector< double > phi_f = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//
//    std::vector< std::vector< double > > grad_phi_0 = { { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0} };
//
//    std::vector< std::vector< double > > grad_phi_f = { { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0},
//                                                        { 0, 0, 0} };
//
//    double dt = 0.05;
//    double t0 = 0.;
//    double tf = 0.25;
//
//    double t = t0;
//
//    //Set up the model parameters
//    std::vector< double > fparams = { 2, 1e3, 1e2,
//                                      2, 7e2, 1e4,
//                                      2, 1e3, 1e4,
//                                      2, 0., 0.0,
//                                      2, 0., 0.0,
//                                      2, 0., 0.0,
//                                      2, 0., 0.0,
//                                      2, 0., 0.0,
//                                      2, 0., 0.0,
//                                      2, 29480, 25480,
//                                      5, 1000, 400, -1500, -1400, -3000,
//                                      11, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0,
//                                      2, 400, -3000,
//                                      0.5, 0.5, 0.5, 1e-09, 1e-09 };
////                                      0.0, 0.0, 0.0, 1e-09, 1e-09 };
//
//    //Initialize the state variable vector
//    std::vector< double > SDVSDefault( 55, 0 );
//
//    //Initialize the additional degree of freedom vectors
//    std::vector< double > current_ADD_DOF;
//    std::vector< std::vector< double > > current_ADD_grad_DOF;
//
//    std::vector< double > previous_ADD_DOF;
//    std::vector< std::vector< double > > previous_ADD_grad_DOF;
//
//    //Initialize the stress measures
//    std::vector< double > current_PK2( 9, 0 );
//
//    std::vector< double > current_SIGMA( 9, 0 );
//
//    std::vector< double > current_M( 27, 0 );
//
//    //Initialize the additional terms vector
//    std::vector< std::vector< double > > ADD_TERMS;
//
//    //Initialize the output message string
//    std::string output_message;
//
//#ifdef DEBUG_MODE
//    solverTools::homotopyMap DEBUG;
//#endif
//
//    std::vector< double > SDVS = SDVSDefault;
//
//    std::vector< double > PK2_result, SIGMA_result, M_result;
//
//    std::vector< double > PK2Answer = {  5.16027494e+03, -4.30718970e-18,  4.23714575e-18,
//                                        -4.30349349e-18,  4.05236467e+03, -1.78599775e-18,
//                                         4.25600330e-18, -1.80293236e-18,  4.05236467e+03 };
//
//    std::vector< double > SIGMAAnswer = {  5051.82,     -4.38296e-18,  4.03864e-18,
//                                          -4.38296e-18,  4116.16,     -1.6831e-18,
//                                           4.03864e-18, -1.6831e-18,   4116.16 };
//
//    std::vector< double > MAnswer = {  1.0019e-46,  -9.52574e-49,  1.00995e-48, -8.8467e-25,   8.11869e-46,
//                                      -5.63453e-34, -1.10865e-32,  2.44433e-34,  1.6875e-49,  -9.38296e-38,
//                                       4.14595e-47, -6.23442e-59, -1.42808e-16,  1.31056e-37, -9.48874e-38,
//                                      -8.77566e-16,  8.05349e-37, -5.83091e-37, -1.10003e-37, -1.80712e-47,
//                                       3.77743e-56, -1.06975e-37,  9.81714e-59, -1.63467e-56, -8.77566e-16,
//                                       8.05349e-37, -1.05373e-36 };
//
//    std::vector< double > SDVSAnswer = { 0.0536436,    0.0316833,    0,            0,            0,
//                                         0.441627,     2.79534e-10,  0,            0,            0,
//                                         0.0406588,   -4.90538e-23,  6.74381e-23,  2.31867e-22, -0.0197043,
//                                        -2.04791e-23, -2.35533e-22,  8.57369e-23, -0.0197043,    0.0150881,
//                                        -5.35887e-23, -5.0789e-24,   2.923e-22,   -0.00745964,  -6.81252e-24,
//                                        -8.98761e-24,  1.5285e-24,  -0.00745964,   1.47911e-31, -4.70981e-22,
//                                        -1.60422e-39, -4.82898e-51,  1.54956e-54,  5.15605e-62,  2.70794e-38,
//                                         1.60531e-36,  1.84864e-40,  1.93673e-43, -1.471e-31,    2.93537e-22,
//                                         1.85391e-40,  4.83968e-64,  1.1391e-40,   5.43161e-38, -1.48104e-62,
//                                        -1.08569e-49, -2.16264e-24, -7.51219e-32,  1.97215e-31, -5.93004e-39,
//                                        -9.83265e-57,  2.94509e-61,  1.67967e-41, -1.14207e-61, -7.96338e-64 };
//
//    std::vector< std::vector< double > > grad_u_prev   = grad_u_0;
//    std::vector< double > phi_prev                     = phi_0;
//    std::vector< std::vector< double > > grad_phi_prev = grad_phi_0;
//
//    std::vector< std::vector< double > > grad_u_curr;
//    std::vector< double > phi_curr;
//    std::vector< std::vector< double > > grad_phi_curr;
//
//    std::vector< double > time;
//
//    double current_grad_u[ 3 ][ 3 ], current_phi[ 9 ], current_grad_phi[ 9 ][ 3 ];
//    double previous_grad_u[ 3 ][ 3 ], previous_phi[ 9 ], previous_grad_phi[ 9 ][ 3 ];
//
//    //Initial state
//    //Update the arrays
//    for ( unsigned int i = 0; i < 3; i++ ){
//        for ( unsigned int j = 0; j < 3; j++ ){
//            previous_grad_u[ i ][ j ] = grad_u_prev[ i ][ j ];
//        }
//    }
//
//    for ( unsigned int i = 0; i < 9; i++ ){
//        previous_phi[ i ] = phi_prev[ i ];
//
//        for ( unsigned int j = 0; j < 3; j++ ){
//            previous_grad_phi[ i ][ j ] = grad_phi_prev[ i ][ j ];
//        }
//    }
//    //Evaluate the model
//    time = { 0., 0. };
//    int errorCode = material->evaluate_model( time, fparams,
//                                              previous_grad_u, previous_phi, previous_grad_phi,
//                                              previous_grad_u, previous_phi, previous_grad_phi,
//                                              SDVS,
//                                              current_ADD_DOF, current_ADD_grad_DOF,
//                                              previous_ADD_DOF, previous_ADD_grad_DOF,
//                                              PK2_result, SIGMA_result, M_result,
//                                              ADD_TERMS,
//                                              output_message
//#ifdef DEBUG_MODE
//                                              , DEBUG
//#endif
//                                            );
//   
//    BOOST_CHECK( errorCode <= 0);
//#ifdef DEBUG_MODE
//    if( errorCode > 0 ){ 
//        output_file.close();
//    }
//#endif
//
//#ifdef DEBUG_MODE
//    output_file << "NEW_INCREMENT\n";
//
//    //Output the time
//    output_file << time[ 0 ] << ", " << time[ 1 ] << "\n";
//
//    //Output the current gradient of u
//    output_file << current_grad_u[ 0 ][ 0 ] << ", " << current_grad_u[ 0 ][ 1 ] << ", " << current_grad_u[ 0 ][ 2 ] << ", "
//                << current_grad_u[ 1 ][ 0 ] << ", " << current_grad_u[ 1 ][ 1 ] << ", " << current_grad_u[ 1 ][ 2 ] << ", "
//                << current_grad_u[ 2 ][ 0 ] << ", " << current_grad_u[ 2 ][ 1 ] << ", " << current_grad_u[ 2 ][ 2 ] << "\n";
//
//    //Output the current micro displacement
//    output_file << current_phi[ 0 ] << ", " << current_phi[ 1 ] << ", " << current_phi[ 2 ] << ", "
//                << current_phi[ 3 ] << ", " << current_phi[ 4 ] << ", " << current_phi[ 5 ] << ", "
//                << current_phi[ 6 ] << ", " << current_phi[ 7 ] << ", " << current_phi[ 8 ] << "\n";
//
//    //Output the current gradient of the micro displacement
//    output_file << current_grad_phi[ 0 ][ 0 ] << ", " << current_grad_phi[ 0 ][ 1 ] << ", " << current_grad_phi[ 0 ][ 2 ] << ", "
//                << current_grad_phi[ 1 ][ 0 ] << ", " << current_grad_phi[ 1 ][ 1 ] << ", " << current_grad_phi[ 1 ][ 2 ] << ", "
//                << current_grad_phi[ 2 ][ 0 ] << ", " << current_grad_phi[ 2 ][ 1 ] << ", " << current_grad_phi[ 2 ][ 2 ] << ", "
//                << current_grad_phi[ 3 ][ 0 ] << ", " << current_grad_phi[ 3 ][ 1 ] << ", " << current_grad_phi[ 3 ][ 2 ] << ", "
//                << current_grad_phi[ 4 ][ 0 ] << ", " << current_grad_phi[ 4 ][ 1 ] << ", " << current_grad_phi[ 4 ][ 2 ] << ", "
//                << current_grad_phi[ 5 ][ 0 ] << ", " << current_grad_phi[ 5 ][ 1 ] << ", " << current_grad_phi[ 5 ][ 2 ] << ", "
//                << current_grad_phi[ 6 ][ 0 ] << ", " << current_grad_phi[ 6 ][ 1 ] << ", " << current_grad_phi[ 6 ][ 2 ] << ", "
//                << current_grad_phi[ 7 ][ 0 ] << ", " << current_grad_phi[ 7 ][ 1 ] << ", " << current_grad_phi[ 7 ][ 2 ] << ", "
//                << current_grad_phi[ 8 ][ 0 ] << ", " << current_grad_phi[ 8 ][ 1 ] << ", " << current_grad_phi[ 8 ][ 2 ] << "\n";
//
//    //Output the PK2 stress
//    output_file << PK2_result[ 0 ] << ", " << PK2_result[ 1 ] << ", " << PK2_result[ 2 ] << ", "
//                << PK2_result[ 3 ] << ", " << PK2_result[ 4 ] << ", " << PK2_result[ 5 ] << ", "
//                << PK2_result[ 6 ] << ", " << PK2_result[ 7 ] << ", " << PK2_result[ 8 ] << "\n";
//
//    //Output the SIGMA stress
//    output_file << SIGMA_result[ 0 ] << ", " << SIGMA_result[ 1 ] << ", " << SIGMA_result[ 2 ] << ", "
//                << SIGMA_result[ 3 ] << ", " << SIGMA_result[ 4 ] << ", " << SIGMA_result[ 5 ] << ", "
//                << SIGMA_result[ 6 ] << ", " << SIGMA_result[ 7 ] << ", " << SIGMA_result[ 8 ] << "\n";
//
//    //Output the M stress
//    output_file << M_result[  0 ] << ", " << M_result[  1 ] << ", " << M_result[  2 ] << ", "
//                << M_result[  3 ] << ", " << M_result[  4 ] << ", " << M_result[  5 ] << ", "
//                << M_result[  6 ] << ", " << M_result[  7 ] << ", " << M_result[  8 ] << ", "
//                << M_result[  9 ] << ", " << M_result[ 10 ] << ", " << M_result[ 11 ] << ", "
//                << M_result[ 12 ] << ", " << M_result[ 13 ] << ", " << M_result[ 14 ] << ", "
//                << M_result[ 15 ] << ", " << M_result[ 16 ] << ", " << M_result[ 17 ] << ", "
//                << M_result[ 18 ] << ", " << M_result[ 19 ] << ", " << M_result[ 20 ] << ", "
//                << M_result[ 21 ] << ", " << M_result[ 22 ] << ", " << M_result[ 23 ] << ", "
//                << M_result[ 24 ] << ", " << M_result[ 25 ] << ", " << M_result[ 26 ] << "\n";
//
//    std::vector< double > PK2_intermediate, SIGMA_intermediate, M_intermediate;
//
//    //Determine if there were non-linear iterations
//    auto inc = DEBUG.find( "converged_values" );
//    if ( inc != DEBUG.end() ){
//        PK2_intermediate   = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediatePK2Stress" ];
//        SIGMA_intermediate = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceMicroStress" ];
//        M_intermediate     = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceHigherOrderStress" ];
//    }
//    else{
//        PK2_intermediate   = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediatePK2Stress" ];
//        SIGMA_intermediate = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceMicroStress" ];
//        M_intermediate     = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceHigherOrderStress" ];
//    }
//
//    //Output the intermediate PK2 stress
//    output_file << PK2_intermediate[ 0 ] << ", " << PK2_intermediate[ 1 ] << ", " << PK2_intermediate[ 2 ] << ", "
//                << PK2_intermediate[ 3 ] << ", " << PK2_intermediate[ 4 ] << ", " << PK2_intermediate[ 5 ] << ", "
//                << PK2_intermediate[ 6 ] << ", " << PK2_intermediate[ 7 ] << ", " << PK2_intermediate[ 8 ] << "\n";
//
//    //Output the intermediate SIGMA stress
//    output_file << SIGMA_intermediate[ 0 ] << ", " << SIGMA_intermediate[ 1 ] << ", " << SIGMA_intermediate[ 2 ] << ", "
//                << SIGMA_intermediate[ 3 ] << ", " << SIGMA_intermediate[ 4 ] << ", " << SIGMA_intermediate[ 5 ] << ", "
//                << SIGMA_intermediate[ 6 ] << ", " << SIGMA_intermediate[ 7 ] << ", " << SIGMA_intermediate[ 8 ] << "\n";
//
//    //Output the intermediate M stress
//    output_file << M_intermediate[  0 ] << ", " << M_intermediate[  1 ] << ", " << M_intermediate[  2 ] << ", "
//                << M_intermediate[  3 ] << ", " << M_intermediate[  4 ] << ", " << M_intermediate[  5 ] << ", "
//                << M_intermediate[  6 ] << ", " << M_intermediate[  7 ] << ", " << M_intermediate[  8 ] << ", "
//                << M_intermediate[  9 ] << ", " << M_intermediate[ 10 ] << ", " << M_intermediate[ 11 ] << ", "
//                << M_intermediate[ 12 ] << ", " << M_intermediate[ 13 ] << ", " << M_intermediate[ 14 ] << ", "
//                << M_intermediate[ 15 ] << ", " << M_intermediate[ 16 ] << ", " << M_intermediate[ 17 ] << ", "
//                << M_intermediate[ 18 ] << ", " << M_intermediate[ 19 ] << ", " << M_intermediate[ 20 ] << ", "
//                << M_intermediate[ 21 ] << ", " << M_intermediate[ 22 ] << ", " << M_intermediate[ 23 ] << ", "
//                << M_intermediate[ 24 ] << ", " << M_intermediate[ 25 ] << ", " << M_intermediate[ 26 ] << "\n";
//
//    //Output the state variables
//    for ( unsigned int i = 0; i < SDVS.size()-1; i++ ){
//        output_file << SDVS[ i ] << ", ";
//    }
//    output_file << SDVS[ SDVS.size() - 1 ] << "\n";
//
//#endif
//
//
//    //Begin iteration
//    while ( t + dt < tf ){
//
//        time = { t + dt, dt };
//
//        //Increment the displacements
//        grad_u_curr   = grad_u_prev   + dt * ( grad_u_f - grad_u_0 );
//        phi_curr      = phi_prev      + dt * ( phi_f - phi_0 );
//        grad_phi_curr = grad_phi_prev + dt * ( grad_phi_f - grad_phi_0 );
//
//        //Update the arrays
//        for ( unsigned int i = 0; i < 3; i++ ){
//            for ( unsigned int j = 0; j < 3; j++ ){
//                current_grad_u[ i ][ j ]  = grad_u_curr[ i ][ j ];
//                previous_grad_u[ i ][ j ] = grad_u_prev[ i ][ j ];
//            }
//        }
//
//        for ( unsigned int i = 0; i < 9; i++ ){
//            current_phi[ i ] = phi_curr[ i ];
//            previous_phi[ i ] = phi_prev[ i ];
//
//            for ( unsigned int j = 0; j < 3; j++ ){
//                current_grad_phi[ i ][ j ] = grad_phi_curr[ i ][ j ];
//                previous_grad_phi[ i ][ j ] = grad_phi_prev[ i ][ j ];
//            }
//        }
//
//        //Evaluate the model
//#ifdef DEBUG_MODE
//        DEBUG.clear();
//#endif
//        int errorCode = material->evaluate_model( time, fparams,
//                                                  current_grad_u, current_phi, current_grad_phi,
//                                                  previous_grad_u, previous_phi, previous_grad_phi,
//                                                  SDVS,
//                                                  current_ADD_DOF, current_ADD_grad_DOF,
//                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                  PK2_result, SIGMA_result, M_result,
//                                                  ADD_TERMS,
//                                                  output_message
//#ifdef DEBUG_MODE
//                                                  , DEBUG
//#endif
//                                                );
//
////        std::cout << "SDVS:\n"; vectorTools::print( SDVS );
//
//#ifdef DEBUG_MODE
//
////        if ( ( fabs( SDVS[ 0 ] ) > 1e-8 ) || ( errorCode != 0 ) ){
////            for ( auto inc = DEBUG.begin(); inc != DEBUG.end(); inc++ ){
////                std::cout << inc->first << "\n";
////                for ( auto itr = inc->second.begin(); itr != inc->second.end(); itr++ ){
////                    if ( itr->first.compare( "converged_values" ) != 0 ){
////                        std::cout << "    " << itr->first << "\n";
////                        std::cout << "        currentMacroGamma:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroGamma" ] );
////                        std::cout << "        currentMicroGamma:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGamma" ] );
////                        std::cout << "        currentMicroGradientGamma:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGradientGamma" ] );
////                        std::cout << "        currentDeformationGradient:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentDeformationGradient" ] );
////                        std::cout << "        currentElasticDeformationGradient:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentElasticDeformationGradient" ] );
////                        std::cout << "        currentPlasticDeformationGradient:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentPlasticDeformationGradient" ] );
////                        std::cout << "        currentPK2Stress:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentPK2Stress" ] );
////                        std::cout << "        currentMacroCohesion:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroCohesion" ] );
////                        std::cout << "        dMacroCohesiondMacroStrainISV:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "dMacroCohesiondMacroStrainISV" ] );
////                        std::cout << "        previousMacroFlowDirection:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "previousMacroFlowDirection" ] );
////                        std::cout << "        previousMicroFlowDirection:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "previousMicroFlowDirection" ] );
////                        std::cout << "        previousMicroGradientFlowDirection:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "previousMicroGradientFlowDirection" ] );
////                        std::cout << "        currentMacroFlowDirection:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroFlowDirection" ] );
////                        std::cout << "        currentMicroFlowDirection:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroFlowDirection" ] );
////                        std::cout << "        currentMicroGradientFlowDirection:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGradientFlowDirection" ] );
////                        std::cout << "        currentMacroStrainISV\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroStrainISV" ] );
////                        std::cout << "        currentMicroStrainISV\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroStrainISV" ] );
////                        std::cout << "        currentMicroGradientStrainISV\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGradientStrainISV" ] );
////                        std::cout << "        macroYieldFunction:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "macroYieldFunction" ] );
////                        std::cout << "        microYieldFunction:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "microYieldFunction" ] );
////                        std::cout << "        microGradientYieldFunction:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "microGradientYieldFunction" ] );
////                    }
////                    else{
////                        std::cout << "        convergedPlasticDeformationGradient:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "convergedPlasticDeformationGradient" ] );
////                        std::cout << "        convergedPlasticMicroDeformation:\n";
////                        std::cout << "        "; vectorTools::print( itr->second[ "convergedPlasticMicroDeformation" ] );
////                    }
////                }
////            }
////        }
////
//        output_file << "NEW_INCREMENT\n";
//
//        //Output the time
//        output_file << time[ 0 ] << ", " << time[ 1 ] << "\n";
//
//        //Output the current gradient of u
//        output_file << current_grad_u[ 0 ][ 0 ] << ", " << current_grad_u[ 0 ][ 1 ] << ", " << current_grad_u[ 0 ][ 2 ] << ", "
//                    << current_grad_u[ 1 ][ 0 ] << ", " << current_grad_u[ 1 ][ 1 ] << ", " << current_grad_u[ 1 ][ 2 ] << ", "
//                    << current_grad_u[ 2 ][ 0 ] << ", " << current_grad_u[ 2 ][ 1 ] << ", " << current_grad_u[ 2 ][ 2 ] << "\n";
//
//        //Output the current micro displacement
//        output_file << current_phi[ 0 ] << ", " << current_phi[ 1 ] << ", " << current_phi[ 2 ] << ", "
//                    << current_phi[ 3 ] << ", " << current_phi[ 4 ] << ", " << current_phi[ 5 ] << ", "
//                    << current_phi[ 6 ] << ", " << current_phi[ 7 ] << ", " << current_phi[ 8 ] << "\n";
//
//        //Output the current gradient of the micro displacement
//        output_file << current_grad_phi[ 0 ][ 0 ] << ", " << current_grad_phi[ 0 ][ 1 ] << ", " << current_grad_phi[ 0 ][ 2 ] << ", "
//                    << current_grad_phi[ 1 ][ 0 ] << ", " << current_grad_phi[ 1 ][ 1 ] << ", " << current_grad_phi[ 1 ][ 2 ] << ", "
//                    << current_grad_phi[ 2 ][ 0 ] << ", " << current_grad_phi[ 2 ][ 1 ] << ", " << current_grad_phi[ 2 ][ 2 ] << ", "
//                    << current_grad_phi[ 3 ][ 0 ] << ", " << current_grad_phi[ 3 ][ 1 ] << ", " << current_grad_phi[ 3 ][ 2 ] << ", "
//                    << current_grad_phi[ 4 ][ 0 ] << ", " << current_grad_phi[ 4 ][ 1 ] << ", " << current_grad_phi[ 4 ][ 2 ] << ", "
//                    << current_grad_phi[ 5 ][ 0 ] << ", " << current_grad_phi[ 5 ][ 1 ] << ", " << current_grad_phi[ 5 ][ 2 ] << ", "
//                    << current_grad_phi[ 6 ][ 0 ] << ", " << current_grad_phi[ 6 ][ 1 ] << ", " << current_grad_phi[ 6 ][ 2 ] << ", "
//                    << current_grad_phi[ 7 ][ 0 ] << ", " << current_grad_phi[ 7 ][ 1 ] << ", " << current_grad_phi[ 7 ][ 2 ] << ", "
//                    << current_grad_phi[ 8 ][ 0 ] << ", " << current_grad_phi[ 8 ][ 1 ] << ", " << current_grad_phi[ 8 ][ 2 ] << "\n";
//
//        //Output the PK2 stress
//        output_file << PK2_result[ 0 ] << ", " << PK2_result[ 1 ] << ", " << PK2_result[ 2 ] << ", "
//                    << PK2_result[ 3 ] << ", " << PK2_result[ 4 ] << ", " << PK2_result[ 5 ] << ", "
//                    << PK2_result[ 6 ] << ", " << PK2_result[ 7 ] << ", " << PK2_result[ 8 ] << "\n";
//
//        //Output the SIGMA stress
//        output_file << SIGMA_result[ 0 ] << ", " << SIGMA_result[ 1 ] << ", " << SIGMA_result[ 2 ] << ", "
//                    << SIGMA_result[ 3 ] << ", " << SIGMA_result[ 4 ] << ", " << SIGMA_result[ 5 ] << ", "
//                    << SIGMA_result[ 6 ] << ", " << SIGMA_result[ 7 ] << ", " << SIGMA_result[ 8 ] << "\n";
//
//        //Output the M stress
//        output_file << M_result[  0 ] << ", " << M_result[  1 ] << ", " << M_result[  2 ] << ", "
//                    << M_result[  3 ] << ", " << M_result[  4 ] << ", " << M_result[  5 ] << ", "
//                    << M_result[  6 ] << ", " << M_result[  7 ] << ", " << M_result[  8 ] << ", "
//                    << M_result[  9 ] << ", " << M_result[ 10 ] << ", " << M_result[ 11 ] << ", "
//                    << M_result[ 12 ] << ", " << M_result[ 13 ] << ", " << M_result[ 14 ] << ", "
//                    << M_result[ 15 ] << ", " << M_result[ 16 ] << ", " << M_result[ 17 ] << ", "
//                    << M_result[ 18 ] << ", " << M_result[ 19 ] << ", " << M_result[ 20 ] << ", "
//                    << M_result[ 21 ] << ", " << M_result[ 22 ] << ", " << M_result[ 23 ] << ", "
//                    << M_result[ 24 ] << ", " << M_result[ 25 ] << ", " << M_result[ 26 ] << "\n";
//
//        //Determine if there were non-linear iterations
//        auto inc = DEBUG.find( "converged_values" );
//        if ( inc != DEBUG.end() ){
//            PK2_intermediate   = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediatePK2Stress" ];
//            SIGMA_intermediate = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceMicroStress" ];
//            M_intermediate     = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceHigherOrderStress" ];
//        }
//        else{
//            PK2_intermediate   = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediatePK2Stress" ];
//            SIGMA_intermediate = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceMicroStress" ];
//            M_intermediate     = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceHigherOrderStress" ];
//        }
//
//        //Output the intermediate PK2 stress
//        output_file << PK2_intermediate[ 0 ] << ", " << PK2_intermediate[ 1 ] << ", " << PK2_intermediate[ 2 ] << ", "
//                    << PK2_intermediate[ 3 ] << ", " << PK2_intermediate[ 4 ] << ", " << PK2_intermediate[ 5 ] << ", "
//                    << PK2_intermediate[ 6 ] << ", " << PK2_intermediate[ 7 ] << ", " << PK2_intermediate[ 8 ] << "\n";
//
//        //Output the intermediate SIGMA stress
//        output_file << SIGMA_intermediate[ 0 ] << ", " << SIGMA_intermediate[ 1 ] << ", " << SIGMA_intermediate[ 2 ] << ", "
//                    << SIGMA_intermediate[ 3 ] << ", " << SIGMA_intermediate[ 4 ] << ", " << SIGMA_intermediate[ 5 ] << ", "
//                    << SIGMA_intermediate[ 6 ] << ", " << SIGMA_intermediate[ 7 ] << ", " << SIGMA_intermediate[ 8 ] << "\n";
//
//        //Output the intermediate M stress
//        output_file << M_intermediate[  0 ] << ", " << M_intermediate[  1 ] << ", " << M_intermediate[  2 ] << ", "
//                    << M_intermediate[  3 ] << ", " << M_intermediate[  4 ] << ", " << M_intermediate[  5 ] << ", "
//                    << M_intermediate[  6 ] << ", " << M_intermediate[  7 ] << ", " << M_intermediate[  8 ] << ", "
//                    << M_intermediate[  9 ] << ", " << M_intermediate[ 10 ] << ", " << M_intermediate[ 11 ] << ", "
//                    << M_intermediate[ 12 ] << ", " << M_intermediate[ 13 ] << ", " << M_intermediate[ 14 ] << ", "
//                    << M_intermediate[ 15 ] << ", " << M_intermediate[ 16 ] << ", " << M_intermediate[ 17 ] << ", "
//                    << M_intermediate[ 18 ] << ", " << M_intermediate[ 19 ] << ", " << M_intermediate[ 20 ] << ", "
//                    << M_intermediate[ 21 ] << ", " << M_intermediate[ 22 ] << ", " << M_intermediate[ 23 ] << ", "
//                    << M_intermediate[ 24 ] << ", " << M_intermediate[ 25 ] << ", " << M_intermediate[ 26 ] << "\n";
//        
//        //Output the state variables
//        for ( unsigned int i = 0; i < SDVS.size()-1; i++ ){
//            output_file << SDVS[ i ] << ", ";
//        }
//        output_file << SDVS[ SDVS.size() - 1 ] << "\n";
//
//#endif
//
//        BOOST_CHECK( errorCode <= 0 );
//#ifdef DEBUG_MODE
//        if ( errorCode > 0 ){
//            output_file.close();
//        }
//#endif
//
//        t += dt;
//
//        grad_u_prev   = grad_u_curr;
//        phi_prev      = phi_curr;
//        grad_phi_prev = grad_phi_curr;
//
//#ifdef DEBUG_MODE    
//        if ( t > 0.25 ){ output_file.close(); return 1; }
//#endif
//    }
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( SDVSAnswer, SDVS ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( PK2Answer, PK2_result ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( SIGMAAnswer, SIGMA_result ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( MAnswer, M_result ) );
//    
//#ifdef DEBUG_MODE
//    output_file.close();
//#endif
//}

BOOST_AUTO_TEST_CASE( testComputeDruckerPragerInternalParameters ){
    /*!
     * Test the computation of the internal parameters for the Drucker-Prager plasticity.
     *
     */

    parameterType frictionAngle = 0.25;
    parameterType beta = .423;

    parameterType Aanswer = 1.528893501990677;
    parameterType Banswer = 0.39039060414065774;

    parameterType Aresult, Bresult;

    errorOut error = micromorphicElastoPlasticity::computeDruckerPragerInternalParameters( frictionAngle, beta, Aresult, Bresult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( Aresult, Aanswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( Bresult, Banswer ) );
}

BOOST_AUTO_TEST_CASE( testComputePlasticMultiplierResidual ){
    /*!
     * Test the computation of the plastic multiplier residual.
     *
     */

    std::vector< double > fparams = { 2, 1.0e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.5e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e2, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    parameterVector macroHardeningParameters;
    parameterVector microHardeningParameters;
    parameterVector microGradientHardeningParameters;

    parameterVector macroFlowParameters;
    parameterVector microFlowParameters;
    parameterVector microGradientFlowParameters;

    parameterVector macroYieldParameters;
    parameterVector microYieldParameters;
    parameterVector microGradientYieldParameters;

    parameterVector Amatrix, Bmatrix, Cmatrix, Dmatrix;
    parameterType alphaMacro, alphaMicro, alphaMicroGradient;
    parameterType relativeTolerance, absoluteTolerance;

    errorOut error = micromorphicElastoPlasticity::extractMaterialParameters( fparams,
                                                                              macroHardeningParameters, microHardeningParameters,
                                                                              microGradientHardeningParameters,
                                                                              macroFlowParameters, microFlowParameters,
                                                                              microGradientFlowParameters,
                                                                              macroYieldParameters, microYieldParameters,
                                                                              microGradientYieldParameters,
                                                                              Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                              alphaMacro, alphaMicro, alphaMicroGradient,
                                                                              relativeTolerance, absoluteTolerance );

    BOOST_CHECK( !error );

    constantType   Dt = 2.5;
    variableType   currentMacroGamma = 0.01;
    variableType   currentMicroGamma = 0.02;
    variableVector currentMicroGradientGamma = {0.011, 0.021, 0.031};
    variableType   previousMacroGamma = 0.;
    variableType   previousMicroGamma = 0.;
    variableVector previousMicroGradientGamma( 3, 0 );
    variableType   previousMacroStrainISV = 0.;
    variableType   previousMicroStrainISV = 0.;
    variableVector previousMicroGradientStrainISV( 3., 0. );
    variableVector currentDeformationGradient = { 0.04482969,  0.88562312, -0.38710144,
                                                 -0.93722716,  0.19666568,  0.41677155,
                                                  0.46929057,  0.33779672,  0.81392228 };
    variableVector currentMicroDeformation = {  0.51930689,  0.27954023, -0.85955731,
                                                0.09469279, -0.99381243, -0.23218079,
                                               -0.82281393,  0.09643296, -0.54637704 };
    variableVector currentGradientMicroDeformation = { 0.04176306, -0.0151958 , -0.00090558, -0.01844751,  0.04512391,
                                                       0.02174263, -0.00508239, -0.01827377,  0.00541031, -0.01330239,
                                                      -0.02479987,  0.02914825,  0.00168841,  0.00230506,  0.00994845,
                                                       0.00413116,  0.04555686, -0.00431862, -0.0138286 , -0.04412473,
                                                       0.02016718, -0.03868735,  0.03842166, -0.0009337 ,  0.02977617,
                                                       0.02310445,  0.02827616 };
    variableVector previousPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticGradientMicroDeformation = variableVector( 27, 0 );
    variableVector previousPlasticMacroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroGradientVelocityGradient = variableVector( 27, 0 );

    variableVector currentElasticDeformationGradient = currentDeformationGradient;
    variableVector currentElasticMicroDeformation = currentMicroDeformation;
    variableVector currentElasticGradientMicroDeformation = currentGradientMicroDeformation;
    variableVector currentPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticGradientMicroDeformation = variableVector( 27, 0 );

    variableType   previousdMacroGdMacroCohesion = 0;
    variableType   previousdMicroGdMicroCohesion = 0;
    variableMatrix previousdMicroGradientGdMicroGradientCohesion( 3, variableVector( 3, 0 ) );

    variableType currentMacroStrainISV = 0;
    variableType currentMicroStrainISV = 0;
    variableVector currentMicroGradientStrainISV( 3, 0 );

    BOOST_CHECK( !error );

    solverTools::floatMatrix floatArgsDefault =
        {
            { Dt },
            currentDeformationGradient,
            currentMicroDeformation,
            currentGradientMicroDeformation,
//            { currentMacroGamma },
//            { currentMicroGamma },
//            currentMicroGradientGamma,
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

    variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

    solverTools::floatMatrix floatOutsDefault = 
        {
            currentPK2Stress,
            currentReferenceMicroStress,
            currentReferenceHigherOrderStress,
            { currentMacroStrainISV },
            { currentMicroStrainISV },
            currentMicroGradientStrainISV,
            {}, {}, {}
        };

    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intArgs = { { 0 } };
    solverTools::intMatrix intOutsDefault = { { 0, 0 } };
    
    solverTools::intMatrix intOuts = intOutsDefault;

#ifdef DEBUG_MODE
    solverTools::debugMap DEBUG;
#endif

    solverTools::floatVector x =
        {
            currentMacroGamma,
            currentMicroGamma,
            currentMicroGradientGamma[ 0 ],
            currentMicroGradientGamma[ 1 ],
            currentMicroGradientGamma[ 2 ]
        };

    solverTools::floatVector residualAnswer = { 319.6905968 , 391.17264214, -2.92818166, -5.60542895, -8.29895044 };

    solverTools::floatVector residual;
    solverTools::floatMatrix jacobian;

    error = micromorphicElastoPlasticity::computePlasticMultiplierResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                            floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                            , DEBUG
#endif
                                                                          );
    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( residual, residualAnswer ) ); 

    //Tests of the Jacobians
    solverTools::floatType eps = 1e-6;
    for ( unsigned int i = 0; i < x.size(); i++ ){
        solverTools::floatVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector rP, rM;
        solverTools::floatMatrix _J;

        solverTools::floatMatrix fOP, fOM;
        fOP = floatOutsDefault;
        fOM = floatOutsDefault;

        solverTools::intMatrix iOP, iOM;
        iOP = intOutsDefault;
        iOM = intOutsDefault;

#ifdef DEBUG_MODE
        solverTools::debugMap DEBUG_P, DEBUG_M;
#endif

        error =  micromorphicElastoPlasticity::computePlasticMultiplierResidual( x + delta, floatArgs, intArgs, rP, _J,
                                                                                 fOP, iOP
#ifdef DEBUG_MODE
                                                                                 , DEBUG_P
#endif
                                                                               );

        BOOST_CHECK( !error );

        error =  micromorphicElastoPlasticity::computePlasticMultiplierResidual( x - delta, floatArgs, intArgs, rM, _J,
                                                                                 fOM, iOM
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        solverTools::debugMap numericGradients;
        for ( auto itP = DEBUG_P.begin(); itP != DEBUG_P.end(); itP++ ){
            auto itM = DEBUG_M.find( itP->first );
            BOOST_CHECK( itM != DEBUG_M.end() );

            numericGradients.emplace( itP->first, ( itP->second - itM->second ) / ( 2 * delta[ i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "currentPlasticDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticDeformation" ][ j ],
                                            DEBUG[ "dCurrentPlasticDeformationdGammas" ][ 5 * j + i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "stresses" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "stresses" ][ j ],
                                            DEBUG[ "dStressdGammas" ][ 5 * j + i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                            DEBUG[ "dElasticRightCauchyGreendGammas" ][ 5 * j + i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "yieldFunctionValues" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "yieldFunctionValues" ][ j ],
                                            DEBUG[ "dYieldFunctionValuesdGammas" ][ 5 * j + i ] ) );
        }
#endif
        //Test of the residual and comparison to the Jacobian
        solverTools::floatVector gradCol = ( rP - rM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticMultiplierResidual2 ){
    /*!
     * Test the computation of the plastic multiplier residual.
     *
     */

    std::vector< double > fparams = { 2, 1.0e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.5e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e2, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    parameterVector macroHardeningParameters;
    parameterVector microHardeningParameters;
    parameterVector microGradientHardeningParameters;

    parameterVector macroFlowParameters;
    parameterVector microFlowParameters;
    parameterVector microGradientFlowParameters;

    parameterVector macroYieldParameters;
    parameterVector microYieldParameters;
    parameterVector microGradientYieldParameters;

    parameterVector Amatrix, Bmatrix, Cmatrix, Dmatrix;
    parameterType alphaMacro, alphaMicro, alphaMicroGradient;
    parameterType relativeTolerance, absoluteTolerance;

    errorOut error = micromorphicElastoPlasticity::extractMaterialParameters( fparams,
                                                                              macroHardeningParameters, microHardeningParameters,
                                                                              microGradientHardeningParameters,
                                                                              macroFlowParameters, microFlowParameters,
                                                                              microGradientFlowParameters,
                                                                              macroYieldParameters, microYieldParameters,
                                                                              microGradientYieldParameters,
                                                                              Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                              alphaMacro, alphaMicro, alphaMicroGradient,
                                                                              relativeTolerance, absoluteTolerance );

    BOOST_CHECK( !error );

    constantType   Dt = 2.5;
    variableType   currentMacroGamma = 0.01;
    variableType   currentMicroGamma = 0.02;
    variableVector currentMicroGradientGamma = {0.011, 0.021, 0.031};
    variableType   previousMacroGamma = 0.;
    variableType   previousMicroGamma = 0.;
    variableVector previousMicroGradientGamma( 3, 0 );
    variableType   previousMacroStrainISV = 0.;
    variableType   previousMicroStrainISV = 0.;
    variableVector previousMicroGradientStrainISV( 3., 0. );
    variableVector currentDeformationGradient = { 0.04482969,  0.88562312, -0.38710144,
                                                 -0.93722716,  0.19666568,  0.41677155,
                                                  0.46929057,  0.33779672,  0.81392228 };
    variableVector currentMicroDeformation = {  0.51930689,  0.27954023, -0.85955731,
                                                0.09469279, -0.99381243, -0.23218079,
                                               -0.82281393,  0.09643296, -0.54637704 };
    variableVector currentGradientMicroDeformation = { 0.04176306, -0.0151958 , -0.00090558, -0.01844751,  0.04512391,
                                                       0.02174263, -0.00508239, -0.01827377,  0.00541031, -0.01330239,
                                                      -0.02479987,  0.02914825,  0.00168841,  0.00230506,  0.00994845,
                                                       0.00413116,  0.04555686, -0.00431862, -0.0138286 , -0.04412473,
                                                       0.02016718, -0.03868735,  0.03842166, -0.0009337 ,  0.02977617,
                                                       0.02310445,  0.02827616 };
    variableVector previousPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticGradientMicroDeformation = variableVector( 27, 0 );
    variableVector previousPlasticMacroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroGradientVelocityGradient = variableVector( 27, 0 );

    variableVector currentElasticDeformationGradient = currentDeformationGradient;
    variableVector currentElasticMicroDeformation = currentMicroDeformation;
    variableVector currentElasticGradientMicroDeformation = currentGradientMicroDeformation;
    variableVector currentPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticGradientMicroDeformation = variableVector( 27, 0 );

    variableType   previousdMacroGdMacroCohesion = 0;
    variableType   previousdMicroGdMicroCohesion = 0;
    variableMatrix previousdMicroGradientGdMicroGradientCohesion( 3, variableVector( 3, 0 ) );

    variableType currentMacroStrainISV = 0;
    variableType currentMicroStrainISV = 0;
    variableVector currentMicroGradientStrainISV( 3, 0 );

    BOOST_CHECK( !error );

    solverTools::floatMatrix floatArgsDefault =
        {
            { Dt },
            currentDeformationGradient,
            currentMicroDeformation,
            currentGradientMicroDeformation,
//            { currentMacroGamma },
//            { currentMicroGamma },
//            currentMicroGradientGamma,
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

    variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

    solverTools::floatMatrix floatOutsDefault = 
        {
            currentPK2Stress,
            currentReferenceMicroStress,
            currentReferenceHigherOrderStress,
            { currentMacroStrainISV },
            { currentMicroStrainISV },
            currentMicroGradientStrainISV,
            {}, {}, {},
            {}, {}, {}, {}, {}
        };

    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intArgs = { { 1 } };
    solverTools::intMatrix intOutsDefault = { { 0, 0 } };
    
    solverTools::intMatrix intOuts = intOutsDefault;

#ifdef DEBUG_MODE
    solverTools::debugMap DEBUG;
#endif

    solverTools::floatVector x =
        {
            currentMacroGamma,
            currentMicroGamma,
            currentMicroGradientGamma[ 0 ],
            currentMicroGradientGamma[ 1 ],
            currentMicroGradientGamma[ 2 ]
        };

    solverTools::floatVector residualAnswer = { 319.6905968 , 391.17264214, -2.92818166, -5.60542895, -8.29895044 };

    solverTools::floatVector residual;
    solverTools::floatMatrix jacobian;

    error = micromorphicElastoPlasticity::computePlasticMultiplierResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                            floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                            , DEBUG
#endif
                                                                          );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( residual, residualAnswer ) ); 

    //Tests of the Jacobians
    solverTools::floatType eps = 1e-6;
    for ( unsigned int i = 0; i < x.size(); i++ ){
        solverTools::floatVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector rP, rM;
        solverTools::floatMatrix _J;

        solverTools::floatMatrix fOP, fOM;
        fOP = floatOutsDefault;
        fOM = floatOutsDefault;

        solverTools::intMatrix iOP, iOM;
        iOP = intOutsDefault;
        iOM = intOutsDefault;

#ifdef DEBUG_MODE
        solverTools::debugMap DEBUG_P, DEBUG_M;
#endif

        error =  micromorphicElastoPlasticity::computePlasticMultiplierResidual( x + delta, floatArgs, intArgs, rP, _J,
                                                                                 fOP, iOP
#ifdef DEBUG_MODE
                                                                                 , DEBUG_P
#endif
                                                                               );

        BOOST_CHECK( !error );

        error =  micromorphicElastoPlasticity::computePlasticMultiplierResidual( x - delta, floatArgs, intArgs, rM, _J,
                                                                                 fOM, iOM
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        solverTools::debugMap numericGradients;
        for ( auto itP = DEBUG_P.begin(); itP != DEBUG_P.end(); itP++ ){
            auto itM = DEBUG_M.find( itP->first );
            BOOST_CHECK( itM != DEBUG_M.end() );

            numericGradients.emplace( itP->first, ( itP->second - itM->second ) / ( 2 * delta[ i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "currentPlasticDeformation" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticDeformation" ][ j ],
                                            DEBUG[ "dCurrentPlasticDeformationdGammas" ][ 5 * j + i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "stresses" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "stresses" ][ j ],
                                            DEBUG[ "dStressdGammas" ][ 5 * j + i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                            DEBUG[ "dElasticRightCauchyGreendGammas" ][ 5 * j + i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "yieldFunctionValues" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "yieldFunctionValues" ][ j ],
                                            DEBUG[ "dYieldFunctionValuesdGammas" ][ 5 * j + i ] ) );
        }
#endif
        //Test of the residual and comparison to the Jacobian
        solverTools::floatVector gradCol = ( rP - rM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ] ) );
        }

        //Test of the partial derivatives of the plastic deformation w.r.t. the plastic multipliers
        gradCol = vectorTools::appendVectors( { fOP[ 6 ] - fOM[ 6 ], fOP[ 7 ] - fOM[ 7 ], fOP[ 8 ] - fOM[ 8 ] } ) / ( 2 * delta[ i ] );
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 9 ][ 5 * j + i ] ) );
        }

        //Test of the partial derivatives of the Stresses w.r.t. the plastic multipliers
        gradCol = vectorTools::appendVectors( { fOP[ 0 ] - fOM[ 0 ], fOP[ 1 ] - fOM[ 1 ], fOP[ 2 ] - fOM[ 2 ] } ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 10 ][ 5 * j + i ] ) );
        }
    }

    //Test the Jacobians w.r.t. the deformation measures
    solverTools::floatVector deformation = vectorTools::appendVectors( { currentDeformationGradient,
                                                                         currentMicroDeformation,
                                                                         currentGradientMicroDeformation } );

    for ( unsigned int i = 0; i < deformation.size(); i++ ){
        solverTools::floatVector delta( deformation.size(), 0 );
        delta[ i ] = eps * fabs( deformation[ i ] ) + eps;

        solverTools::floatVector rP, rM;
        solverTools::floatMatrix _J;

        solverTools::floatMatrix fA_P, fA_M;
        fA_P = floatArgs;
        fA_M = floatArgs;

        fA_P[ 1 ] += solverTools::floatVector( delta.begin() +  0, delta.begin() +  9 );
        fA_P[ 2 ] += solverTools::floatVector( delta.begin() +  9, delta.begin() + 18 );
        fA_P[ 3 ] += solverTools::floatVector( delta.begin() + 18, delta.begin() + 45 );

        fA_M[ 1 ] -= solverTools::floatVector( delta.begin() +  0, delta.begin() +  9 );
        fA_M[ 2 ] -= solverTools::floatVector( delta.begin() +  9, delta.begin() + 18 );
        fA_M[ 3 ] -= solverTools::floatVector( delta.begin() + 18, delta.begin() + 45 );

        solverTools::floatMatrix fOP, fOM;
        fOP = floatOutsDefault;
        fOM = floatOutsDefault;

        solverTools::intMatrix iOP, iOM;
        iOP = intOutsDefault;
        iOM = intOutsDefault;

#ifdef DEBUG_MODE
        solverTools::debugMap DEBUG_P, DEBUG_M;
#endif

        error =  micromorphicElastoPlasticity::computePlasticMultiplierResidual( x, fA_P, intArgs, rP, _J,
                                                                                 fOP, iOP
#ifdef DEBUG_MODE
                                                                                 , DEBUG_P
#endif
                                                                               );

        BOOST_CHECK( !error );

        error =  micromorphicElastoPlasticity::computePlasticMultiplierResidual( x, fA_M, intArgs, rM, _J,
                                                                                 fOM, iOM
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        BOOST_CHECK( !error );

#ifdef DEBUG_MODE
        solverTools::debugMap numericGradients;
        for ( auto itP = DEBUG_P.begin(); itP != DEBUG_P.end(); itP++ ){
            auto itM = DEBUG_M.find( itP->first );

            BOOST_CHECK( itM != DEBUG_M.end() );

            numericGradients.emplace( itP->first, ( itP->second - itM->second ) / ( 2 * delta[ i ] ) );
        }

        //Test the gradient of the plastic deformation w.r.t. the deformation
        for ( unsigned int j = 0; j < numericGradients[ "currentPlasticDeformation" ].size(); j++ ){
            
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentPlasticDeformation" ][ j ],
                                            DEBUG[ "dCurrentPlasticDeformationdDeformation" ][ 45 * j + i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                            DEBUG[ "dElasticRightCauchyGreendDeformation" ][ 45 * j + i ] ) );
        }

        for ( unsigned int j = 0; j < numericGradients[ "yieldFunctionValues" ].size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( numericGradients[ "yieldFunctionValues" ][ j ],
                                            DEBUG[ "dYieldFunctionValuesdDeformation" ][ 45 * j + i ] ) );
        }
#endif

        solverTools::floatVector gradCol;
        //Check the gradient of the plastic deformation w.r.t. the deformation measures
        gradCol = vectorTools::appendVectors( { fOP[ 6 ] - fOM[ 6 ], fOP[ 7 ] - fOM[ 7 ], fOP[ 8 ] - fOM[ 8 ] } );
        gradCol /= 2 * delta[ i ];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 11 ][ 45 * j + i ] ) );
        }

        //Check the gradient of the stresses w.r.t. the deformation measures
        gradCol = vectorTools::appendVectors( { fOP[ 0 ] - fOM[ 0 ], fOP[ 1 ] - fOM[ 1 ], fOP[ 2 ] - fOM[ 2 ] } );
        gradCol /= 2 * delta[ i ];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 12 ][ 45 * j + i ] ) );
        }

        //Check the gradient of the residual w.r.t. the deformation measures
        gradCol = ( rP - rM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 13 ][ 45 * j + i ] ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testHomotopySolve ){
    /*!
     * Test of the homotopy solver with a problematic 
     * timestep observed in application.
     *
     */

    solverTools::floatVector x0 = { 0, 0, 0, 0, 0 };

    solverTools::floatMatrix floatArgs =
       {
           { 0.005000, },
           { 0.988293, 0.000000, -0.000000, 0.000000, 1.037000, -0.000000, 0.000000, 0.000000, 0.988293, },
           { 0.984408, 0.000000, -0.000000, -0.000000, 1.030228, -0.000000, 0.000000, 0.000000, 0.984408, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -0.000000, -0.000000, 0.000000, 0.000000, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, },
           { 0.000000, },
           { 0.000000, },
           { 0.000000, 0.000000, 0.000000, },
           { 1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000, },
           { 1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, },
           { 0.000000, },
           { 0.000000, 0.000000, 0.000000, },
           { 0.000000, },
           { 0.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 1000.000000, 100.000000, },
           { 700.000000, 10000.000000, },
           { 1000.000000, 10000.000000, },
           { 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 80440.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 25480.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 80440.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 25480.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 80440.000000, },
           { 3700.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 1500.000000, 0.000000, 1600.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1500.000000, 0.000000, 0.000000, 0.000000, 1600.000000, 0.000000, 0.000000, 0.000000, 1600.000000, 0.000000, 1500.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 3700.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1500.000000, 0.000000, 1600.000000, 0.000000, 0.000000, 0.000000, 1600.000000, 0.000000, 0.000000, 0.000000, 1500.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1600.000000, 0.000000, 1500.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 3700.000000, },
           { 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, },
           { -5600.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, -3000.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, -5600.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, -3000.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, -5600.000000, },
           { 0.500000, },
           { 0.500000, },
           { 0.500000, }
       };

    solverTools::intMatrix intArgs = { { 1 } };

    solverTools::floatMatrix floatOuts =
        {
            {}, {}, {},
            {}, {}, {},
            {}, {}, {},
            {}, {}, {}, {}, {}            
        };

    solverTools::intMatrix intOuts = { { } };

    solverTools::solverType linearSolver;
    solverTools::floatMatrix J;
    solverTools::intVector boundVariableIndices = {  0,  1,  2,  3,  4 };
    solverTools::intVector boundSigns           = {  0,  0,  0,  0,  0 };
    solverTools::floatVector boundValues        = {  0,  0,  0,  0,  0 };

    solverTools::stdFncNLFJ func
                = static_cast< solverTools::NonLinearFunctionWithJacobian >( micromorphicElastoPlasticity::computePlasticMultiplierResidual );

    solverTools::floatVector solutionVector;

    bool convergeFlag, fatalErrorFlag;
    solverTools::floatType relativeTolerance = 1e-9;
    solverTools::floatType absoluteTolerance = 1e-9;

#ifdef DEBUG_MODE
    solverTools::homotopyMap DEBUG;
#endif

    solverTools::floatVector answer = { 0, 2.91946515, 0, 0, 0 };

    errorOut error = solverTools::homotopySolver( func, x0, solutionVector, convergeFlag, fatalErrorFlag,
                                                  floatOuts, intOuts, floatArgs, intArgs, linearSolver, J,
                                                  boundVariableIndices, boundSigns, boundValues, true,
#ifdef DEBUG_MODE
                                                  DEBUG,
#endif
                                                  20, relativeTolerance, absoluteTolerance,
                                                  1e-4, 5, 1.0, .01 );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, solutionVector ) );

}

BOOST_AUTO_TEST_CASE( testComputePlasticMultiplierLagrangian ){
    /*!
     * Test the computation of the Lagrangian of the plastic multiplier.
     *
     */

    solverTools::floatVector x0 = { 0.1, 0.01, 0.02, 0.03, 0.04, 1, 1, 1, 1, 1 };

    solverTools::floatMatrix floatArgs =
       {
           { 0.005000, },
           { 0.988293, 0.000000, -0.000000, 0.000000, 1.037000, -0.000000, 0.000000, 0.000000, 0.988293, },
           { 0.984408, 0.000000, -0.000000, -0.000000, 1.030228, -0.000000, 0.000000, 0.000000, 0.984408, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -0.000000, -0.000000, 0.000000, 0.000000, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, },
           { 0.000000, },
           { 0.000000, },
           { 0.000000, 0.000000, 0.000000, },
           { 1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000, },
           { 1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, },
           { 0.000000, },
           { 0.000000, 0.000000, 0.000000, },
           { 0.000000, },
           { 0.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 1000.000000, 100.000000, },
           { 700.000000, 10000.000000, },
           { 1000.000000, 10000.000000, },
           { 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 0.000000, 0.000000, },
           { 80440.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 25480.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 80440.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 25480.000000, 0.000000, 25480.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 29480.000000, 0.000000, 0.000000, 0.000000, 80440.000000, },
           { 3700.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 1500.000000, 0.000000, 1600.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1500.000000, 0.000000, 0.000000, 0.000000, 1600.000000, 0.000000, 0.000000, 0.000000, 1600.000000, 0.000000, 1500.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 3700.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1500.000000, 0.000000, 1600.000000, 0.000000, 0.000000, 0.000000, 1600.000000, 0.000000, 0.000000, 0.000000, 1500.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1600.000000, 0.000000, 1500.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 600.000000, 0.000000, 0.000000, 0.000000, 3700.000000, },
           { 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1000000.000000, },
           { -5600.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, -3000.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, -5600.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -3000.000000, 0.000000, -3000.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, 400.000000, 0.000000, 0.000000, 0.000000, -5600.000000, },
           { 0.500000, },
           { 0.500000, },
           { 0.500000, }
       };

    solverTools::intMatrix intArgs = { };

    solverTools::floatMatrix floatOuts = { };

    solverTools::intMatrix intOuts = { { } };

    solverTools::floatType lagrangian;
    solverTools::floatType lagrangian_answer = 117267.160955;
    solverTools::floatVector jacobian;

#ifdef DEBUG_MODE
    solverTools::debugMap DEBUG;
#endif

    errorOut error = micromorphicElastoPlasticity::computePlasticMultiplierLagrangian( x0, floatArgs, intArgs, lagrangian, jacobian,
                                                                                       floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                                       , DEBUG
#endif
                                                                                     );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( lagrangian, lagrangian_answer ) );

    //Check the Jacobians
    solverTools::floatType eps = 1e-6;
    for ( unsigned int i = 0; i < x0.size(); i++ ){
        solverTools::floatVector delta( x0.size(), 0 );
        delta[ i ] = eps * fabs( x0[ i ] ) + eps;

        solverTools::floatType lP, lM;
        solverTools::floatVector JP, JM;

        error = micromorphicElastoPlasticity::computePlasticMultiplierLagrangian( x0 + delta, floatArgs, intArgs, lP, JP,
                                                                                  floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                                  , DEBUG
#endif
                                                                                );

        BOOST_CHECK( !error );

        error = micromorphicElastoPlasticity::computePlasticMultiplierLagrangian( x0 - delta, floatArgs, intArgs, lM, JM,
                                                                                  floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                                  , DEBUG
#endif
                                                                                );

        BOOST_CHECK( !error );

        solverTools::floatType grad = ( lP - lM ) / ( 2 * delta[ i ] );
        BOOST_CHECK( vectorTools::fuzzyEquals( grad, jacobian[ i ] ) );
    }

}
