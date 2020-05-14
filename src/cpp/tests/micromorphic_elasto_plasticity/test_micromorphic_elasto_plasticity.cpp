//Tests for constitutive_tools

#include<micromorphic_elasto_plasticity.h>
#include<sstream>
#include<fstream>
#include<iostream>
#include<iomanip>

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

int test_computeSecondOrderDruckerPragerYieldEquation( std::ofstream &results ){
    /*!
     * Test the computation of the second order stress Drucker-Prager yield equation.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_computeSecondOrderDruckerPragerYieldEquation (test 1) & False\n";
        return 1;
    }

    //Test the Jacobian
    variableType resultJ, dFdCohesion;
    variableVector dFdS, dFdC;

    error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, C,
                                                                                        frictionAngle, beta, resultJ,
                                                                                        dFdS, dFdCohesion, dFdC );

    if ( error ){
        error->print();
        results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_computeSecondOrderDruckerPragerYieldEquation (test 2) & False\n";
        return 1;
    }


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

    if ( error ){
        error->print();
        results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ2, answer ) ){
        results << "test_computeSecondOrderDruckerPragerYieldEquation (test 3) & False\n";
        return 1;
    }

    //Test dFdStress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S + delta, cohesion, C,
                                                                                            frictionAngle, beta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        constantType gradCol = ( resultJ - result ) / delta[i];

        if ( !vectorTools::fuzzyEquals( gradCol, dFdS[i] ) ){
            results << "test_computeSecondOrderDruckerPragerYieldEquation (test 4) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol, dFdSJ2[i] ) ){
            results << "test_computeSecondOrderDruckerPragerYieldEquation (test 5) & False\n";
            return 1;
        }
    }

    //Test dFdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, C + delta,
                                                                                            frictionAngle, beta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        constantType gradCol = ( resultJ - result ) / delta[i];

        if ( !vectorTools::fuzzyEquals( gradCol, dFdC[i], 1e-4 ) ){
            results << "test_computeSecondOrderDruckerPragerYieldEquation (test 6) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol, dFdCJ2[i], 1e-4 ) ){
            results << "test_computeSecondOrderDruckerPragerYieldEquation (test 7) & False\n";
            return 1;
        }
    }

    //Test dFdcohesion
    constantType deltas = eps * fabs( cohesion ) + eps;

    error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion + deltas, C, 
                                                                                        frictionAngle, beta, resultJ );

    if ( error ){
        error->print();
        results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ( resultJ - result ) / deltas, dFdcohesionJ2 ) ){
        results << "test_computeSecondOrderDruckerPragerYieldEquation (test 8) & False\n";
        return 1;
    }

    //Test d2FdStress2
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        variableVector dFdSp, dFdSm, dFdCp, dFdCm;
        variableType dFdcohesionp, dFdcohesionm;

        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S + delta, cohesion, C,
                                                                                            frictionAngle, beta, resultJ,
                                                                                            dFdSp, dFdcohesionp, dFdCp );

        if ( error ){
            error->print();
            results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }
                
        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S - delta, cohesion, C,
                                                                                            frictionAngle, beta, resultJ,
                                                                                            dFdSm, dFdcohesionm, dFdCm );

        if ( error ){
            error->print();
            results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        constantVector gradCol = ( dFdSp - dFdSm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], d2FdS2J2[j][i] ) ){
                results << "test_computeSecondOrderDruckerPragerYieldEquation (test 9) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }
                
        error = micromorphicElastoPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, C - delta,
                                                                                            frictionAngle, beta, resultJ,
                                                                                            dFdSm, dFdcohesionm, dFdCm );

        if ( error ){
            error->print();
            results << "test_computeSecondOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        constantVector gradCol = ( dFdSp - dFdSm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], d2FdSdCJ2[j][i] ) ){
                results << "test_computeSecondOrderDruckerPragerYieldEquation (test 10) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeSecondOrderDruckerPragerYieldEquation & True\n";
    return 0;
}

int test_computeHigherOrderDruckerPragerYieldEquation( std::ofstream &results ){
    /*!
     * Test the computation of the higher order stress Drucker-Prager yield equation.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_computeHigherOrderDruckerPragerYieldEquation (test 1) & False\n";
        return 1;
    }

    //Test the Jacobians

    variableVector resultJ;
    variableMatrix dFdStress, dFdc, dFdRCG;

    error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C,
                                                                                        frictionAngle, beta, resultJ,
                                                                                        dFdStress, dFdc, dFdRCG );

    if ( error ){
        error->print();
        results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_computeHigherOrderDruckerPragerYieldEquation (test 2) & False\n";
        return 1;
    }

    variableVector resultJ2;
    variableMatrix dFdStressJ2, dFdcJ2, dFdRCGJ2;
    variableMatrix d2FdStress2J2, d2FdStressdRCGJ2;

    error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C,
                                                                                        frictionAngle, beta, resultJ2,
                                                                                        dFdStressJ2, dFdcJ2, dFdRCGJ2,
                                                                                        d2FdStress2J2, d2FdStressdRCGJ2 );

    if ( error ){
        error->print();
        results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ2, answer ) ){
        results << "test_computeHigherOrderDruckerPragerYieldEquation (test 3) & False\n";
        return 1;
    }
    
    //Test derivatives w.r.t stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M + delta, cohesion, C,
                                                                                            frictionAngle, beta, resultP );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M - delta, cohesion, C,
                                                                                            frictionAngle, beta, resultM );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        constantVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dFdStress[j][i] ) ){
                std::cout << "i, j: " << i << ", " << j << "\n";
                std::cout << "gradCol:\n"; vectorTools::print( gradCol );
                std::cout << "dFdStress:\n"; vectorTools::print( dFdStress );
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 4) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dFdStressJ2[j][i] ) ){
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 5) & False\n";
                return 1;
            }
        }

        variableMatrix dFdStressP, dFdcP, dFdRCGP;
        variableMatrix dFdStressM, dFdcM, dFdRCGM;

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M + delta, cohesion, C,
                                                                                            frictionAngle, beta, resultP,
                                                                                            dFdStressP, dFdcP, dFdRCGP );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M - delta, cohesion, C,
                                                                                            frictionAngle, beta, resultM,
                                                                                            dFdStressM, dFdcM, dFdRCGM );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }


        constantMatrix gradMat = ( dFdStressP - dFdStressM ) / ( 2 * delta[i] );

        unsigned int n, o, p;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        n = ( int )( i / 9 );
                        o = ( int )( (i - 9 * n ) / 3 );
                        p = ( i - 9 * n - 3 * o ) % 3;
                        if ( !vectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                        d2FdStress2J2[ j ][ 243 * k + 81 * l + 27 * m + 9 * n + 3 * o + p ] ) ){
                            results << "test_computeHigherOrderDruckerPragerYieldEquation (test 6) & False\n";
                            return 1;
                        }
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

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion - delta, C,
                                                                                            frictionAngle, beta, resultM );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        constantVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dFdc[j][i] ) ){
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 7) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dFdcJ2[j][i] ) ){
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 8) & False\n";
                return 1;
            }
        }
    }

    //Test derivatives w.r.t. the right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C + delta,
                                                                                            frictionAngle, beta, resultP );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C - delta,
                                                                                            frictionAngle, beta, resultM );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        constantVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dFdRCG[j][i] ) ){
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 9) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dFdRCG[j][i] ) ){
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 10) & False\n";
                return 1;
            }
        }

        variableMatrix dFdStressP, dFdcP, dFdRCGP;
        variableMatrix dFdStressM, dFdcM, dFdRCGM;

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C + delta,
                                                                                            frictionAngle, beta, resultP,
                                                                                            dFdStressP, dFdcP, dFdRCGP );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, C - delta,
                                                                                            frictionAngle, beta, resultM,
                                                                                            dFdStressM, dFdcM, dFdRCGM );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderDruckerPragerYieldEquation & False\n";
            return 1;
        }


        constantMatrix gradMat = ( dFdStressP - dFdStressM ) / ( 2 * delta[i] );

        unsigned int n, o;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        n = ( int )( i / 3 );
                        o = ( i % 3 );
                        if ( !vectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                        d2FdStressdRCGJ2[ j ][ 81 * k + 27 * l + 9 * m + 3 * n + o ] ) ){
                            results << "test_computeHigherOrderDruckerPragerYieldEquation (test 11) & False\n";
                            return 1;
                        }
                    }
                }
            }
        }
    }

    results << "test_computeHigherOrderDruckerPragerYieldEquation & True\n";
    return 0;
}

int test_computeElasticPartOfDeformation( std::ofstream &results ){
    /*!
     * Test of the computation of the elastic part of the various deformation measures.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_computeElasticPartOfDeformation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultFe, answerFe ) ){
        results << "test_computeElasticPartOfDeformation (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultChie, answerChie ) ){
        results << "test_computeElasticPartOfDeformation (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultGradChie, answerGradChie ) ){
        results << "test_computeElasticPartOfDeformation (test 3) & False\n";
        return 1;
    }

    variableVector resultFe2, resultChie2, resultGradChie2, resultInvFp, resultInvChip;

    error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi, 
                                                                            Fp, chip, gradChip,
                                                                            resultInvFp, resultInvChip,
                                                                            resultFe2, resultChie2,
                                                                            resultGradChie2 );

    if ( error ){
        error->print();
        results << "test_computeElasticPartOfDeformation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultFe2, answerFe ) ){
        results << "test_computeElasticPartOfDeformation (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultChie2, answerChie ) ){
        results << "test_computeElasticPartOfDeformation (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultGradChie2, answerGradChie ) ){
        results << "test_computeElasticPartOfDeformation (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultInvFp, vectorTools::inverse( Fp, 3, 3 ) ) ){
        results << "test_computeElasticPartOfDeformation (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultInvChip, vectorTools::inverse( chip, 3, 3 ) ) ){
        results << "test_computeElasticPartOfDeformation (test 8) & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_computeElasticPartOfDeformation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultFeJ, answerFe ) ){
        results << "test_computeElasticPartOfDeformation (test 9) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultChieJ, answerChie ) ){
        results << "test_computeElasticPartOfDeformation (test 10) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultGradChieJ, answerGradChie ) ){
        results << "test_computeElasticPartOfDeformation (test 11) & False\n";
        return 1;
    }

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

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F - delta,  chi,  gradChi,
                                                                                Fp, chip, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        //Test dFedF
        variableVector gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dFedF[ j ][ i ] ) ){
                results << "test_computeElasticPartOfDeformation (test 12) & False\n";
                return 1;
            }
        }

        //Test dChiedF
        gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticPartOfDeformation (test 13) & False\n";
            return 1;
        }

        //Test dGradChiedF
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticPartOfDeformation (test 14) & False\n";
            return 1;
        }
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

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi,
                                                                                Fp - delta, chip, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        variableVector gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dFedFp[ j ][ i ] ) ){
                results << "test_computeElasticPartOfDeformation (test 15) & False\n";
                return 1;
            }
        }

        //Test dChiedFp
        gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticPartOfDeformation (test 16) & False\n";
            return 1;
        }

        //Test dGradChiedFp
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dGradChiedFp[j][i] ) ){
                results << "test_computeElasticPartOfDeformation (test 17) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi - delta,  gradChi,
                                                                                Fp, chip, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        //Test dChiedChi
        variableVector gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dChiedChi[ j ][ i ] ) ){
                results << "test_computeElasticPartOfDeformation (test 18) & False\n";
                return 1;
            }
        }

        //Test dFedChi
        gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticPartOfDeformation (test 19) & False\n";
            return 1;
        }

        //Test dGradChiedChi
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dGradChiedChi[j][i] ) ){
                results << "test_computeElasticPartOfDeformation (test 20) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi,
                                                                                Fp, chip - delta, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        variableVector gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dChiedChip[ j ][ i ] ) ){
                results << "test_computeElasticPartOfDeformation (test 21) & False\n";
                return 1;
            }
        }

        //Test dFedFp
        gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticPartOfDeformation (test 22) & False\n";
            return 1;
        }

        //Test dGradChiedChip
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dGradChiedChip[j][i] ) ){
                results << "test_computeElasticPartOfDeformation (test 23) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi - delta,
                                                                                Fp, chip, gradChip,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        //Test dFedGradChi
        variableVector gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticPartOfDeformation (test 24) & False\n";
            return 1;
        }

        //Test dChiedGradChi
        gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticPartOfDeformation (test 25) & False\n";
            return 1;
        }

        //Test dGradChiedGradChi
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dGradChiedGradChi[j][i] ) ){
                results << "test_computeElasticPartOfDeformation (test 26) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeElasticPartOfDeformation(  F,  chi,  gradChi,
                                                                                Fp, chip, gradChip - delta,
                                                                                resultFeM, resultChieM,
                                                                                resultGradChieM );

        if ( error ){
            error->print();
            results << "test_computeElasticPartOfDeformation & False\n";
            return 1;
        }

        //Test dFedGradChi
        variableVector gradCol = ( resultFeP - resultFeM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticPartOfDeformation (test 27) & False\n";
            return 1;
        }

        //Test dChiedGradChi
        gradCol = ( resultChieP - resultChieM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticPartOfDeformation (test 28) & False\n";
            return 1;
        }

        //Test dGradChiedGradChi
        gradCol = ( resultGradChieP - resultGradChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dGradChiedGradChip[j][i] ) ){
                results << "test_computeElasticPartOfDeformation (test 29) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeElasticPartOfDeformation & True\n";
    return 0;
}

int test_computeElasticDeformationMeasures( std::ofstream &results ){
    /*!
     * Test the computation of the elastic deformation measures 
     * required for the computation of the plasticity.
     *
     * :param std::ofstream &results: The output file
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

    if ( error ){
        error->print();
        results << "test_computeElasticDeformationMeasures & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultCe, answerCe ) ){
        results << "test_computeElasticDeformationMeasures (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultCChie, answerCChie ) ){
        results << "test_computeElasticDeformationMeasures (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultPsie, answerPsie ) ){
        results << "test_computeElasticDeformationMeasures (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultGammae, answerGammae ) ){
        results << "test_computeElasticDeformationMeasures (test 4) & False\n";
        return 1;
    }

    //Tests of Jacobians
    variableVector resultCeJ, resultCChieJ, resultPsieJ, resultGammaeJ;
    variableMatrix dRCGedFe, dMicroRCGedChie, dPsiedFe, dPsiedChie, dGammaedFe, dGammaedGradChie;

    error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie, gradChie, resultCeJ,
                                                                             resultCChieJ, resultPsieJ,
                                                                             resultGammaeJ, dRCGedFe,
                                                                             dMicroRCGedChie, dPsiedFe, dPsiedChie,
                                                                             dGammaedFe, dGammaedGradChie );

    if ( error ){
        error->print();
        results << "test_computeElasticDeformationMeasures & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultCeJ, answerCe ) ){
        results << "test_computeElasticDeformationMeasures (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultCChieJ, answerCChie ) ){
        results << "test_computeElasticDeformationMeasures (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultPsieJ, answerPsie ) ){
        results << "test_computeElasticDeformationMeasures (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultGammaeJ, answerGammae ) ){
        results << "test_computeElasticDeformationMeasures (test 8) & False\n";
        return 1;
    }

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
        if ( error ){
            error->print();
            results << "test_computeElasticDeformationMeasures & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe - delta, chie, gradChie,
                                                                                 resultCeM, resultCChieM, resultPsieM,
                                                                                 resultGammaeM );
        if ( error ){
            error->print();
            results << "test_computeElasticDeformationMeasures & False\n";
            return 1;
        }

        variableVector gradCol = ( resultCeP - resultCeM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dRCGedFe[j][i] ) ){
                results << "test_computeElasticDeformationMeasures (test 9) & False\n";
                return 1;
            }
        }

        gradCol = ( resultCChieP - resultCChieM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticDeformationMeasures (test 10) & False\n";
            return 1;
        }

        gradCol = ( resultPsieP - resultPsieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dPsiedFe[j][i] ) ){
                results << "test_computeElasticDeformationMeasures (test 11) & False\n";
                return 1;
            }
        }

        gradCol = ( resultGammaeP - resultGammaeM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dGammaedFe[j][i] ) ){
                results << "test_computeElasticDeformationMeasures (test 12) & False\n";
                return 1;
            }
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
        if ( error ){
            error->print();
            results << "test_computeElasticDeformationMeasures & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie - delta, gradChie,
                                                                                 resultCeM, resultCChieM, resultPsieM,
                                                                                 resultGammaeM );
        if ( error ){
            error->print();
            results << "test_computeElasticDeformationMeasures & False\n";
            return 1;
        }

        variableVector gradCol = ( resultCeP - resultCeM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticDeformationMeasures (test 13) & False\n";
            return 1;
        }

        gradCol = ( resultCChieP - resultCChieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dMicroRCGedChie[j][i] ) ){
                results << "test_computeElasticDeformationMeasures (test 14) & False\n";
                return 1;
            }
        }

        gradCol = ( resultPsieP - resultPsieM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dPsiedChie[j][i] ) ){
                results << "test_computeElasticDeformationMeasures (test 15) & False\n";
                return 1;
            }
        }

        gradCol = ( resultGammaeP - resultGammaeM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticDeformationMeasures (test 16) & False\n";
            return 1;
        }
    }

    for ( unsigned int i = 0; i < gradChie.size(); i++ ){
        constantVector delta( gradChie.size(), 0 );
        delta[i] = eps * fabs( gradChie[ i ] ) + eps;

        variableVector resultCeP, resultCChieP, resultPsieP, resultGammaeP;
        variableVector resultCeM, resultCChieM, resultPsieM, resultGammaeM;

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie, gradChie + delta,
                                                                                 resultCeP, resultCChieP, resultPsieP,
                                                                                 resultGammaeP );
        if ( error ){
            error->print();
            results << "test_computeElasticDeformationMeasures & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeElasticDeformationMeasures( Fe, chie, gradChie - delta,
                                                                                 resultCeM, resultCChieM, resultPsieM,
                                                                                 resultGammaeM );
        if ( error ){
            error->print();
            results << "test_computeElasticDeformationMeasures & False\n";
            return 1;
        }

        variableVector gradCol = ( resultCeP - resultCeM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticDeformationMeasures (test 17) & False\n";
            return 1;
        }

        gradCol = ( resultCChieP - resultCChieM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticDeformationMeasures (test 18) & False\n";
            return 1;
        }

        gradCol = ( resultPsieP - resultPsieM ) / ( 2 * delta[i] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computeElasticDeformationMeasures (test 19) & False\n";
            return 1;
        }

        gradCol = ( resultGammaeP - resultGammaeM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dGammaedGradChie[j][i] ) ){
                results << "test_computeElasticDeformationMeasures (test 20) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeElasticDeformationMeasures & True\n";
    return 0;
}

int test_computePlasticVelocityGradients( std::ofstream &results ){
    /*!
     * Test the computation of the plastic velocity gradients.
     *
     * :param std::ofstream &results: The output file.
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

    variableVector answerMicroLp = { -85.67983387, -16.91839826, 127.3318347 ,
                                       0.65035144,   0.1459189 ,  -0.71988301,
                                     -36.05794838,  -7.86041652,  52.33737079 };

    variableVector answerMicroGradientLp = { -83.50670143,  -11.14831106,  -39.09529065,  -16.58901227,
                                              -8.71869628,   19.76266969,   64.74989223,    3.55092183,
                                              59.30524632,  126.06867513,   26.68572242,   83.03661657,
                                              25.93560443,    6.13134257,   16.20263835, -185.220318  ,
                                             -38.17442863, -123.4797927 ,  -56.62958098,   -7.64862742,
                                             -16.00279814,  -11.26027965,   -4.19458173,    8.34519958,
                                              58.1455205 ,    5.42232797,   23.44933638 };

    variableVector resultLp, resultMicroLp, resultMicroGradientLp;

    errorOut error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                                    Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                                    microFlowDirection, microGradientFlowDirection,
                                                                                    resultLp, resultMicroLp, resultMicroGradientLp );

    if ( error ){
        error->print();
        results << "test_computePlasticVelocityGradients & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultLp, answerLp ) ){
        results << "test_computePlasticVelocityGradients (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroLp, answerMicroLp ) ){
        results << "test_computePlasticVelocityGradients (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroGradientLp, answerMicroGradientLp ) ){
        results << "test_computePlasticVelocityGradients (test 3) & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_computePlasticVelocityGradients & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultLpJ, answerLp ) ){
        results << "test_computePlasticVelocityGradients (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroLpJ, answerMicroLp ) ){
        results << "test_computePlasticVelocityGradients (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroGradientLpJ, answerMicroGradientLp ) ){
        results << "test_computePlasticVelocityGradients (test 6) & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_computePlasticVelocityGradients & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultLpJ2, answerLp ) ){
        results << "test_computePlasticVelocityGradients (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroLpJ2, answerMicroLp ) ){
        results << "test_computePlasticVelocityGradients (test 8) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroGradientLpJ2, answerMicroGradientLp ) ){
        results << "test_computePlasticVelocityGradients (test 9) & False\n";
        return 1;
    }

    //Tests of Jacobians w.r.t. macroGamma
    constantType eps = 1e-6;
    constantType scalarDelta = eps * fabs( macroGamma ) + eps;

    variableVector resultLpP, resultMicroLpP, resultMicroGradientLpP;
    variableVector resultLpM, resultMicroLpM, resultMicroGradientLpM;

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma + scalarDelta, microGamma, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpP, resultMicroLpP, resultMicroGradientLpP );

    if ( error ){
        error->print();
        results << "test_computePlasticVelocityGradients & False\n";
        return 1;
    }

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma - scalarDelta, microGamma, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpM, resultMicroLpM, resultMicroGradientLpM );

    if ( error ){
        error->print();
        results << "test_computePlasticVelocityGradients & False\n";
        return 1;
    }

    variableVector gradCol = ( resultLpP - resultLpM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradCol, dMacroLpdMacroGammaJ ) ){
        results << "test_computePlasticVelocityGradients (test 8) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( gradCol, dMacroLpdMacroGammaJ2 ) ){
        results << "test_computePlasticVelocityGradients (test 9) & False\n";
        return 1;
    }

    gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
        results << "test_computePlasticVelocityGradients (test 10) & False\n";
        return 1;
    }

    gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
        results << "test_computePlasticVelocityGradients (test 11) & False\n";
        return 1;
    }

    //Test of Jacobians w.r.t. microGamma
    scalarDelta = eps * fabs( microGamma ) + eps;

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma + scalarDelta, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpP, resultMicroLpP, resultMicroGradientLpP );

    if ( error ){
        error->print();
        results << "test_computePlasticVelocityGradients & False\n";
        return 1;
    }

    error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma - scalarDelta, microGradientGamma,
                                                                           Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                           microFlowDirection, microGradientFlowDirection,
                                                                           resultLpM, resultMicroLpM, resultMicroGradientLpM );

    if ( error ){
        error->print();
        results << "test_computePlasticVelocityGradients & False\n";
        return 1;
    }

    gradCol = ( resultLpP - resultLpM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradCol,  dMacroLpdMicroGammaJ ) ){
        results << "test_computePlasticVelocityGradients (test 13) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( gradCol,  dMacroLpdMicroGammaJ2 ) ){
        results << "test_computePlasticVelocityGradients (test 14) & False\n";
        return 1;
    }

    gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGammaJ ) ){
        results << "test_computePlasticVelocityGradients (test 14) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGammaJ2 ) ){
        results << "test_computePlasticVelocityGradients (test 15) & False\n";
        return 1;
    }

    gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradCol,  dMicroGradientLpdMicroGammaJ ) ){
        results << "test_computePlasticVelocityGradients (test 15) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( gradCol,  dMicroGradientLpdMicroGammaJ2 ) ){
        results << "test_computePlasticVelocityGradients (test 16) & False\n";
        return 1;
    }

    //Test of Jacobians w.r.t. the microGradientGamma
    for ( unsigned int i = 0; i < microGradientGamma.size(); i++ ){
        constantVector delta( microGradientGamma.size(), 0 );
        delta[i] = eps * fabs( microGradientGamma[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma + delta,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma - delta,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computePlasticVelocityGradients (test 16) & False\n";
            return 1;
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computePlasticVelocityGradients (test 17) & False\n";
            return 1;
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientLpdMicroGradientGammaJ[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 18) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientLpdMicroGradientGammaJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 19) & False\n";
                return 1;
            }
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
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce - delta, microCe, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMacroLdElasticRCGJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 20) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computePlasticVelocityGradients (test 21) & False\n";
            return 1;
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computePlasticVelocityGradients (test 22) & False\n";
            return 1;
        }
    }

    //Test of Jacobians w.r.t. the micro right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < microCe.size(); i++ ){
        constantVector delta( microCe.size(), 0 );
        delta[i] = eps * fabs( microCe[ i ] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe + delta, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpP, resultMicroLpP, resultMicroGradientLpP );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe - delta, Psie, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computePlasticVelocityGradients (test 23) & False\n";
            return 1;
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroLdMicroElasticRCGJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 24) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroElasticRCGJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 25) & False\n";
                return 1;
            }
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
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie - delta, Gammae, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computePlasticVelocityGradients (test 26) & False\n";
            return 1;
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroLdElasticPsiJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 27) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticPsiJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 28) & False\n";
                return 1;
            }
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
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae - delta, macroFlowDirection, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            results << "test_computePlasticVelocityGradients (test 29) & False\n";
            return 1;
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computePlasticVelocityGradients (test 30) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticGammaJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 31) & False\n";
                return 1;
            }
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
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection - delta, 
                                                                               microFlowDirection, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMacroLdMacroFlowDirectionJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 32) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computePlasticVelocityGradients (test 33) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computePlasticVelocityGradients (test 34) & False\n";
                return 1;
            }
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
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection,
                                                                               microFlowDirection - delta, microGradientFlowDirection,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMacroLdMicroFlowDirectionJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 35) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroLdMicroFlowDirectionJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 36) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroFlowDirectionJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 37) & False\n";
                return 1;
            }
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
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticVelocityGradients( macroGamma, microGamma, microGradientGamma,
                                                                               Ce, microCe, Psie, Gammae, macroFlowDirection,
                                                                               microFlowDirection, microGradientFlowDirection - delta,
                                                                               resultLpM, resultMicroLpM, resultMicroGradientLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticVelocityGradients & False\n";
            return 1;
        }

        gradCol = ( resultLpP - resultLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computePlasticVelocityGradients (test 38) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computePlasticVelocityGradients (test 39) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientLpP - resultMicroGradientLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientFlowDirectionJ2[ j ][ i ] ) ){
                results << "test_computePlasticVelocityGradients (test 40) & False\n";
                return 1;
            }
        }
    }

    results << "test_computePlasticVelocityGradients & True\n";
    return 0;
}

int test_computePlasticMacroVelocityGradient( std::ofstream &results ){
    /*!
     * Test the computation of the plastic macro velocity gradient.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_computePlasticMacroVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMacroLp, resultMacroLp ) ){
        results << "test_computePlasticMacroVelocityGradient (test 1) & False\n";
        return 1;
    }

    //Tests of the Jacobians
    variableVector resultMacroLpJ;
    variableVector dMacroLdMacroGammaJ, dMacroLdMicroGammaJ;

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                               macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpJ, dMacroLdMacroGammaJ,
                                                                               dMacroLdMicroGammaJ );

    if ( error ){
        error->print();
        results << "test_computePlasticMacroVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMacroLp, resultMacroLpJ ) ){
        results << "test_computePlasticMacroVelocityGradient (test 2) & False\n";
        return 1;
    }

    variableVector resultMacroLpJ2;
    variableVector dMacroLdMacroGammaJ2, dMacroLdMicroGammaJ2;
    variableMatrix dMacroLdElasticRCG, dMacroLdMacroFlowDirection, dMacroLdMicroFlowDirection;

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                               macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpJ2, dMacroLdMacroGammaJ2,
                                                                               dMacroLdMicroGammaJ2, dMacroLdElasticRCG,
                                                                               dMacroLdMacroFlowDirection,
                                                                               dMacroLdMicroFlowDirection );

    if ( error ){
        error->print();
        results << "test_computePlasticMacroVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMacroLp, resultMacroLpJ2 ) ){
        results << "test_computePlasticMacroVelocityGradient (test 3) & False\n";
        return 1;
    }

    //Tests Jacobians w.r.t. macroGamma
    constantType eps = 1e-6;
    constantType scalarDelta = eps * fabs( macroGamma) + eps;

    variableVector resultMacroLpP, resultMacroLpM;

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma + scalarDelta, microGamma,
                                                                               inverseCe, macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpP );

    if ( error ){
        error->print();
        results << "test_computePlasticMacroVelocityGradient & False\n";
        return 1;
    }

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma - scalarDelta, microGamma,
                                                                               inverseCe, macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpM );

    if ( error ){
        error->print();
        results << "test_computePlasticMacroVelocityGradient & False\n";
        return 1;
    }

    variableVector gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradCol, dMacroLdMacroGammaJ ) ){
        results << "test_computePlasticMacroVelocityGradient (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( gradCol, dMacroLdMacroGammaJ2 ) ){
        results << "test_computePlasticMacroVelocityGradient (test 5) & False\n";
        return 1;
    }

    //Test Jacobians w.r.t. microGamma
    scalarDelta = eps * fabs( microGamma) + eps;

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma + scalarDelta,
                                                                               inverseCe, macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpP );

    if ( error ){
        error->print();
        results << "test_computePlasticMacroVelocityGradient & False\n";
        return 1;
    }

    error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma - scalarDelta,
                                                                               inverseCe, macroFlowDirection, microFlowDirection,
                                                                               resultMacroLpM );

    if ( error ){
        error->print();
        results << "test_computePlasticMacroVelocityGradient & False\n";
        return 1;
    }

    gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradCol, dMacroLdMicroGammaJ ) ){
        results << "test_computePlasticMacroVelocityGradient (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( gradCol, dMacroLdMicroGammaJ2 ) ){
        results << "test_computePlasticMacroVelocityGradient (test 7) & False\n";
        return 1;
    }

    //Test Jacobians w.r.t. the right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < Ce.size(); i++ ){
        constantVector delta( Ce.size(), 0 );
        delta[i] = eps * fabs( Ce[i] ) + eps;

        variableVector inverseCeTemp = vectorTools::inverse( Ce + delta, 3, 3 );

        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCeTemp, macroFlowDirection,
                                                                                   microFlowDirection, resultMacroLpP );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMacroVelocityGradient & False\n";
            return 1;
        }

        inverseCeTemp = vectorTools::inverse( Ce - delta, 3, 3 );
    
        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCeTemp, macroFlowDirection,
                                                                                   microFlowDirection, resultMacroLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMacroVelocityGradient & False\n";
            return 1;
        }
    
        gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMacroLdElasticRCG[ j ][ i ] ) ){
                results << "test_computePlasticMacroVelocityGradient (test 8) & False\n";
                return 1;
            }
        }
    }

    //Test Jacobians w.r.t. the macro flow direction
    for ( unsigned int i = 0; i < macroFlowDirection.size(); i++ ){
        constantVector delta( macroFlowDirection.size(), 0 );
        delta[i] = eps * fabs( macroFlowDirection[i] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCe, macroFlowDirection + delta,
                                                                                   microFlowDirection, resultMacroLpP );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMacroVelocityGradient & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCe, macroFlowDirection - delta,
                                                                                   microFlowDirection, resultMacroLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMacroVelocityGradient & False\n";
            return 1;
        }
    
        gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMacroLdMacroFlowDirection[ j ][ i ] ) ){
                results << "test_computePlasticMacroVelocityGradient (test 9) & False\n";
                return 1;
            }
        }
    }

    //Test Jacobians w.r.t. the micro flow direction
    for ( unsigned int i = 0; i < microFlowDirection.size(); i++ ){
        constantVector delta( microFlowDirection.size(), 0 );
        delta[i] = eps * fabs( microFlowDirection[i] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCe, macroFlowDirection,
                                                                                   microFlowDirection + delta, resultMacroLpP );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMacroVelocityGradient & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                                   inverseCe, macroFlowDirection,
                                                                                   microFlowDirection - delta, resultMacroLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMacroVelocityGradient & False\n";
            return 1;
        }
    
        gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMacroLdMicroFlowDirection[ j ][ i ] ) ){
                results << "test_computePlasticMacroVelocityGradient (test 10) & False\n";
                return 1;
            }
        }
    }

    results << "test_computePlasticMacroVelocityGradient & True\n";
    return 0;
}


int test_computePlasticMicroVelocityGradient( std::ofstream &results ){
    /*!
     * Test the computation of the plastic micro velocity gradient.
     *
     * :param std::ofstream &results: The output file.
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

    variableVector answerMicroLp = { -85.67983387, -16.91839826, 127.3318347 ,
                                       0.65035144,   0.1459189 ,  -0.71988301,
                                     -36.05794838,  -7.86041652,  52.33737079 };

    variableVector resultMicroLp;

    errorOut error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                        microFlowDirection, resultMicroLp );

    if ( error ){
        error->print();
        results << "test_computePlasticMicroVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroLp, resultMicroLp ) ){
        results << "test_computePlasticMicroVelocityGradient (test 1) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableVector resultMicroLpJ;
    variableVector dMicroLpdMicroGamma;

    error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                               microFlowDirection, resultMicroLpJ,
                                                                               dMicroLpdMicroGamma );

    if ( error ){
        error->print();
        results << "test_computePlasticMicroVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroLp, resultMicroLpJ ) ){
        results << "test_computePlasticMicroVelocityGradient (test 2) & False\n";
        return 1;
    }

    variableVector resultMicroLpJ2;
    variableVector dMicroLpdMicroGamma2;
    variableMatrix dMicroLpdMicroRCG, dMicroLpdPsie, dMicroLpdMicroFlowDirection;

    error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                               microFlowDirection, resultMicroLpJ2,
                                                                               dMicroLpdMicroGamma2,
                                                                               dMicroLpdMicroRCG, dMicroLpdPsie,
                                                                               dMicroLpdMicroFlowDirection );

    if ( error ){
        error->print();
        results << "test_computePlasticMicroVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroLp, resultMicroLpJ ) ){
        results << "test_computePlasticMicroVelocityGradient (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroLp, resultMicroLpJ2 ) ){
        results << "test_computePlasticMicroVelocityGradient (test 3) & False\n";
        return 1;
    }

    constantType eps = 1e-6;
    constantType scalarDelta = eps * fabs( microGamma ) + eps;

    variableVector resultMicroLpP, resultMicroLpM;

    error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma + scalarDelta, Ce, Psie, invPsie,
                                                                               microFlowDirection, resultMicroLpP );

    if ( error ){
        error->print();
        results << "test_computePlasticMicroVelocityGradient & False\n";
        return 1;
    }

    error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma - scalarDelta, Ce, Psie, invPsie,
                                                                               microFlowDirection, resultMicroLpM );

    if ( error ){
        error->print();
        results << "test_computePlasticMicroVelocityGradient & False\n";
        return 1;
    }

    variableVector gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGamma ) ){
        results << "test_computePlasticMicroVelocityGradient (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGamma2 ) ){
        results << "test_computePlasticMicroVelocityGradient (test 5) & False\n";
        return 1;
    }

    //Test Jacobian w.r.t. the elastic micro right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < Ce.size(); i++ ){
        constantVector delta( Ce.size(), 0 );
        delta[i] = eps * fabs( Ce[i] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce + delta, Psie, invPsie,
                                                                                   microFlowDirection, resultMicroLpP );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMicroVelocityGradient & False\n";
            return 1;
        }
    
        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce - delta, Psie, invPsie,
                                                                                   microFlowDirection, resultMicroLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMicroVelocityGradient & False\n";
            return 1;
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dMicroLpdMicroRCG[ j ][ i ] ) ){
                results << "test_computePlasticMicroVelocityGradient (test 6) & False\n";
                return 1;
            }
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
    
        if ( error ){
            error->print();
            results << "test_computePlasticMicroVelocityGradient & False\n";
            return 1;
        }
    
        PsieTemp = Psie - delta;
        invPsieTemp = vectorTools::inverse( PsieTemp, 3, 3 );

        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, PsieTemp, invPsieTemp,
                                                                                   microFlowDirection, resultMicroLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMicroVelocityGradient & False\n";
            return 1;
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dMicroLpdPsie[ j ][ i ] ) ){
                results << "test_computePlasticMicroVelocityGradient (test 8) & False\n";
                return 1;
            }
        }
    }

    for ( unsigned int i = 0; i < microFlowDirection.size(); i++ ){
        constantVector delta( microFlowDirection.size(), 0 );
        delta[i] = eps * fabs( microFlowDirection[i] ) + eps;

        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                   microFlowDirection + delta, resultMicroLpP );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMicroVelocityGradient & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                   microFlowDirection - delta, resultMicroLpM );
    
        if ( error ){
            error->print();
            results << "test_computePlasticMicroVelocityGradient & False\n";
            return 1;
        }

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dMicroLpdMicroFlowDirection[ j ][ i ] ) ){
                results << "test_computePlasticMicroVelocityGradient (test 9) & False\n";
                return 1;
            }
        }
    }

    results << "test_computePlasticMicroVelocityGradient & True\n";
    return 0;
}

int test_computePlasticMicroGradientVelocityGradient( std::ofstream &results ){
    /*!
     * Test the computation of the plastic micro gradient velocity gradient.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_computePlasticMicroGradientVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLp ) ){
        results << "test_computePlasticMicroGradientVelocityGradient (test 1) & False\n";
        return 1;
    }

    variableVector resultMicroGradLp1;
    variableVector resultSkewTerm;

    error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie, 
                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                       microLp, resultMicroGradLp1,
                                                                                       resultSkewTerm );

    if ( error ){
        error->print();
        results << "test_computePlasticMicroGradientVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLp1 ) ){
        results << "test_computePlasticMicroGradientVelocityGradient (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerSkewTerm, resultSkewTerm ) ){
        results << "test_computePlasticMicroGradientVelocityGradient (test 3) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableVector resultMicroGradLpJ;
    variableMatrix dPlasticMicroGradientLdMicroGradientGamma;
    variableMatrix dPlasticMicroGradientLdPlasticMicroL;

    error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                       microLp, resultMicroGradLpJ,
                                                                                       dPlasticMicroGradientLdMicroGradientGamma,
                                                                                       dPlasticMicroGradientLdPlasticMicroL );
    if ( error ){
        error->print();
        results << "test_computePlasticMicroGradientVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLpJ ) ){
        results << "test_computePlasticMicroGradientVelocityGradient (test 4) & False\n";
        return 1;
    }

    variableVector resultMicroGradLpJ2, resultSkewTerm2;
    variableMatrix dPlasticMicroGradientLdMicroGradientGamma2;
    variableMatrix dPlasticMicroGradientLdPlasticMicroL2;

    error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                       microLp, resultMicroGradLpJ2, resultSkewTerm2,
                                                                                       dPlasticMicroGradientLdMicroGradientGamma2,
                                                                                       dPlasticMicroGradientLdPlasticMicroL2 );
    if ( error ){
        error->print();
        results << "test_computePlasticMicroGradientVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLpJ2 ) ){
        results << "test_computePlasticMicroGradientVelocityGradient (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerSkewTerm, resultSkewTerm2 ) ){
        results << "test_computePlasticMicroGradientVelocityGradient (test 6) & False\n";
        return 1;
    }

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
    if ( error ){
        error->print();
        results << "test_computePlasticMicroGradientVelocityGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLpJ3 ) ){
        results << "test_computePlasticMicroGradientVelocityGradient (test 5) & False\n";
        return 1;
    }

    //Test computation of Jacobians w.r.t. microGradientGamma
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < microGradientGamma.size(); i++ ){
        constantVector delta( microGradientGamma.size(), 0 );
        delta[i] = eps * fabs( microGradientGamma[i] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma + delta, Psie, invPsie,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpP );

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma - delta, Psie, invPsie,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpM );

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientGamma[ j ][ i ] ) ){
                results << "test_computePlasticMicroGradientVelocityGradient (test 7) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientGamma2[ j ][ i ] ) ){
                results << "test_computePlasticMicroGradientVelocityGradient (test 8) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientGamma3[ j ][ i ] ) ){
                results << "test_computePlasticMicroGradientVelocityGradient (test 9) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp - delta, resultMicroGradLpM );

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdPlasticMicroL[ j ][ i ] ) ){
                results << "test_computePlasticMicroGradientVelocityGradient (test 10) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdPlasticMicroL2[ j ][ i ] ) ){
                results << "test_computePlasticMicroGradientVelocityGradient (test 11) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdPlasticMicroL3[ j ][ i ] ) ){
                results << "test_computePlasticMicroGradientVelocityGradient (test 12) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        PsieTemp = Psie - delta;
        invPsieTemp = vectorTools::inverse( PsieTemp, 3, 3 );

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, PsieTemp, invPsieTemp,
                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpM );

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticPsi[ j ][ i ] ) ){
                results << "test_computePlasticMicroGradientVelocityGradient (test 13) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                           elasticGamma - delta,
                                                                                           microGradientFlowDirection,
                                                                                           microLp, resultMicroGradLpM );

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticGamma[ j ][ i ] ) ){
                results << "test_computePlasticMicroGradientVelocityGradient (test 14) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                           elasticGamma,
                                                                                           microGradientFlowDirection - delta,
                                                                                           microLp, resultMicroGradLpM );

        if ( error ){
            error->print();
            results << "test_computePlasticMicroGradientVelocityGradient & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientFlowDirection[ j ][ i ] ) ){
                results << "test_computePlasticMicroGradientVelocityGradient (test 15) & False\n";
                return 1;
            }
        }
    }

    results << "test_computePlasticMicroGradientVelocityGradient & True\n";
    return 0;
}

int test_evolvePlasticMicroGradChi( std::ofstream &results ){
    /*!
     * Test the evolution of the plastic micro gradient deformation.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_evolvePlasticMicroGradChi & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerCurrentPlasticMicroGradient, resultCurrentPlasticMicroGradient ) ){
        results << "test_evolvePlasticMicroGradChi (test 1)  & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_evolvePlasticMicroGradChi & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerCurrentPlasticMicroGradient, resultCurrentPlasticMicroGradient2 ) ){
        results << "test_evolvePlasticMicroGradChi (test 2)  & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerLHS, LHS2 ) ){
        results << "test_evolvePlasticMicroGradChi (test 3)  & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_evolvePlasticMicroGradChi & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerCurrentPlasticMicroGradient, resultCurrentPlasticMicroGradientJ ) ){
        results << "test_evolvePlasticMicroGradChi (test 4)  & False\n";
        return 1;
    }

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

        if ( error ){
            error->print();
            results << "test_evolvePlasticMicroGradChi & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evolvePlasticMicroGradChi & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMicroDeformation[ j ][ i ] ) ){
                results << "tet_evolvePlasticMicroGradChi (test 5) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_evolvePlasticMicroGradChi & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evolvePlasticMicroGradChi & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient[ j ][ i ] ) ){
                results << "tet_evolvePlasticMicroGradChi (test 6) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_evolvePlasticMicroGradChi & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evolvePlasticMicroGradChi & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient[ j ][ i ] ) ){
                results << "tet_evolvePlasticMicroGradChi (test 7) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_evolvePlasticMicroGradChi & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evolvePlasticMicroGradChi & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient[ j ][ i ] ) ){
                results << "tet_evolvePlasticMicroGradChi (test 8) & False\n";
                return 1;
            }
        }
    }

    results << "test_evolvePlasticMicroGradChi & True\n";
    return 0;
}

int test_evolvePlasticDeformation( std::ofstream &results ){
    /*!
     * Evolve the plastic deformation.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        results << "test_evolvePlasticDeformation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMacro, answerMacro ) ){
        results << "test_evolvePlasticDeformation (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicro, answerMicro ) ){
        results << "test_evolvePlasticDeformation (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroGrad, answerMicroGrad ) ){
        results << "test_evolvePlasticDeformation (test 3) & False\n";
        return 1;
    }

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

    if ( error ){
        results << "test_evolvePlasticDeformation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMacroJ, answerMacro ) ){
        results << "test_evolvePlasticDeformation (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroJ, answerMicro ) ){
        results << "test_evolvePlasticDeformation (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroGradJ, answerMicroGrad ) ){
        results << "test_evolvePlasticDeformation (test 6) & False\n";
        return 1;
    }

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

        if ( error ){
            results << "test_evolvePlasticDeformation & False\n";
            return 1;
        }

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

        if ( error ){
            results << "test_evolvePlasticDeformation & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dFdMacroL[ j ][ i ] ) ){
                results << "test_evolvePlasticDeformation (test 7) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_evolvePlasticDeformation (test 8) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dGradChidMacroL[ j ][ i ] ) ){
                results << "test_evolvePlasticDeformation (test 9) & False\n";
                return 1;
            }
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

        if ( error ){
            results << "test_evolvePlasticDeformation & False\n";
            return 1;
        }

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

        if ( error ){
            results << "test_evolvePlasticDeformation & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_evolvePlasticDeformation (test 10) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dChidMicroL[ j ][ i ] ) ){
                results << "test_evolvePlasticDeformation (test 11) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dGradChidMicroL[ j ][ i ] ) ){
                results << "test_evolvePlasticDeformation (test 12) & False\n";
                return 1;
            }
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

        if ( error ){
            results << "test_evolvePlasticDeformation & False\n";
            return 1;
        }

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

        if ( error ){
            results << "test_evolvePlasticDeformation & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_evolvePlasticDeformation (test 13) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_evolvePlasticDeformation (test 14) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dGradChidMicroGradL[ j ][ i ] ) ){
                results << "test_evolvePlasticDeformation (test 15) & False\n";
                return 1;
            }
        }
    }

    results << "test_evolvePlasticDeformation & True\n";
    return 0;
}

int test_evolveStrainStateVariables( std::ofstream &results){
    /*!
     * Test the evolution of the strain-like state variables.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMacroISV, answerMacroISV ) ){
        results << "test_evolveStrainStateVariables (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroISV, answerMicroISV ) ){
        results << "test_evolveStrainStateVariables (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroGradISV, answerMicroGradISV ) ){
        results << "test_evolveStrainStateVariables (test 3) & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMacroISVJ, answerMacroISV ) ){
        results << "test_evolveStrainStateVariables (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroISVJ, answerMicroISV ) ){
        results << "test_evolveStrainStateVariables (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroGradISVJ, answerMicroGradISV ) ){
        results << "test_evolveStrainStateVariables (test 6) & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

    variableType gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradS, dMacroISVdMacroGamma ) ){
        results << "test_evolveStrainStateVariables (test 7) & False\n";
        return 1;
    }

    gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradS, 0. ) ){
        results << "test_evolveStrainStateVariables (test 8) & False\n";
        return 1;
    }

    variableVector gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * scalarDelta );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_evolveStrainStateVariables (test 9) * False\n";
            return 1;
        }
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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

    gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradS, 0. ) ){
        results << "test_evolveStrainStateVariables (test 10) & False\n";
        return 1;
    }

    gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradS, dMicroISVdMicroGamma ) ){
        results << "test_evolveStrainStateVariables (test 11) & False\n";
        return 1;
    }

    gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * scalarDelta );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_evolveStrainStateVariables (test 12) * False\n";
            return 1;
        }
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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

    gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradS, dMacroISVddGdMacroC ) ){
        results << "test_evolveStrainStateVariables (test 13) & False\n";
        return 1;
    }

    gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradS, 0. ) ){
        results << "test_evolveStrainStateVariables (test 14) & False\n";
        return 1;
    }

    gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * scalarDelta );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_evolveStrainStateVariables (test 15) * False\n";
            return 1;
        }
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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_evolveStrainStateVariables & False\n";
        return 1;
    }

    gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradS, 0. ) ){
        results << "test_evolveStrainStateVariables (test 16) & False\n";
        return 1;
    }

    gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * scalarDelta );

    if ( !vectorTools::fuzzyEquals( gradS, dMicroISVddGdMicroC ) ){
        results << "test_evolveStrainStateVariables (test 17) & False\n";
        return 1;
    }

    gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * scalarDelta );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_evolveStrainStateVariables (test 18) * False\n";
            return 1;
        }
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
    
        if ( error ){
            error->print();
            results << "test_evolveStrainStateVariables & False\n";
            return 1;
        }
    
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
    
        if ( error ){
            error->print();
            results << "test_evolveStrainStateVariables & False\n";
            return 1;
        }
    
        gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * delta[ i ] );
    
        if ( !vectorTools::fuzzyEquals( gradS, 0. ) ){
            results << "test_evolveStrainStateVariables (test 19) & False\n";
            return 1;
        }
    
        gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * delta[ i ] );
    
        if ( !vectorTools::fuzzyEquals( gradS, 0. ) ){
            results << "test_evolveStrainStateVariables (test 20) & False\n";
            return 1;
        }
    
        gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * delta[ i ] );
    
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientISVdMicroGradGamma[ j ][ i ] ) ){
                results << "test_evolveStrainStateVariables (test 21) * False\n";
                return 1;
            }
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
    
        if ( error ){
            error->print();
            results << "test_evolveStrainStateVariables & False\n";
            return 1;
        }
    
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
    
        if ( error ){
            error->print();
            results << "test_evolveStrainStateVariables & False\n";
            return 1;
        }
    
        gradS = ( resultMacroISVP - resultMacroISVM ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );
    
        if ( !vectorTools::fuzzyEquals( gradS, 0. ) ){
            results << "test_evolveStrainStateVariables (test 22) & False\n";
            return 1;
        }
    
        gradS = ( resultMicroISVP - resultMicroISVM ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );
    
        if ( !vectorTools::fuzzyEquals( gradS, 0. ) ){
            results << "test_evolveStrainStateVariables (test 23) & False\n";
            return 1;
        }
    
        gradCol = ( resultMicroGradISVP - resultMicroGradISVM ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );
    
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientISVddGdMicroGradGamma[ j ][ i ] ) ){
                results << "test_evolveStrainStateVariables (test 24) * False\n";
                return 1;
            }
        }
    }

    results << "test_evolveStrainStateVariables & True\n";
    return 0;
}

int test_computeFlowDirections( std::ofstream &results ){
    /*!
     * Test the computation of the flow directions
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_computeFlowDirections & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMacroFlowDirection, resultMacroFlowDirection ) ){
        results << "test_computeFlowDirections (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroFlowDirection, resultMicroFlowDirection ) ){
        results << "test_computeFlowDirections (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroGradientFlowDirection, resultMicroGradientFlowDirection ) ){
        results << "test_computeFlowDirections (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerdGdMacroCohesion, resultdGdMacroCohesion ) ){
        results << "test_computeFlowDirection (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerdGdMicroCohesion, resultdGdMicroCohesion ) ){
        results << "test_computeFlowDirection (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerdGdMicroGradientCohesion, resultdGdMicroGradientCohesion ) ){
        results << "test_computeFlowDirection (test 6) & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_computeFlowDirections & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMacroFlowDirection, resultMacroFlowDirectionJ ) ){
        results << "test_computeFlowDirections (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroFlowDirection, resultMicroFlowDirectionJ ) ){
        results << "test_computeFlowDirections (test 8) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroGradientFlowDirection, resultMicroGradientFlowDirectionJ ) ){
        results << "test_computeFlowDirections (test 9) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerdGdMacroCohesion, resultdGdMacroCohesionJ ) ){
        results << "test_computeFlowDirection (test 10) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerdGdMicroCohesion, resultdGdMicroCohesionJ ) ){
        results << "test_computeFlowDirection (test 11) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerdGdMicroGradientCohesion, resultdGdMicroGradientCohesionJ ) ){
        results << "test_computeFlowDirection (test 12) & False\n";
        return 1;
    }

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

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress - delta, referenceMicroStress, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMacroFlowDirectiondPK2Stress[ j ][ i ] ) ){
                std::cout << "i, j: " << i << ", " << j << "\n";
                std::cout << "result:\n"; vectorTools::print( dMacroFlowDirectiondPK2Stress );
                std::cout << "answer:\n"; vectorTools::print( gradCol );
                results << "test_computeFlowDirections (test 13) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeFlowDirections (test 14) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeFlowDirections (test 15) & False\n";
                return 1;
            }
        }

        variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 16) & False\n";
            return 1;
        }

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 17) & False\n";
            return 1;
        }

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                if ( !vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) ){
                    results << "test_computeFlowDirections (test 18) & False\n";
                    return 1;
                }
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

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress - delta, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeFlowDirections (test 19) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroFlowDirectiondSigma[ j ][ i ] ) ){
                std::cout << "i, j: " << i << ", " << j << "\n";
                std::cout << "result:\n"; vectorTools::print( dMicroFlowDirectiondSigma );
                std::cout << "answer:\n"; vectorTools::print( gradCol );
                results << "test_computeFlowDirections (test 20) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeFlowDirections (test 21) & False\n";
                return 1;
            }
        }

        variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 22) & False\n";
            return 1;
        }

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 23) & False\n";
            return 1;
        }

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                if ( !vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) ){
                    results << "test_computeFlowDirections (test 24) & False\n";
                    return 1;
                }
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

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress - delta,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeFlowDirections (test 25) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeFlowDirections (test 26) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFlowDirectiondM[ j ][ i ] ) ){
                results << "test_computeFlowDirections (test 27) & False\n";
                return 1;
            }
        }

        variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 28) & False\n";
            return 1;
        }

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 29) & False\n";
            return 1;
        }

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                if ( !vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) ){
                    results << "test_computeFlowDirections (test 30) & False\n";
                    return 1;
                }
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

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion,
                                                                     elasticRightCauchyGreen - delta, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMacroFlowDirectiondElasticRCG[ j ][ i ] ) ){
                results << "test_computeFlowDirections (test 31) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroFlowDirectiondElasticRCG[ j ][ i ] ) ){
                results << "test_computeFlowDirections (test 32) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFlowDirectiondElasticRCG[ j ][ i ] ) ){
                results << "test_computeFlowDirections (test 33) & False\n";
                return 1;
            }
        }

        variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 34) & False\n";
            return 1;
        }

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 35) & False\n";
            return 1;
        }

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                if ( !vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) ){
                    results << "test_computeFlowDirections (test 36) & False\n";
                    return 1;
                }
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

    if ( error ){
        error->print();
        results << "test_computeFlowDirections & False\n";
        return 1;
    }

    error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                 macroCohesion - deltaScalar, microCohesion, microGradientCohesion,
                                                                 elasticRightCauchyGreen, macroFlowParameters,
                                                                 microFlowParameters, microGradientFlowParameters,
                                                                 resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                 resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                 resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

    if ( error ){
        error->print();
        results << "test_computeFlowDirections & False\n";
        return 1;
    }

    variableVector gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_computeFlowDirections (test 37) & False\n";
            return 1;
        }
    }

    gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_computeFlowDirections (test 38) & False\n";
            return 1;
        }
    }

    gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_computeFlowDirections (test 39) & False\n";
            return 1;
        }
    }

    variableType gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * deltaScalar );

    if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
        results << "test_computeFlowDirections (test 40) & False\n";
        return 1;
    }

    gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * deltaScalar );

    if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
        results << "test_computeFlowDirections (test 41) & False\n";
        return 1;
    }

    variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradMat.size(); j++ ){
        for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
            if ( !vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) ){
                results << "test_computeFlowDirections (test 42) & False\n";
                return 1;
            }
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

    if ( error ){
        error->print();
        results << "test_computeFlowDirections & False\n";
        return 1;
    }

    error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                 macroCohesion, microCohesion - deltaScalar, microGradientCohesion,
                                                                 elasticRightCauchyGreen, macroFlowParameters,
                                                                 microFlowParameters, microGradientFlowParameters,
                                                                 resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                 resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                 resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

    if ( error ){
        error->print();
        results << "test_computeFlowDirections & False\n";
        return 1;
    }

    gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_computeFlowDirections (test 43) & False\n";
            return 1;
        }
    }

    gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_computeFlowDirections (test 44) & False\n";
            return 1;
        }
    }

    gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradCol.size(); j++ ){
        if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
            results << "test_computeFlowDirections (test 45) & False\n";
            return 1;
        }
    }

    gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * deltaScalar );

    if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
        results << "test_computeFlowDirections (test 46) & False\n";
        return 1;
    }

    gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * deltaScalar );

    if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
        results << "test_computeFlowDirections (test 47) & False\n";
        return 1;
    }

    gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * deltaScalar );

    for ( unsigned int j = 0; j < gradMat.size(); j++ ){
        for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
            if ( !vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) ){
                results << "test_computeFlowDirections (test 48) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeFlowDirections( PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                     macroCohesion, microCohesion, microGradientCohesion - delta,
                                                                     elasticRightCauchyGreen, macroFlowParameters,
                                                                     microFlowParameters, microGradientFlowParameters,
                                                                     resultMacroFlowDirectionM, resultMicroFlowDirectionM,
                                                                     resultMicroGradientFlowDirectionM, resultdGdMacroCohesionM,
                                                                     resultdGdMicroCohesionM, resultdGdMicroGradientCohesionM );

        if ( error ){
            error->print();
            results << "test_computeFlowDirections & False\n";
            return 1;
        }

        gradCol = ( resultMacroFlowDirectionP - resultMacroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeFlowDirections (test 49) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroFlowDirectionP - resultMicroFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeFlowDirections (test 50) & False\n";
                return 1;
            }
        }

        gradCol = ( resultMicroGradientFlowDirectionP - resultMicroGradientFlowDirectionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeFlowDirections (test 51) & False\n";
                return 1;
            }
        }

        gradScalar = ( resultdGdMacroCohesionP - resultdGdMacroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 52) & False\n";
            return 1;
        }

        gradScalar = ( resultdGdMicroCohesionP - resultdGdMicroCohesionM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradScalar, 0. ) ){
            results << "test_computeFlowDirections (test 53) & False\n";
            return 1;
        }

        variableMatrix gradMat = ( resultdGdMicroGradientCohesionP - resultdGdMicroGradientCohesionM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradMat.size(); j++ ){
            for ( unsigned int k = 0; k < gradMat[ j ].size(); k++ ){
                if ( !vectorTools::fuzzyEquals( gradMat[ j ][ k ], 0. ) ){
                    results << "test_computeFlowDirections (test 54) & False\n";
                    return 1;
                }
            }
        }
    }

    results << "test_computeFlowDirections & True\n";
    return 0;
}

int test_extractMaterialParameters( std::ofstream &results ){
    /*!
     * Test the extraction of the material parameters.
     *
     * :param std::ofstream &results: The output file.
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
    if ( error ){
        error->print();
        results << "test_extractMaterialParameters & False\n";
        return 1;
    }

    error = micromorphicLinearElasticity::formIsotropicB( 2.8, 0.76, 0.15, 9.8, 5.4, answerBmatrix );
    if ( error ){
        error->print();
        results << "test_extractMaterialParameters & False\n";
        return 1;
    }

    error = micromorphicLinearElasticity::formIsotropicC( { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.}, answerCmatrix );
    if ( error ){
        error->print();
        results << "test_extractMaterialParameters & False\n";
        return 1;
    }

    error = micromorphicLinearElasticity::formIsotropicD( 0.76, 5.4, answerDmatrix );
    if ( error ){
        error->print();
        results << "test_extractMaterialParameters & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_extractMaterialParameters & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( macroHardeningParameters, answerMacroHardeningParameters ) ){
        results << "test_extractMaterialParameters (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microHardeningParameters, answerMicroHardeningParameters ) ){
        results << "test_extractMaterialParameters (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microGradientHardeningParameters, answerMicroGradientHardeningParameters ) ){
        results << "test_extractMaterialParameters (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( macroFlowParameters, answerMacroFlowParameters ) ){
        results << "test_extractMaterialParameters (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microFlowParameters, answerMicroFlowParameters ) ){
        results << "test_extractMaterialParameters (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microGradientFlowParameters, answerMicroGradientFlowParameters ) ){
        results << "test_extractMaterialParameters (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( macroYieldParameters, answerMacroYieldParameters ) ){
        results << "test_extractMaterialParameters (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microYieldParameters, answerMicroYieldParameters ) ){
        results << "test_extractMaterialParameters (test 8) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microGradientYieldParameters, answerMicroGradientYieldParameters ) ){
        results << "test_extractMaterialParameters (test 9) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( Amatrix, answerAmatrix ) ){
        results << "test_extractMaterialParameters (test 10) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( Bmatrix, answerBmatrix ) ){
        results << "test_extractMaterialParameters (test 11) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( Cmatrix, answerCmatrix ) ){
        results << "test_extractMaterialParameters (test 12) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( Dmatrix, answerDmatrix ) ){
        results << "test_extractMaterialParameters (test 13) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( alphaMacro, answerAlphaMacro ) ){
        results << "test_extractMaterialParameters (test 14) & False\n";
    }

    if ( !vectorTools::fuzzyEquals( alphaMicro, answerAlphaMicro ) ){
        results << "test_extractMaterialParameters (test 15) & False\n";
    }

    if ( !vectorTools::fuzzyEquals( alphaMicroGradient, answerAlphaMicroGradient ) ){
        results << "test_extractMaterialParameters (test 16) & False\n";
    }

    if ( !vectorTools::fuzzyEquals( alphaMicro, answerAlphaMicro ) ){
        results << "test_extractMaterialParameters (test 16) & False\n";
    }

    if ( !vectorTools::fuzzyEquals( relativeTolerance, answerRelativeTolerance ) ){
        results << "test_extractMaterialParameters (test 17) & False\n";
    }

    if ( !vectorTools::fuzzyEquals( absoluteTolerance, answerAbsoluteTolerance ) ){
        results << "test_extractMaterialParameters (test 18) & False\n";
    }

    results << "test_extractMaterialParameters & True\n";
    return 0;
}

int test_extractStateVariables( std::ofstream &results ){
    /*!
     * Test the extraction of the state variables.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_extractStateVariables & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerPreviousMacroStrainISV, previousMacroStrainISV ) ){
        results << "test_extractStateVariables (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerPreviousMicroStrainISV, previousMicroStrainISV ) ){
        results << "test_extractStateVariables (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerPreviousMicroGradientStrainISV, previousMicroGradientStrainISV ) ){
        results << "test_extractStateVariables (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerPreviousMacroGamma, previousMacroGamma ) ){
        results << "test_extractStateVariables (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerPreviousMicroGamma, previousMicroGamma ) ){
        results << "test_extractStateVariables (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerPreviousMicroGradientGamma, previousMicroGradientGamma ) ){
        results << "test_extractStateVariables (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerPreviousPlasticDeformationGradient, previousPlasticDeformationGradient ) ){
        vectorTools::print( answerPreviousPlasticDeformationGradient );
        vectorTools::print( previousPlasticDeformationGradient );
        results << "test_extractStateVariables (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerPreviousPlasticMicroDeformation, previousPlasticMicroDeformation ) ){
        results << "test_extractStateVariables (test 8) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerPreviousPlasticGradientMicroDeformation, previousPlasticGradientMicroDeformation ) ){
        results << "test_extractStateVariables (test 9) & False\n";
        return 1;
    }

    results << "test_extractStateVariables & True\n";
    return 0;
}

int test_cout_redirect( std::ofstream &results ){
    /*!
     * Test the utility function which redirects cout to a string buffer.
     *
     * :param std::ofstream &results: The output file.
     */

    std::stringbuf buffer;
    micromorphicElastoPlasticity::cout_redirect rd( &buffer );
    
    std::string answer = "hello world\n";

    std::cout << answer;

    if ( answer.compare( buffer.str() ) != 0 ){
        results << "test_cout_redirect (test 1) & False\n";
        return 1;
    }

    results << "test_cout_redirect & True\n";
    return 0;
}

int test_cerr_redirect( std::ofstream &results ){
    /*!
     * Test the utility function which redirects cerr to a string buffer.
     *
     * :param std::ofstream &results: The output file.
     */

    std::stringbuf buffer;
    micromorphicElastoPlasticity::cerr_redirect rd( &buffer );
    
    std::string answer = "hello world\n";

    std::cerr << answer;

    if ( answer.compare( buffer.str() ) != 0 ){
        results << "test_cerr_redirect (test 1) & False\n";
        return 1;
    }

    results << "test_cerr_redirect & True\n";
    return 0;
}

int test_assembleFundamentalDeformationMeasures( std::ofstream &results ){
    /*!
     * Assemble the fundamental deformation measures from the degrees of freedom.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        results << "test_assembleFundamentalDeformationMeasures & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultF, answerDeformationGradient ) ){
        results << "test_assembleFundamentalDeformationMeasures (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultChi, answerMicroDeformation ) ){
        results << "test_assembleFundamentalDeformationMeasures (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultGradChi, answerGradientMicroDeformation ) ){
        results << "test_assembleFundamentalDeformationMeasures (test 3) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableVector resultFJ, resultChiJ, resultGradChiJ;
    variableMatrix dFdGradU, dChidPhi, dGradChidGradPhi;

    error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, grad_phi,
                                                                                  resultFJ, resultChiJ, resultGradChiJ,
                                                                                  dFdGradU, dChidPhi, dGradChidGradPhi );

    if ( error ){
        results << "test_assembleFundamentalDeformationMeasures & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultFJ, answerDeformationGradient ) ){
        results << "test_assembleFundamentalDeformationMeasures (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultChiJ, answerMicroDeformation ) ){
        results << "test_assembleFundamentalDeformationMeasures (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultGradChiJ, answerGradientMicroDeformation ) ){
        results << "test_assembleFundamentalDeformationMeasures (test 6) & False\n";
        return 1;
    }

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

        if ( error ){
            results << "test_assembleFundamentalDeformationMeasures & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( negative_perturb, phi, grad_phi,
                                                                                      FM, chiM, gradChiM );

        if ( error ){
            results << "test_assembleFundamentalDeformationMeasures & False\n";
            return 1;
        }

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dFdGradU[ j ][ i ] ) ){
                results << "test_assembleFundamentalDeformationMeasures (test 7) & False\n";
                return 1;
            }
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_assembleFundamentalDeformationMeasures (test 8) & False\n";
            }
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_assembleFundamentalDeformationMeasures (test 9) & False\n";
            }
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

        if ( error ){
            results << "test_assembleFundamentalDeformationMeasures & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, negative_perturb, grad_phi,
                                                                                      FM, chiM, gradChiM );

        if ( error ){
            results << "test_assembleFundamentalDeformationMeasures & False\n";
            return 1;
        }

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_assembleFundamentalDeformationMeasures (test 10) & False\n";
                return 1;
            }
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dChidPhi[ j ][ i ] ) ){
                results << "test_assembleFundamentalDeformationMeasures (test 11) & False\n";
            }
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_assembleFundamentalDeformationMeasures (test 12) & False\n";
            }
        }
    }

    for ( unsigned int i = 0; i < 27; i++ ){
        constantMatrix delta( 9, constantVector( 3, 0 ) );
        unsigned int ii, ij;
        ii = ( int )( i / 3 );
        ij = i % 3;
        delta[ ii ][ ij ] = eps * fabs( grad_u[ ii ][ ij ] ) + eps;

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

        if ( error ){
            results << "test_assembleFundamentalDeformationMeasures & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, negative_perturb,
                                                                                      FM, chiM, gradChiM );

        if ( error ){
            results << "test_assembleFundamentalDeformationMeasures & False\n";
            return 1;
        }

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_assembleFundamentalDeformationMeasures (test 13) & False\n";
                return 1;
            }
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_assembleFundamentalDeformationMeasures (test 14) & False\n";
            }
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dGradChidGradPhi[ j ][ i ] ) ){
                results << "test_assembleFundamentalDeformationMeasures (test 15) & False\n";
            }
        }
    }

    results << "test_assembleFundamentalDeformationMeasures & True\n";
    return 0;
}

int test_evaluateYieldFunctions( std::ofstream &results ){
    /*!
     * Test the evaluation of the yield functions
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_evaluateYieldFunctions & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_evaluateYieldFunctions (test 1) & False\n";
        return 1;
    }

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

    if ( error ){
        error->print();
        results << "test_evaluateYieldFunctions & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_evaluateYieldFunctions (test 2) & False\n";
        return 1;
    }

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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol[ 0 ], dMacroFdPK2[ i ] ) ){
            results << "test_evaluateYieldFunctions (test 3) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol[ 1 ], 0. ) ){
            results << "test_evaluateYieldFunctions (test 4) & False\n";
            return 1;
        }

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_evaluateYieldFunctions (test 5) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol[ 0 ], 0. ) ){
            results << "test_evaluateYieldFunctions (test 6) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol[ 1 ], dMicroFdSigma[ i ] ) ){
            results << "test_evaluateYieldFunctions (test 7) & False\n";
            return 1;
        }

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_evaluateYieldFunctions (test 8) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol[ 0 ], 0. ) ){
            results << "test_evaluateYieldFunctions (test 9) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol[ 1 ], 0. ) ){
            results << "test_evaluateYieldFunctions (test 10) & False\n";
            return 1;
        }

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFdM[ j - 2 ][ i ] ) ){
                results << "test_evaluateYieldFunctions (test 11) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta );

        if ( !vectorTools::fuzzyEquals( gradCol[ 0 ], dMacroFdMacroC ) ){
            results << "test_evaluateYieldFunctions (test 12) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol[ 1 ], 0. ) ){
            results << "test_evaluateYieldFunctions (test 13) & False\n";
            return 1;
        }

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_evaluateYieldFunctions (test 14) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta );

        if ( !vectorTools::fuzzyEquals( gradCol[ 0 ], 0. ) ){
            results << "test_evaluateYieldFunctions (test 15) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol[ 1 ], dMicroFdMicroC ) ){
            results << "test_evaluateYieldFunctions (test 16) & False\n";
            return 1;
        }

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_evaluateYieldFunctions (test 17) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol[ 0 ], 0. ) ){
            results << "test_evaluateYieldFunctions (test 18) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol[ 1 ], 0. ) ){
            results << "test_evaluateYieldFunctions (test 19) & False\n";
            return 1;
        }

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFdMicroGradientC[ j - 2 ][ i ] ) ){
                results << "test_evaluateYieldFunctions (test 20) & False\n";
                return 1;
            }
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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

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

        if ( error ){
            error->print();
            results << "test_evaluateYieldFunctions & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( gradCol[ 0 ], dMacroFdElasticRCG[ i ] ) ){
            results << "test_evaluateYieldFunctions (test 21) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol[ 1 ], dMicroFdElasticRCG[ i ] ) ){
            results << "test_evaluateYieldFunctions (test 22) & False\n";
            return 1;
        }

        for ( unsigned int j = 2; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientFdElasticRCG[ j - 2 ][ i ] ) ){
                results << "test_evaluateYieldFunctions (test 23) & False\n";
                return 1;
            }
        }
    }

    results << "test_evaluateYieldFunctions & True\n";
    return 0;
}

int test_computeCohesion( std::ofstream &results ){
    /*!
     * Test the computation of the cohesion
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        results << "test_computeCohesion & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( macroC, answerMacroC ) ){
        results << "test_computeCohesion (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microC, answerMicroC ) ){
        results << "test_computeCohesion (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microGradientC, answerMicroGradientC ) ){
        results << "test_computeCohesion (test 3) & False\n";
        return 1;
    }

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

    if ( error ){
        results << "test_computeCohesion & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( macroCJ, answerMacroC ) ){
        results << "test_computeCohesion (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microCJ, answerMicroC ) ){
        results << "test_computeCohesion (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( microGradientCJ, answerMicroGradientC ) ){
        results << "test_computeCohesion (test 6) & False\n";
        return 1;
    }

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

        if ( error ){
            results << "test_computeCohesion & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV - delta, microStrainISV, microGradientStrainISV,
                                                               macroHardeningParameters, microHardeningParameters,
                                                               microGradientHardeningParameters, macroCM, microCM,
                                                               microGradientCM );

        if ( error ){
            results << "test_computeCohesion & False\n";
            return 1;
        }

        variableType grad = ( macroCP - macroCM ) / ( 2 * delta );

        if ( !vectorTools::fuzzyEquals( grad, dMacroCdMacroStrainISV ) ){
            results << "test_computeCohesion (test 7) & False\n";
            return 1;
        }

        grad = ( microCP - microCM ) / ( 2 * delta );

        if ( !vectorTools::fuzzyEquals( grad, 0. ) ){
            results << "test_computeCohesion (test 8) & False\n";
            return 1;
        }

        variableVector gradCol = ( microGradientCP - microGradientCM ) / ( 2 * delta );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeCohesion (test 9) & False\n";
                return 1;
            }
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

        if ( error ){
            results << "test_computeCohesion & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV, microStrainISV - delta, microGradientStrainISV,
                                                               macroHardeningParameters, microHardeningParameters,
                                                               microGradientHardeningParameters, macroCM, microCM,
                                                               microGradientCM );

        if ( error ){
            results << "test_computeCohesion & False\n";
            return 1;
        }

        variableType grad = ( macroCP - macroCM ) / ( 2 * delta );

        if ( !vectorTools::fuzzyEquals( grad, 0. ) ){
            results << "test_computeCohesion (test 10) & False\n";
            return 1;
        }

        grad = ( microCP - microCM ) / ( 2 * delta );

        if ( !vectorTools::fuzzyEquals( grad, dMicroCdMicroStrainISV ) ){
            results << "test_computeCohesion (test 11) & False\n";
            return 1;
        }

        variableVector gradCol = ( microGradientCP - microGradientCM ) / ( 2 * delta );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                results << "test_computeCohesion (test 12) & False\n";
                return 1;
            }
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

        if ( error ){
            results << "test_computeCohesion & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computeCohesion( macroStrainISV, microStrainISV, microGradientStrainISV - delta,
                                                               macroHardeningParameters, microHardeningParameters,
                                                               microGradientHardeningParameters, macroCM, microCM,
                                                               microGradientCM );

        if ( error ){
            results << "test_computeCohesion & False\n";
            return 1;
        }

        variableType grad = ( macroCP - macroCM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( grad, 0. ) ){
            results << "test_computeCohesion (test 13) & False\n";
            return 1;
        }

        grad = ( microCP - microCM ) / ( 2 * delta[ i ] );

        if ( !vectorTools::fuzzyEquals( grad, 0. ) ){
            results << "test_computeCohesion (test 14) & False\n";
            return 1;
        }

        variableVector gradCol = ( microGradientCP - microGradientCM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMicroGradientCdMicroStrainISV[ j ][ i ] ) ){
                results << "test_computeCohesion (test 15) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeCohesion & True\n";
    return 0;
}

int test_evaluate_model( std::ofstream &results){
    /*!
     * Test the evaluation of the constitutive model.
     *
     * :param std::ofstream &results: The output file.
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

    solverTools::floatVector PK2Answer = { 177.067  , 13.287 ,  -0.577489,
                                            10.7222 , 152.621,  -0.288668,
                                            -1.38898, 1.61632, 150.85 };
    solverTools::floatVector SigmaAnswer = { 178.961  ,  13.5623 ,  -2.43027,
                                              13.5623 , 151.785  ,   1.62465,
                                              -2.43027,   1.62465, 149.31 };
    solverTools::floatVector MAnswer = { 0.541081, -0.533639,  0.640843,  2.92886 ,  1.1505  ,
                                         1.14946 ,  0.605546, -2.62366 ,  1.54942 , -2.40585 ,
                                        -0.670181, -0.848562,  0.678723,  0.433934, -0.206886,
                                        -2.7026  ,  0.964407,  1.68227 , -0.480138,  2.68016 ,
                                        -0.628373,  1.14616 , -0.10665 , -2.2393  , -0.765349,
                                         0.746722,  0.994415 };

    solverTools::floatVector SDVSAnswer = { -0.0394628  ,  0.102019   ,  0.0453858  ,  0.0262009 ,  0.0232167 ,
                                            -0.0196885  ,  0.0357449  ,  0.0453858  ,  0.0262009 ,  0.0232167 ,
                                             0.00813011 ,  0.00956498 , -0.000691731,  0.00885806, -0.0103998 ,
                                             0.000451969, -0.00115188 ,  0.00088041 , -0.00882951,  0.0622934 ,
                                             0.0274584  , -0.00113356 ,  0.0403186  , -0.00945181,  0.00119733,
                                            -0.00227071 ,  0.000880704, -0.00792983 ,  0.0554737 , -0.0202861 ,
                                            -0.0208186  ,  0.0122486  ,  0.0119279  ,  0.037049  ,  0.0168925 ,
                                            -0.0430597  , -0.0242213  ,  0.0244271  , -0.00880588,  0.0383104 ,
                                             0.00423276 ,  0.0111449  ,  0.0155982  ,  0.0136987 , -0.00813611,
                                            -0.0305008  ,  0.0171903  , -0.0126259  ,  0.0500443 , -0.0276203 ,
                                             0.0193062  ,  0.0310543  ,  0.0157931  ,  0.0168109 ,  0.0196658 };

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

    if ( errorCode != 0 ){
        std::cout << output_message;
        results << "test_evaluate_model & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SDVS, SDVSAnswer ) ){
        results << "test_evaluate_model (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( PK2Answer, current_PK2, 1e-5, 1e-5 ) ){
        std::cout << "error: "; vectorTools::print( PK2Answer - current_PK2 );
        results << "test_evaluate_model (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SigmaAnswer, current_SIGMA, 1e-5, 1e-5 ) ){
        results << "test_evaluate_model (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( MAnswer, current_M, 1e-5 ) ){
        results << "test_evaluate_model (test 4) & False\n";
        return 1;
    }

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
    if ( DEBUG.size() == 0 ){
        std::cout << "debug failure\n";
        results << "test_evaluate_model & False\n";
        return 1;
    }
#endif

    if ( errorCode != 0 ){
        std::cout << output_message;
        results << "test_evaluate_model & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SDVS, SDVSAnswer ) ){
        results << "test_evaluate_model (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( PK2Answer, current_PK2, 1e-5, 1e-5 ) ){
        results << "test_evaluate_model (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SigmaAnswer, current_SIGMA, 1e-5, 1e-5 ) ){
        results << "test_evaluate_model (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( MAnswer, current_M, 1e-5 ) ){
        results << "test_evaluate_model (test 8) & False\n";
        return 1;
    }

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

        if ( errorCode != 0 ){
            std::cout << output_message;
            results << "test_evaluate_model & False\n";
            return 1;
        }

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

        if ( errorCode != 0 ){
            std::cout << output_message;
            results << "test_evaluate_model & False\n";
            return 1;
        }

#ifdef DEBUG_MODE
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P[ "converged_values" ][ "converged_values" ].begin();
                   it_P != DEBUG_P[ "converged_values" ][ "converged_values" ].end(); it_P++ ){
            
            auto it_M = DEBUG_M[ "converged_values" ][ "converged_values" ].find( it_P->first );
            if ( it_M == DEBUG_M[ "converged_values" ][ "converged_values" ].end() ){
                std::cerr << "ERROR: A KEY EXISTS IN DEBUG_P THAT DOESNT EXIST IN DEBUG_M\n";
                results << "test_evaluate_model & False\n";
                return 1;
            }

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i / 3 ][ i % 3 ] ) );

        }

        if ( numericGradients.size() == 0 ){
            results << "test_evaluate_model & False\n";
            return 1;
        }

        //Check the total Jacobians of the plastic deformation measures
        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticDeformationGradient" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticDeformationGradient" ][ j ],
                                            DEBUG[ "totaldPlasticDeformationGradientdDeformationGradient" ][ 9 * j + i ],
                                            1e-6, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "convergedPlasticDeformationGradient" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPlasticDeformationGradientdDeformationGradient" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "convergedPlasticDeformationGradient" ][ j ]
                                        - DEBUG[ "totaldPlasticDeformationGradientdDeformationGradient" ][ 9 * j + i ] << "\n";
                results << "test_evaluate_model (dPlasticDeformationGradientdDeformationGradient) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticMicroDeformation" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticMicroDeformationdDeformationGradient" ][ 9 * j + i ],
                                            1e-6, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "convergedPlasticMicroDeformation" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPlasticMicroDeformationdDeformationGradient" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "convergedPlasticDeformationGradient" ][ j ]
                                        - DEBUG[ "totaldPlasticMicroDeformationdDeformationGradient" ][ 9 * j + i ] << "\n";
                results << "test_evaluate_model (dPlasticMicroDeformationdDeformationGradient) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticGradientMicroDeformation" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticGradientMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticGradientMicroDeformationdDeformationGradient" ][ 9 * j + i ],
                                            1e-6, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "convergedPlasticGradientMicroDeformation" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPlasticGradientMicroDeformationdDeformationGradient" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "convergedPlasticGradientDeformationGradient" ][ j ]
                                        - DEBUG[ "totaldPlasticGradientMicroDeformationdDeformationGradient" ][ 9 * j + i ] << "\n";
                results << "test_evaluate_model (dPlasticGradientMicroDeformationdDeformationGradient) & False\n";
                return 1;
            }
        }

        //Check the total Jacobians of the intermediate stresses
        for ( unsigned int j = 0; j < numericGradients[ "intermediatePK2Stress" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "intermediatePK2Stress" ][ j ],
                                            DEBUG[ "totaldPK2StressdDeformationGradient" ][ 9 * j + i ],
                                            1e-5, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "intermediatePK2Stress" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPK2StressdDeformationGradient" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "intermediatePK2Stress" ][ j ]
                                        - DEBUG[ "totaldPK2StressdDeformationGradient" ][ 9 * j + i ] << "\n";
                std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "totaldPK2StressdDeformationGradient" ], 9, 9 ) );
                results << "test_evaluate_model (dPK2StressdDeformationGradient) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceMicroStress" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceMicroStress" ][ j ],
                                            DEBUG[ "totaldReferenceMicroStressdDeformationGradient" ][ 9 * j + i ],
                                            1e-5, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "intermediateReferenceMicroStress" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldReferenceMicroStressdDeformationGradient" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "intermediateReferenceMicroStress" ][ j ]
                                        - DEBUG[ "totaldReferenceMicroStressdDeformationGradient" ][ 9 * j + i ] << "\n";
                std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "totaldReferenceMicroStressdDeformationGradient" ], 9, 9 ) );
                results << "test_evaluate_model (dReferenceMicroStressdDeformationGradient) & False\n";
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceHigherOrderStress" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceHigherOrderStress" ][ j ],
                                            DEBUG[ "totaldReferenceHigherOrderStressdDeformationGradient" ][ 9 * j + i ],
                                            1e-5, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "intermediateReferenceMicroStress" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldReferenceHigherOrderStressdDeformationGradient" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "intermediateReferenceMicroStress" ][ j ]
                                        - DEBUG[ "totaldReferenceHigherOrderStressdDeformationGradient" ][ 9 * j + i ] << "\n";
                std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "totaldReferenceHigherOrderStressdDeformationGradient" ], 27, 9 ) );
                results << "test_evaluate_model (dReferenceHigherOrderStressdDeformationGradient) & False\n";
            }
        }
#endif

        solverTools::floatVector gradCol = ( PK2_P - PK2_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DPK2Dgrad_u[ j ][ i ], 1e-5 ) ){
                std::cout << "row, column: " << j << ", " << i << "\n";
                std::cout << "num:   " << gradCol[ j ] << "\n";
                std::cout << "ana:   " << DPK2Dgrad_u[ j ][ i ] << "\n";
                std::cout << "error: " << gradCol[ j ] - DPK2Dgrad_u[ j ][ i ] << "\n";
                results << "test_evaluate_model (test 9) & False\n";
                return 1;
            }
        }

        gradCol = ( SIGMA_P - SIGMA_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DSIGMADgrad_u[ j ][ i ], 1e-5 ) ){
                std::cout << "row, column: " << j << ", " << i << "\n";
                std::cout << "num:   " << gradCol[ j ] << "\n";
                std::cout << "ana:   " << DSIGMADgrad_u[ j ][ i ] << "\n";
                std::cout << "error: " << gradCol[ j ] - DSIGMADgrad_u[ j ][ i ] << "\n";
                results << "test_evaluate_model (test 10) & False\n";
                return 1;
            }
        }

        gradCol = ( M_P - M_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DMDgrad_u[ j ][ i ], 1e-5 ) ){
                std::cout << "row, column: " << j << ", " << i << "\n";
                std::cout << "num:   " << gradCol[ j ] << "\n";
                std::cout << "ana:   " << DSIGMADgrad_u[ j ][ i ] << "\n";
                std::cout << "error: " << gradCol[ j ] - DMDgrad_u[ j ][ i ] << "\n";
                results << "test_evaluate_model (test 11) & False\n";
                return 1;
            }
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

        if ( errorCode != 0 ){
            std::cout << output_message;
            results << "test_evaluate_model & False\n";
            return 1;
        }

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

        if ( errorCode != 0 ){
            std::cout << output_message;
            results << "test_evaluate_model & False\n";
            return 1;
        }

#ifdef DEBUG_MODE
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P[ "converged_values" ][ "converged_values" ].begin();
                   it_P != DEBUG_P[ "converged_values" ][ "converged_values" ].end(); it_P++ ){
            
            auto it_M = DEBUG_M[ "converged_values" ][ "converged_values" ].find( it_P->first );
            if ( it_M == DEBUG_M[ "converged_values" ][ "converged_values" ].end() ){
                std::cerr << "ERROR: A KEY EXISTS IN DEBUG_P THAT DOESNT EXIST IN DEBUG_M\n";
                results << "test_evaluate_model & False\n";
                return 1;
            }

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );

        }

        if ( numericGradients.size() == 0 ){
            results << "test_evaluate_model & False\n";
            return 1;
        }

        //Check the total Jacobians of the plastic deformation measures
        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticDeformationGradient" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticDeformationGradient" ][ j ],
                                            DEBUG[ "totaldPlasticDeformationGradientdMicroDeformation" ][ 9 * j + i ],
                                            1e-6, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "convergedPlasticDeformationGradient" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPlasticDeformationGradientdMicroDeformation" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "convergedPlasticDeformationGradient" ][ j ]
                                        - DEBUG[ "totaldPlasticDeformationGradientdMicroDeformation" ][ 9 * j + i ] << "\n";
                results << "test_evaluate_model (dPlasticDeformationGradientdMicroDeformation) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticMicroDeformation" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticMicroDeformationdMicroDeformation" ][ 9 * j + i ],
                                            1e-6, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "convergedPlasticMicroDeformation" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPlasticMicroDeformationdMicroDeformation" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "convergedPlasticDeformationGradient" ][ j ]
                                        - DEBUG[ "totaldPlasticMicroDeformationdMicroDeformation" ][ 9 * j + i ] << "\n";
                results << "test_evaluate_model (dPlasticMicroDeformationdMicroDeformation) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticGradientMicroDeformation" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticGradientMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticGradientMicroDeformationdMicroDeformation" ][ 9 * j + i ],
                                            1e-6, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "convergedPlasticGradientMicroDeformation" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPlasticGradientMicroDeformationdMicroDeformation" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "convergedPlasticGradientDeformationGradient" ][ j ]
                                        - DEBUG[ "totaldPlasticGradientMicroDeformationdMicroDeformation" ][ 9 * j + i ] << "\n";
                results << "test_evaluate_model (dPlasticGradientMicroDeformationdMicroDeformation) & False\n";
                return 1;
            }
        }

        //Check the total Jacobians of the intermediate stresses
        for ( unsigned int j = 0; j < numericGradients[ "intermediatePK2Stress" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "intermediatePK2Stress" ][ j ],
                                            DEBUG[ "totaldPK2StressdMicroDeformation" ][ 9 * j + i ],
                                            1e-5, 1e-5 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "intermediatePK2Stress" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPK2StressdMicroDeformation" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "intermediatePK2Stress" ][ j ]
                                        - DEBUG[ "totaldPK2StressdMicroDeformation" ][ 9 * j + i ] << "\n";
                std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "totaldPK2StressdMicroDeformation" ], 9, 9 ) );
                results << "test_evaluate_model (dPK2StressdMicroDeformation) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceMicroStress" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceMicroStress" ][ j ],
                                            DEBUG[ "totaldReferenceMicroStressdMicroDeformation" ][ 9 * j + i ],
                                            1e-5, 1e-5 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "intermediateReferenceMicroStress" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldReferenceMicroStressdMicroDeformation" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "intermediateReferenceMicroStress" ][ j ]
                                        - DEBUG[ "totaldReferenceMicroStressdMicroDeformation" ][ 9 * j + i ] << "\n";
                std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "totaldReferenceMicroStressdMicroDeformation" ], 9, 9 ) );
                results << "test_evaluate_model (dReferenceMicroStressdMicroDeformation) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceHigherOrderStress" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceHigherOrderStress" ][ j ],
                                            DEBUG[ "totaldReferenceHigherOrderStressdMicroDeformation" ][ 9 * j + i ],
                                            1e-5, 1e-5 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "intermediateReferenceMicroStress" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldReferenceHigherOrderStressdMicroDeformation" ][ 9 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "intermediateReferenceMicroStress" ][ j ]
                                        - DEBUG[ "totaldReferenceHigherOrderStressdMicroDeformation" ][ 9 * j + i ] << "\n";
                std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "totaldReferenceHigherOrderStressdMicroDeformation" ], 27, 9 ) );
                results << "test_evaluate_model (dReferenceHigherOrderStressdMicroDeformation) & False\n";
                return 1;
            }
        }

#endif

        solverTools::floatVector gradCol = ( PK2_P - PK2_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DPK2Dphi[ j ][ i ], 1e-5, 1e-4 ) ){
                std::cout << "row, column: " << j << ", " << i << "\n";
                std::cout << "num:   " << gradCol[ j ] << "\n";
                std::cout << "ana:   " << DPK2Dphi[ j ][ i ] << "\n";
                std::cout << "error: " << gradCol[ j ] - DPK2Dphi[ j ][ i ] << "\n";
                results << "test_evaluate_model (test 12) & False\n";
                return 1;
            }
        }

        gradCol = ( SIGMA_P - SIGMA_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DSIGMADphi[ j ][ i ], 1e-5, 1e-4 ) ){
                std::cout << "row, column: " << j << ", " << i << "\n";
                std::cout << "num:   " << gradCol[ j ] << "\n";
                std::cout << "ana:   " << DSIGMADphi[ j ][ i ] << "\n";
                std::cout << "error: " << gradCol[ j ] - DSIGMADphi[ j ][ i ] << "\n";
                results << "test_evaluate_model (test 13) & False\n";
                return 1;
            }
        }

        gradCol = ( M_P - M_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DMDphi[ j ][ i ], 1e-5, 1e-4 ) ){
                std::cout << "row, column: " << j << ", " << i << "\n";
                std::cout << "num:   " << gradCol[ j ] << "\n";
                std::cout << "ana:   " << DMDphi[ j ][ i ] << "\n";
                std::cout << "error: " << gradCol[ j ] - DMDphi[ j ][ i ] << "\n";
                results << "test_evaluate_model (test 14) & False\n";
                return 1;
            }
        }
    }

    //Test the jacobians w.r.t. the gradient of the micro displacement
    for ( unsigned int i = 0; i < 27; i++ ){
        constantMatrix delta( 9, constantVector( 3, 0 ) );
        delta[ i / 3 ][ i % 3 ] = eps * fabs( current_grad_u[ i / 3 ][ i % 3 ] ) + eps;

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

        if ( errorCode != 0 ){
            std::cout << output_message;
            results << "test_evaluate_model & False\n";
            return 1;
        }

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

        if ( errorCode != 0 ){
            std::cout << output_message;
            results << "test_evaluate_model & False\n";
            return 1;
        }

#ifdef DEBUG_MODE
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P[ "converged_values" ][ "converged_values" ].begin();
                   it_P != DEBUG_P[ "converged_values" ][ "converged_values" ].end(); it_P++ ){
            
            auto it_M = DEBUG_M[ "converged_values" ][ "converged_values" ].find( it_P->first );
            if ( it_M == DEBUG_M[ "converged_values" ][ "converged_values" ].end() ){
                std::cerr << "ERROR: A KEY EXISTS IN DEBUG_P THAT DOESNT EXIST IN DEBUG_M\n";
                results << "test_evaluate_model & False\n";
                return 1;
            }

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i / 3 ][ i % 3 ] ) );

        }

        if ( numericGradients.size() == 0 ){
            results << "test_evaluate_model & False\n";
            return 1;
        }

        //Check the total Jacobians of the plastic deformation measures
        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticDeformationGradient" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticDeformationGradient" ][ j ],
                                            DEBUG[ "totaldPlasticDeformationGradientdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-6, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "convergedPlasticDeformationGradient" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPlasticDeformationGradientdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "convergedPlasticDeformationGradient" ][ j ]
                                        - DEBUG[ "totaldPlasticDeformationGradientdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                results << "test_evaluate_model (dPlasticDeformationGradientdGradientMicroDeformation) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticMicroDeformation" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticMicroDeformationdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-6, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "convergedPlasticMicroDeformation" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPlasticMicroDeformationdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "convergedPlasticDeformationGradient" ][ j ]
                                        - DEBUG[ "totaldPlasticMicroDeformationdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                results << "test_evaluate_model (dPlasticMicroDeformationdGradientMicroDeformation) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "convergedPlasticGradientMicroDeformation" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "convergedPlasticGradientMicroDeformation" ][ j ],
                                            DEBUG[ "totaldPlasticGradientMicroDeformationdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-6, 1e-8 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "convergedPlasticGradientMicroDeformation" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPlasticGradientMicroDeformationdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "convergedPlasticGradientDeformationGradient" ][ j ]
                                        - DEBUG[ "totaldPlasticGradientMicroDeformationdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                results << "test_evaluate_model (dPlasticGradientMicroDeformationdDeformationGradient) & False\n";
                return 1;
            }
        }

        //Check the total Jacobians of the intermediate stresses
        for ( unsigned int j = 0; j < numericGradients[ "intermediatePK2Stress" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "intermediatePK2Stress" ][ j ],
                                            DEBUG[ "totaldPK2StressdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-5, 1e-5 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "intermediatePK2Stress" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldPK2StressdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "intermediatePK2Stress" ][ j ]
                                        - DEBUG[ "totaldPK2StressdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "totaldPK2StressdGradientMicroDeformation" ], 9, 27 ) );
                results << "test_evaluate_model (dPK2StressdGradientMicroDeformation) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceMicroStress" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceMicroStress" ][ j ],
                                            DEBUG[ "totaldReferenceMicroStressdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-5, 1e-5 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "intermediateReferenceMicroStress" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldReferenceMicroStressdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "intermediateReferenceMicroStress" ][ j ]
                                        - DEBUG[ "totaldReferenceMicroStressdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "totaldReferenceMicroStressdGradientMicroDeformation" ], 9, 27 ) );
                results << "test_evaluate_model (dReferenceMicroStressdGradientMicroDeformation) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < numericGradients[ "intermediateReferenceHigherOrderStress" ].size(); j++ ){
            if ( !vectorTools::fuzzyEquals( numericGradients[ "intermediateReferenceHigherOrderStress" ][ j ],
                                            DEBUG[ "totaldReferenceHigherOrderStressdGradientMicroDeformation" ][ 27 * j + i ],
                                            1e-5, 1e-5 ) ){
                std::cout << "i, j:  " << i << ", " << j << "\n";
                std::cout << "num:   " << numericGradients[ "intermediateReferenceHigherOrderStress" ][ j ] << "\n";
                std::cout << "ana:   " << DEBUG[ "totaldReferenceHigherOrderStressdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                std::cout << "error: " << numericGradients[ "intermediateReferenceHigherOrderStress" ][ j ]
                                        - DEBUG[ "totaldReferenceHigherOrderStressdGradientMicroDeformation" ][ 27 * j + i ] << "\n";
                std::cout << "numeric:  "; vectorTools::print( numericGradients[ "intermediateReferenceHigherOrderStress" ] );
                std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "totaldReferenceHigherOrderStressdGradientMicroDeformation" ], 27, 27 ) );
                results << "test_evaluate_model (dReferenceHigherOrderStressdGradientMicroDeformation) & False\n";
                return 1;
            }
        }
#endif

        solverTools::floatVector gradCol = ( PK2_P - PK2_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DPK2Dgrad_phi[ j ][ i ] ) ){
                results << "test_evaluate_model (test 15) & False\n";
                return 1;
            }
        }

//        std::cout << "numeric DPK2Dgrad_u:\n"; vectorTools::print( gradCol );

        gradCol = ( SIGMA_P - SIGMA_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DSIGMADgrad_phi[ j ][ i ] ) ){
                results << "test_evaluate_model (test 16) & False\n";
                return 1;
            }
        }

//        std::cout << "numeric DSIGMADphi:\n"; vectorTools::print( gradCol );

        gradCol = ( M_P - M_M ) / ( 2 * delta[ i / 3 ][ i % 3 ] );
        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DMDgrad_phi[ j ][ i ] ) ){
                results << "test_evaluate_model (test 17) & False\n";
                return 1;
            }
        }

//        std::cout << "numeric DMDgrad_phi:\n"; vectorTools::print( gradCol );

    }

    results << "test_evaluate_model & True\n";
    return 1;
}

int test_computePlasticDeformationResidual( std::ofstream &results ){
    /*!
     * Test the computation of the plastic deformation residual.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_computePlasticDeformationResidual & False\n";
        return 1;
    }

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

    //Compute the current values of the cohesion
    variableType currentMacroCohesion, currentMicroCohesion;
    variableVector currentMicroGradientCohesion;

    error = micromorphicElastoPlasticity::computeCohesion( currentMacroStrainISV, currentMicroStrainISV,
                                                           currentMicroGradientStrainISV,
                                                           macroHardeningParameters, microHardeningParameters,
                                                           microGradientHardeningParameters,
                                                           currentMacroCohesion, currentMicroCohesion,
                                                           currentMicroGradientCohesion );

    if ( error ){
        error->print();
        results << "test_computePlasticDeformationResidual & False\n";
        return 1;
    }

    solverTools::floatMatrix floatArgsDefault =
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

    variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

    solverTools::floatMatrix floatOutsDefault = 
        {
            currentPK2Stress,
            currentReferenceMicroStress,
            currentReferenceHigherOrderStress
        };

    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intArgs = { { 0 } };
    solverTools::intMatrix intOutsDefault = { { 1, 1, 0, 0, 0 } };
    
    solverTools::intMatrix intOuts = intOutsDefault;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG;
#endif

    solverTools::floatVector x( 55, 0 );
    for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
        x[ i + 0 ] = currentPlasticDeformationGradient[ i ];
        x[ i + 9 ] = currentPlasticMicroDeformation[ i ];
    }
    for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
        x[ i + 18 ] = currentPlasticGradientMicroDeformation[ i ];
    }
    x[ 45 ] = currentMacroStrainISV;
    x[ 46 ] = currentMicroStrainISV;
    for ( unsigned int i = 0; i < currentMicroGradientStrainISV.size(); i++ ){
        x[ 47 + i ] = currentMicroGradientStrainISV[ i ];
    }
    x[ 50 ] = currentMacroGamma;
    x[ 51 ] = currentMicroGamma;
    for ( unsigned int i = 0; i < currentMicroGradientGamma.size(); i++ ){
        x[ 52 + i ] = currentMicroGradientGamma[ i ];
    }

    solverTools::floatVector residual;
    solverTools::floatMatrix jacobian;
    solverTools::floatVector residualAnswer = { -1.81569518e-02,  1.50135543e-02, -8.34170711e-03,  1.43829032e-02,
                                                 3.06498331e-03, -2.53977877e-02, -1.52215981e-02, -1.99094335e-02,
                                                -6.94793308e-03,  4.09965704e-03,  2.40258795e-03, -2.14047929e-02,
                                                 1.89482872e-02,  7.62350440e-03, -2.85891306e-03,  2.79109725e-03,
                                                -2.08777378e-02,  5.99160625e-03, -9.42067987e-04, -1.07382946e-02,
                                                 6.81373533e-03, -6.78480101e-03,  1.35221379e-02,  8.21125674e-03,
                                                 1.96686396e-02,  1.36813595e-02,  1.59185809e-02, -8.00705013e-03,
                                                -1.38683197e-03,  2.28972653e-03, -9.66137216e-03, -1.65581701e-02,
                                                 5.13825136e-03, -1.61062277e-02,  2.69112683e-02,  1.88244553e-02,
                                                 8.65099576e-03,  1.10231885e-03,  1.16071302e-02, -1.42204396e-02,
                                                 2.44683133e-02,  2.67680522e-03,  3.13831236e-02,  1.74318178e-02,
                                                 3.07478347e-03, -2.00436336e-02, -5.70816527e-02, -1.94400838e-02,
                                                -3.71128873e-02, -5.47856908e-02,  3.84911950e+02,  4.73427460e+02,
                                                -2.91632000e+00, -5.56867000e+00, -8.22474000e+00 };

    error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                             floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                             , DEBUG
#endif
                                                                             );

    if ( error ){
        error->print();
        results << "test_computePlasticDeformationResidual & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( residualAnswer, residual ) ){
        std::cout << "residual error\n";
        vectorTools::print( residualAnswer - residual );
        results << "test_computePlasticDeformationResidual (test 1) & False\n";
        return 1;
    }

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

        if ( error ){
            error->print();
            results << "test_computePlasticDeformationResidual & False\n";
            return 1;
        }

        fO = floatOutsDefault;
        iO = intOutsDefault;

        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x - delta, fA, iA, residual_M, jacobian_M,
                                                                                 fO, iO
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        if ( error ){
            error->print();
            results << "test_computePlasticDeformationResidual & False\n";
            return 1;
        }

#ifdef DEBUG_MODE

        //Debug each of the sub-jacobians if required. This can be very slow so it isn't done all the time

        //Assemble the numeric Jacobians
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P.begin(); it_P != DEBUG_P.end(); it_P++ ){
            
            auto it_M = DEBUG_M.find( it_P->first );
            if ( it_M == DEBUG_M.end() ){
                std::cerr << "ERROR: A KEY EXISTS IN DEBUG_P THAT DOESNT EXIST IN DEBUG_M\n";
                results << "test_computePlasticDeformationResidual & False\n";
                return 1;
            }

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );
        }

        if ( i < 9 ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                DEBUG[ "dElasticDeformationGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticDeformationGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroDeformationdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticGradientMicroDeformationdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticRightCauchyGreendPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticRightCauchyGreendPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroRightCauchyGreendPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticPsidPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticGammadPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dPK2StressdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceMicroStressdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceHigherOrderStressdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMacroCohesiondPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMicroCohesiondPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMicroGradientCohesiondPlasticDeformationGradient) & False\n";
                return 1;
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroFlowDirectiondPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroFlowDirectiondPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-8, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroGradientFlowDirectiondPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMacroStrainISVdPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMicroStrainISVdPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMicroGradientStrainISVdPlasticDeformationGradient) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticMacroVelocityGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticMicroVelocityGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-8 ) ){
                    std::cout << "ana: " << DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ] << "\n";
                    std::cout << "res: " << numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ] << "\n";
                    results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedPlasticMicroDeformationdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient ) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                            DEBUG[ "dMacroYielddPlasticDeformationGradient" ][ i ],
                                            1e-5, 1e-8 ) ){
                results << "test_computePlasticDeformationResidual (dMacroYielddPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                            DEBUG[ "dMicroYielddPlasticDeformationGradient" ][ i ],
                                            1e-5, 1e-8 ) ){
                results << "test_computePlasticDeformationResidual (dMicroYielddPlasticDeformationGradient) & False\n";
                return 1;
            }

            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                DEBUG[ "dMicroGradientYielddPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroGradientYielddPlasticDeformationGradient ) & False\n";
                    return 1;
                }
            }
        }

        if ( ( i >= 9 ) && ( i < 18 ) ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticDeformationGradientdPlasticMicroDeformation) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroDeformationdPlasticMicroDeformation) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticGradientMicroDeformationdPlasticMicroDeformation) & False\n";

                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticRightCauchyGreendPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticMicroRightCauchyGreendPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroRightCauchyGreendPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticPsidPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    vectorTools::print( numericGradients[ "currentElasticGamma" ] );
                    results << "test_computePlasticDeformationResidual (dElasticGammadPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dPK2StressdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceMicroStressdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceHigherOrderStressdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMacroCohesiondPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMicroCohesiondPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMicroGradientCohesiondPlasticMicroDeformation) & False\n";
                return 1;
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroFlowDirectiondPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroFlowDirectiondPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-8, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroGradientFlowDirectiondPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMacroStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMicroStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMicroGradientStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticMacroVelocityGradientdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticMicroVelocityGradientdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroVelocityGradient" ] );
                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation" ], 27, 9 ) );
                    results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedPlasticMicroDeformationdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-8 ) ){
                    std::cout << "num: " << numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ] << "\n";
                    std::cout << "ana: " << DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ] << "\n";
                    results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                            DEBUG[ "dMacroYielddPlasticMicroDeformation" ][ i - 9 ],
                                            1e-5, 1e-8 ) ){
                results << "test_computePlasticDeformationResidual (dMacroYielddPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                            DEBUG[ "dMicroYielddPlasticMicroDeformation" ][ i - 9 ],
                                            1e-5, 1e-8 ) ){
                results << "test_computePlasticDeformationResidual (dMicroYielddPlasticMicroDeformation) & False\n";
                return 1;
            }

            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                DEBUG[ "dMicroGradientYielddPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-5, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroGradientYielddPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }
        }

        if ( ( i >= 18 ) && ( i < 45 ) ){
            
            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticDeformationGradientdPlasticGradientMicroDeformation) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroDeformationdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticGradientMicroDeformationdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticRightCauchyGreendPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroRightCauchyGreendPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticPsidPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticGammadPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) ){
                    vectorTools::print( numericGradients[ "currentPK2Stress" ] );
                    vectorTools::print( vectorTools::inflate( DEBUG[ "dPK2StressdPlasticGradientMicroDeformation" ], 27, 9 ) );
                    results << "test_computePlasticDeformationResidual (dPK2StressdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceMicroStressdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceHigherOrderStressdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMacroCohesiondPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMicroCohesiondPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMicroGradientCohesiondPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroFlowDirectiondPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroFlowDirectiondPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroGradientFlowDirectiondPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMacroStrainISVdPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMicroStrainISVdPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMicroGradientStrainISVdPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-7 ) ){
                    std::cout << "ana: " << DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ] << "\n";
                    std::cout << "num: " << numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ] << "\n";
                    std::cout << "error: " << DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ] - numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ] << "\n";
                    results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                            DEBUG[ "dMacroYielddPlasticGradientMicroDeformation" ][ i - 18 ],
                                            1e-5, 1e-8 ) ){
                results << "test_computePlasticDeformationResidual (dMacroYielddPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                            DEBUG[ "dMicroYielddPlasticGradientMicroDeformation" ][ i - 18 ],
                                            1e-5, 1e-8 ) ){
                std::cout << "i, j: " << i << "\n";
                std::cout << numericGradients[ "microYieldFunction" ][ 0 ] << "\n";
                std::cout << DEBUG[ "dMicroYielddPlasticGradientMicroDeformation" ][ i - 18 ] << "\n";
                results << "test_computePlasticDeformationResidual (dMicroYielddPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                DEBUG[ "dMicroGradientYielddPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroGradientYielddPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }
        }

        if ( ( i >= 45 ) && ( i < 50 ) ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticDeformationGradientdStrainLikeISVs) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroDeformationdStrainLikeISVs) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticGradientMicroDeformationdStrainLikeISVs) & False\n";

                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticRightCauchyGreendStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroRightCauchyGreendStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticPsidStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    vectorTools::print( numericGradients[ "currentElasticGamma" ] );
                    results << "test_computePlasticDeformationResidual (dElasticGammadStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPK2StressdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceMicroStressdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceHigherOrderStressdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( i == 45 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                                DEBUG[ "dMacroCohesiondMacroStrainISV" ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dCurrentMacroCohesiondMacroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dCurrentMicroCohesiondMacroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                                variableVector( 3, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dCurrentMicroGradientCohesiondMacroStrainISV) & False\n";
                    return 1;
                }
            }
            else if ( i == 46 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dCurrentMacroCohesiondMicroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                                DEBUG[ "dMicroCohesiondMicroStrainISV" ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dCurrentMicroCohesiondMicroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                                variableVector( 3, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dCurrentMicroGradientCohesiondMicroStrainISV) & False\n";
                    return 1;
                }
            }
            else if ( i > 46 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dCurrentMacroCohesiondMicroGradientStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dCurrentMicroCohesiondMicroGradientStrainISV) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < 3; j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ][ j ],
                                                    DEBUG[ "dMicroGradientCohesiondMicroGradientStrainISV" ][ 3 * j + i - 47 ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dCurrentMicroGradientCohesiondMicroGradientStrainISV) & False\n";
                        return 1;
                    }
                }
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroFlowDirectiondStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroFlowDirectiondStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroGradientFlowDirectiondStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMacroStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMicroStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dExpectedMicroGradientStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticMacroVelocityGradientdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticMicroVelocityGradientdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdStrainLikeISVs & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedPlasticMicroDeformationdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroDeformationdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( i == 45 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
                                                DEBUG[ "dMacroYielddMacroStrainISV" ],
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroYielddMacroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                                0.,
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroYielddMacroStrainISV) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual (dMacroGradientYielddMacroStrainISV) & False\n";
                        return 1;
                    }
                }
            }
            else if ( i == 46 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                                0.,
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroYielddMicroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
                                                DEBUG[ "dMicroYielddMicroStrainISV" ],
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroYielddMicroStrainISV) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual (dMacroGradientYielddMicroStrainISV) & False\n";
                        return 1;
                    }
                }
            }
            else if ( i >= 47 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                                0.,
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroYielddMicroGradientStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                                0.,
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroYielddMicroGradientStrainISV) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    DEBUG[ "dMicroGradientYielddMicroGradientStrainISV" ][ 3 * j + i - 47 ],
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual (dMacroGradientYielddMicroGradientStrainISV) & False\n";
                        return 1;
                    }
                }
            } 
        }

        if ( ( i >= 50 ) && ( i < 55 ) ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticDeformationGradientdGammas) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroDeformationdGammas) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticGradientMicroDeformationdGammas) & False\n";

                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticRightCauchyGreendGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticMicroRightCauchyGreendGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dElasticPsidGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    vectorTools::print( numericGradients[ "currentElasticGamma" ] );
                    results << "test_computePlasticDeformationResidual (dElasticGammadGammas) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dPK2StressdGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceMicroStressdGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dReferenceHigherOrderStressdGammas) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMacroCohesiondGammas) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMicroCohesiondGammas) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual (dCurrentMicroGradientCohesiondGammas) & False\n";
                return 1;
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroFlowDirectiondGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroFlowDirectiondGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroGradientFlowDirectiondGammas) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( i == 50 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                                DEBUG[ "dExpectedMacroISVdMacroGamma" ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedMacroStrainISVdMacroGamma) & False\n";
                    return 1;
                }

                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedMicroStrainISVdMacroGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                                variableVector( 3, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedMicroGradientStrainISVdMacroGamma) & False\n";
                    return 1;
                }
            }

            if ( i == 51 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedMacroStrainISVdMicroGamma) & False\n";
                    return 1;
                }

                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                                DEBUG[ "dExpectedMicroISVdMicroGamma" ],
                                                1e-9, 1e-9 ) ){
                    std::cout << "ana: " << DEBUG[ "dExpectedMicroISVdMicroGamma" ][ 0 ] << "\n";
                    std::cout << "num: " << numericGradients[ "expectedMicroStrainISV" ][ 0 ] << "\n";
                    std::cout << "err: " << DEBUG[ "dExpectedMicroISVdMicroGamma" ][ 0 ] - numericGradients[ "expectedMicroStrainISV" ][ 0 ] << "\n";
                    results << "test_computePlasticDeformationResidual (dExpectedMicroStrainISVdMicroGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                                variableVector( 3, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedMicroGradientStrainISVdMicroGamma) & False\n";
                    return 1;
                }
            }

            if ( i >= 52 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedMacroStrainISVdMicroGradientGamma) & False\n";
                    return 1;
                }

                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual (dExpectedMicroStrainISVdMicroGradientGamma) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < 3; j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ][ j ],
                                                    DEBUG[ "dExpectedMicroGradientISVdMicroGradientGamma" ][ 3 * j + i - 52 ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dExpectedMicroGradientStrainISVdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
            }


            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            if ( i == 50 ){

                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMacroVelocityGradientdMacroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dPlasticMacroVelocityGradientdMacroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dPlasticMicroVelocityGradientdMacroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdMacroGamma) & False\n";
                        return 1;
                    }
                }
            }

            if ( i == 51 ){

                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMacroVelocityGradientdMicroGamma" ][ j ],
                                                    1e-9, 1e-8 ) ){
                        std::cout << "num: " << numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ] << "\n";
                        std::cout << "ana: " << DEBUG[ "dPlasticMacroVelocityGradientdMicroGamma" ][ j ] << "\n";
                        std::cout << "err: " << numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ] - DEBUG[ "dPlasticMacroVelocityGradientdMicroGamma" ][ j ] << "\n";
                        results << "test_computePlasticDeformationResidual (dPlasticMacroVelocityGradientdMicroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMicroVelocityGradientdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        std::cout << "num: " << numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ] << "\n";
                        std::cout << "ana: " << DEBUG[ "dPlasticMicroVelocityGradientdMicroGamma" ][ j ] << "\n";
                        std::cout << "err: " << numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ] - DEBUG[ "dPlasticMicroVelocityGradientdMicroGamma" ][ j ] << "\n";
                        results << "test_computePlasticDeformationResidual (dPlasticMicroVelocityGradientdMicroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdMicroGamma) & False\n";
                        return 1;
                    }
                }
            }

            if ( i >= 52 ){

                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dPlasticMacroVelocityGradientdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dPlasticMicroVelocityGradientdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGradientGamma" ][ 3 * j + i - 52 ],
                                                    1e-9, 1e-9 ) ){
                        std::cout << "num1: " << numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ] << "\n";
                        std::cout << "ana1: " << DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGradientGamma" ][ 3 * j + i - 52 ] << "\n";
                        std::cout << "err1: " << numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ] - DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGradientGamma" ][ 3 * j + i - 52 ] << "\n";
                        results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            if ( i == 50 ){

                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                    DEBUG[ "dExpectedPlasticDeformationGradientdMacroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroVelocityGradientdMacroGamma & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticMicroDeformationdMacroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMacroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroDeformationdMacroGamma) & False\n";
                        return 1;
                    }
                }
            }

            else if ( i == 51 ){

                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                    DEBUG[ "dExpectedPlasticDeformationGradientdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroVelocityGradientdMicroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                    DEBUG[ "dExpectedPlasticMicroDeformationdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticMicroDeformationdMicroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        std::cout << "numeric:  "; vectorTools::print( numericGradients[ "expectedPlasticGradientMicroDeformation" ] );
                        std::cout << "analytic: "; vectorTools::print( DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGamma" ] );
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroDeformationdMicroGamma) & False\n";
                        return 1;
                    }
                }
            }

            else if ( i >= 52 ){

                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroVelocityGradientdMicroGradientGamma & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticMicroDeformationdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGradientGamma" ][ 3 * j + i - 52 ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroDeformationdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( i == 50 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
                                                variableVector( 1, 0 ),
                                                1e-5, 1e-8 ) ){
                    std::cout << "numeric:  "; vectorTools::print( numericGradients[ "macroYieldFunction" ] );
                    std::cout << "analytic: "; vectorTools::print( DEBUG[ "dMacroYielddMacroGamma" ] );
                    results << "test_computePlasticDeformationResidual (dMacroYielddMacroGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
                                                variableVector( 1, 0 ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroYielddMacroMacroGamma) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual (dMacroGradientYielddMacroGamma) & False\n";
                        return 1;
                    }
                }
            }
            else if ( i == 51 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
                                                variableVector( 1, 0. ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroYielddMicroGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
                                                variableVector( 1, 0. ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroYielddMicroGamma) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual (dMacroGradientYielddMicroStrainISV) & False\n";
                        return 1;
                    }
                }
            }
            else if ( i >= 52 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
                                                variableVector( 1, 0. ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMacroYielddMicroGradientGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
                                                variableVector( 1, 0. ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual (dMicroYielddMicroGradientGamma) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual (dMacroGradientYielddMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
            } 

        }

        if ( i >= 55 ){
            std::cout << "ERROR in the test!\n";
            assert( 1 == 0 );
        }

#endif


        //Jacobian test

        solverTools::floatVector gradCol = ( residual_P - residual_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ], 1e-5, 1e-7 ) ){
                std::cout << "i, j: " << j << ", " << i << "\n";
                std::cout << "gradCol[ j ]: " << gradCol[ j ] << "\n";
                std::cout << "jacobian[ j ][ i ]: " << jacobian[ j ][ i ] << "\n";
                std::cout << "error: " << gradCol[ j ] - jacobian[ j ][ i ] << "\n";
                std::cout << "gradCol:\n"; vectorTools::print( gradCol );
                std::cout << "jacobian:\n"; vectorTools::print( jacobian );
                results << "test_computePlasticDeformationResidual (test 2) & False\n";
                return 1;
            }
        }
    }

    results << "test_computePlasticDeformationResidual & True\n";
    return 0;
}

int test_computePlasticDeformationResidual2( std::ofstream &results ){
    /*!
     * Second test the computation of the plastic deformation residual.
     *
     * :param std::ofstream &results: The output file.
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

    if ( error ){
        error->print();
        results << "test_computePlasticDeformationResidual2 & False\n";
        return 1;
    }

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

    //Compute the current values of the cohesion
    variableType currentMacroCohesion, currentMicroCohesion;
    variableVector currentMicroGradientCohesion;

    error = micromorphicElastoPlasticity::computeCohesion( currentMacroStrainISV, currentMicroStrainISV,
                                                           currentMicroGradientStrainISV,
                                                           macroHardeningParameters, microHardeningParameters,
                                                           microGradientHardeningParameters,
                                                           currentMacroCohesion, currentMicroCohesion,
                                                           currentMicroGradientCohesion );

    if ( error ){
        error->print();
        results << "test_computePlasticDeformationResidual2 & False\n";
        return 1;
    }

    solverTools::floatMatrix floatArgsDefault =
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

    variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

    solverTools::floatMatrix floatOutsDefault = 
        {
            currentPK2Stress,
            currentReferenceMicroStress,
            currentReferenceHigherOrderStress,
            {}, {}, {}
        };

    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intArgs = { { 1 } };
    solverTools::intMatrix intOutsDefault = { { 1, 1, 1, 1, 1 } };
    
    solverTools::intMatrix intOuts = intOutsDefault;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG;
#endif

    solverTools::floatVector x( 55, 0 );
    for ( unsigned int i = 0; i < currentPlasticDeformationGradient.size(); i++ ){
        x[ i + 0 ] = currentPlasticDeformationGradient[ i ];
        x[ i + 9 ] = currentPlasticMicroDeformation[ i ];
    }
    for ( unsigned int i = 0; i < currentPlasticGradientMicroDeformation.size(); i++ ){
        x[ i + 18 ] = currentPlasticGradientMicroDeformation[ i ];
    }
    x[ 45 ] = currentMacroStrainISV;
    x[ 46 ] = currentMicroStrainISV;
    for ( unsigned int i = 0; i < currentMicroGradientStrainISV.size(); i++ ){
        x[ 47 + i ] = currentMicroGradientStrainISV[ i ];
    }
    x[ 50 ] = currentMacroGamma;
    x[ 51 ] = currentMicroGamma;
    for ( unsigned int i = 0; i < currentMicroGradientGamma.size(); i++ ){
        x[ 52 + i ] = currentMicroGradientGamma[ i ];
    }

    solverTools::floatVector residual;
    solverTools::floatMatrix jacobian;
    solverTools::floatVector residualAnswer = { 2.40362117e-03,  2.48922816e-02, -5.83893367e-04,  2.65403605e-02,
                                               -3.08798343e-04,  1.69617713e-04,  3.97561325e-03, -3.24284857e-03,
                                               -2.73722736e-02, -2.51136874e-03,  7.78436675e-03, -1.99632443e-02,
                                               -1.32361861e-02,  8.11936091e-03,  7.91525734e-03,  1.38718574e-02,
                                                1.63576816e-02,  8.37538047e-03,  1.47178162e-02, -1.55973801e-03,
                                                2.52508602e-02, -9.90890790e-04, -7.81869562e-03,  6.64541346e-03,
                                                1.21008476e-02,  1.65997331e-02, -5.34466002e-03, -1.31382242e-03,
                                                3.80903918e-02,  2.10651777e-04,  3.08104546e-03,  5.60400074e-03,
                                                2.49141092e-02, -8.94189185e-03,  1.12660558e-02, -1.19963319e-02,
                                                3.40284172e-02,  8.00466552e-03,  5.02773187e-04,  9.96908538e-03,
                                                6.93297045e-03, -3.84414842e-03, -1.07296826e-03,  6.43585721e-04,
                                               -2.58712059e-03, -2.00436336e-02, -5.70816527e-02, -1.94400839e-02,
                                               -3.71128873e-02, -5.47856908e-02,  3.85689555e+02,  4.56790182e+02,
                                               -2.92086025e+00, -5.57253093e+00, -8.22748598e+00 };

    error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                             floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                             , DEBUG
#endif
                                                                             );

    if ( error ){
        error->print();
        results << "test_computePlasticDeformationResidual2 & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( residualAnswer, residual, 1e-5 ) ){
        std::cout << "residual:\n"; vectorTools::print( residual );
        std::cout << "error:\n";
        vectorTools::print( residualAnswer - residual );
        results << "test_computePlasticDeformationResidual2 (test 1) & False\n";
        return 1;
    }

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

        if ( error ){
            error->print();
            results << "test_computePlasticDeformationResidual2 & False\n";
            return 1;
        }

        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x - delta, fA, iA, residual_M, jacobian_M,
                                                                                 fO_M, iO_M
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        if ( error ){
            error->print();
            results << "test_computePlasticDeformationResidual2 & False\n";
            return 1;
        }

#ifdef DEBUG_MODE

        //Debug each of the sub-jacobians if required. This can be very slow so it isn't done all the time

        //Assemble the numeric Jacobians
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P.begin(); it_P != DEBUG_P.end(); it_P++ ){
            
            auto it_M = DEBUG_M.find( it_P->first );
            if ( it_M == DEBUG_M.end() ){
                std::cerr << "ERROR: A KEY EXISTS IN DEBUG_P THAT DOESNT EXIST IN DEBUG_M\n";
                results << "test_computePlasticDeformationResidual2 & False\n";
                return 1;
            }

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );
        }

        if ( i < 9 ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                DEBUG[ "dElasticDeformationGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticDeformationGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroDeformationdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGradientMicroDeformationdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticRightCauchyGreendPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticRightCauchyGreendPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroRightCauchyGreendPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticPsidPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGammadPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPK2StressdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceMicroStressdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceHigherOrderStressdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMacroCohesiondPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMicroCohesiondPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMicroGradientCohesiondPlasticDeformationGradient) & False\n";
                return 1;
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroFlowDirectiondPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroFlowDirectiondPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-8, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientFlowDirectiondPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdPlasticDeformationGradient) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdPlasticDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdPlasticDeformationGradient ) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                            DEBUG[ "dMacroYielddPlasticDeformationGradient" ][ i ],
                                            1e-5, 1e-8 ) ){
                results << "test_computePlasticDeformationResidual2 (dMacroYielddPlasticDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                            DEBUG[ "dMicroYielddPlasticDeformationGradient" ][ i ],
                                            1e-5, 1e-8 ) ){
                results << "test_computePlasticDeformationResidual2 (dMicroYielddPlasticDeformationGradient) & False\n";
                return 1;
            }

            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                DEBUG[ "dMicroGradientYielddPlasticDeformationGradient" ][ 9 * j + i ],
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientYielddPlasticDeformationGradient ) & False\n";
                    return 1;
                }
            }
        }

        if ( ( i >= 9 ) && ( i < 18 ) ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticDeformationGradientdPlasticMicroDeformation) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroDeformationdPlasticMicroDeformation) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGradientMicroDeformationdPlasticMicroDeformation) & False\n";

                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticRightCauchyGreendPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticMicroRightCauchyGreendPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroRightCauchyGreendPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticPsidPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    vectorTools::print( numericGradients[ "currentElasticGamma" ] );
                    results << "test_computePlasticDeformationResidual2 (dElasticGammadPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPK2StressdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceMicroStressdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceHigherOrderStressdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMacroCohesiondPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMicroCohesiondPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMicroGradientCohesiondPlasticMicroDeformation) & False\n";
                return 1;
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroFlowDirectiondPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroFlowDirectiondPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-8, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientFlowDirectiondPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroVelocityGradient" ] );
                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation" ], 27, 9 ) );
                    results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdPlasticMicroDeformation & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    std::cout << "ana: " << DEBUG[ "dExpectedPlasticMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ] << "\n";
                    std::cout << "res: " << numericGradients[ "expectedPlasticMicroDeformation" ][ j ] << "\n";
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                            DEBUG[ "dMacroYielddPlasticMicroDeformation" ][ i - 9 ],
                                            1e-5, 1e-8 ) ){
                results << "test_computePlasticDeformationResidual2 (dMacroYielddPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                            DEBUG[ "dMicroYielddPlasticMicroDeformation" ][ i - 9 ],
                                            1e-5, 1e-8 ) ){
                results << "test_computePlasticDeformationResidual2 (dMicroYielddPlasticMicroDeformation) & False\n";
                return 1;
            }

            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                DEBUG[ "dMicroGradientYielddPlasticMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-5, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroGradientYielddPlasticMicroDeformation) & False\n";
                    return 1;
                }
            }
        }

        if ( ( i >= 18 ) && ( i < 45 ) ){
            
            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticDeformationGradientdPlasticGradientMicroDeformation) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroDeformationdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGradientMicroDeformationdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticRightCauchyGreendPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroRightCauchyGreendPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticPsidPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGammadPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-7 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    std::cout << "num:   " << numericGradients[ "currentPK2Stress" ][ j ] << "\n";
                    std::cout << "ana:   " << DEBUG[ "dPK2StressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ] << "\n";
                    std::cout << "error: " << numericGradients[ "currentPK2Stress" ][ j ] -  DEBUG[ "dPK2StressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ] << "\n";
                    vectorTools::print( numericGradients[ "currentPK2Stress" ] );
                    vectorTools::print( vectorTools::inflate( DEBUG[ "dPK2StressdPlasticGradientMicroDeformation" ], 27, 9 ) );
                    results << "test_computePlasticDeformationResidual2 (dPK2StressdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceMicroStressdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceHigherOrderStressdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMacroCohesiondPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMicroCohesiondPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMicroGradientCohesiondPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroFlowDirectiondPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroFlowDirectiondPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientFlowDirectiondPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-8 ) ){
                    std::cout << "ana: " << DEBUG[ "dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ] << "\n";
                    std::cout << "num: " << numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ] << "\n";
                    results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticDeformationGradientdPlasticGradientMicroDeformation & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                            DEBUG[ "dMacroYielddPlasticGradientMicroDeformation" ][ i - 18 ],
                                            1e-5, 1e-7 ) ){
                results << "test_computePlasticDeformationResidual2 (dMacroYielddPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                            DEBUG[ "dMicroYielddPlasticGradientMicroDeformation" ][ i - 18 ],
                                            1e-5, 1e-7 ) ){
                std::cout << "i, j: " << i << "\n";
                std::cout << numericGradients[ "microYieldFunction" ][ 0 ] << "\n";
                std::cout << DEBUG[ "dMicroYielddPlasticGradientMicroDeformation" ][ i - 18 ] << "\n";
                results << "test_computePlasticDeformationResidual2 (dMicroYielddPlasticGradientMicroDeformation) & False\n";
                return 1;
            }

            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                DEBUG[ "dMicroGradientYielddPlasticGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroGradientYielddPlasticGradientMicroDeformation) & False\n";
                    return 1;
                }
            }
        }

        if ( ( i >= 45 ) && ( i < 50 ) ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticDeformationGradientdStrainLikeISVs) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroDeformationdStrainLikeISVs) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGradientMicroDeformationdStrainLikeISVs) & False\n";

                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticRightCauchyGreendStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroRightCauchyGreendStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticPsidStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    vectorTools::print( numericGradients[ "currentElasticGamma" ] );
                    results << "test_computePlasticDeformationResidual2 (dElasticGammadStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPK2StressdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceMicroStressdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceHigherOrderStressdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( i == 45 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                                DEBUG[ "dMacroCohesiondMacroStrainISV" ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dCurrentMacroCohesiondMacroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dCurrentMicroCohesiondMacroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                                variableVector( 3, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dCurrentMicroGradientCohesiondMacroStrainISV) & False\n";
                    return 1;
                }
            }
            else if ( i == 46 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dCurrentMacroCohesiondMicroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                                DEBUG[ "dMicroCohesiondMicroStrainISV" ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dCurrentMicroCohesiondMicroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                                variableVector( 3, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dCurrentMicroGradientCohesiondMicroStrainISV) & False\n";
                    return 1;
                }
            }
            else if ( i > 46 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dCurrentMacroCohesiondMicroGradientStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dCurrentMicroCohesiondMicroGradientStrainISV) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < 3; j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ][ j ],
                                                    DEBUG[ "dMicroGradientCohesiondMicroGradientStrainISV" ][ 3 * j + i - 47 ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dCurrentMicroGradientCohesiondMicroGradientStrainISV) & False\n";
                        return 1;
                    }
                }
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroFlowDirectiondStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroFlowDirectiondStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientFlowDirectiondStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdPlasticMicroDeformation) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdStrainLikeISVs & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdStrainLikeISVs) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( i == 45 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
                                                DEBUG[ "dMacroYielddMacroStrainISV" ],
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroYielddMacroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                                0.,
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroYielddMacroStrainISV) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual2 (dMacroGradientYielddMacroStrainISV) & False\n";
                        return 1;
                    }
                }
            }
            else if ( i == 46 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                                0.,
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroYielddMicroStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
                                                DEBUG[ "dMicroYielddMicroStrainISV" ],
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroYielddMicroStrainISV) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual2 (dMacroGradientYielddMicroStrainISV) & False\n";
                        return 1;
                    }
                }
            }
            else if ( i >= 47 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ 0 ],
                                                0.,
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroYielddMicroGradientStrainISV) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ 0 ],
                                                0.,
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroYielddMicroGradientStrainISV) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    DEBUG[ "dMicroGradientYielddMicroGradientStrainISV" ][ 3 * j + i - 47 ],
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual2 (dMacroGradientYielddMicroGradientStrainISV) & False\n";
                        return 1;
                    }
                }
            } 
        }

        if ( ( i >= 50 ) && ( i < 55 ) ){

            /*==============================
            | Elastic Deformation Measures |
            ==============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticDeformationGradientdGammas) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroDeformationdGammas) & False\n";

                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGradientMicroDeformationdGammas) & False\n";

                    return 1;
                }
            }

            /*!=====================================
            | Elastic Derived Deformation Measures |
            ======================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticRightCauchyGreendGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroRightCauchyGreendGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticPsidGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    vectorTools::print( numericGradients[ "currentElasticGamma" ] );
                    results << "test_computePlasticDeformationResidual2 (dElasticGammadGammas) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPK2StressdGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceMicroStressdGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                0.,
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceHigherOrderStressdGammas) & False\n";
                    return 1;
                }
            }

            /*!================
            | Cohesion values |
            =================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMacroCohesiondGammas) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            variableVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMicroCohesiondGammas) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            variableVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dCurrentMicroGradientCohesiondGammas) & False\n";
                return 1;
            }

            /*=================
            | Flow Directions |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroFlowDirectiondGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroFlowDirectiondGammas) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientFlowDirectiondGammas) & False\n";
                    return 1;
                }
            }

            /*!==========================
            | Expected strain-like ISVs |
            ===========================*/

            if ( i == 50 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                                DEBUG[ "dExpectedMacroISVdMacroGamma" ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdMacroGamma) & False\n";
                    return 1;
                }

                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdMacroGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                                variableVector( 3, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdMacroGamma) & False\n";
                    return 1;
                }
            }

            if ( i == 51 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdMicroGamma) & False\n";
                    return 1;
                }

                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                                DEBUG[ "dExpectedMicroISVdMicroGamma" ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdMicroGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                                variableVector( 3, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdMicroGamma) & False\n";
                    return 1;
                }
            }

            if ( i >= 52 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdMicroGradientGamma) & False\n";
                    return 1;
                }

                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                                variableVector( 1, 0 ),
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdMicroGradientGamma) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < 3; j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ][ j ],
                                                    DEBUG[ "dExpectedMicroGradientISVdMicroGradientGamma" ][ 3 * j + i - 52 ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
            }


            /*!===========================
            | Plastic velocity gradients |
            ============================*/

            if ( i == 50 ){

                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMacroVelocityGradientdMacroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdMacroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdMacroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdMacroGamma) & False\n";
                        return 1;
                    }
                }
            }

            if ( i == 51 ){

                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMacroVelocityGradientdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdMicroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMicroVelocityGradientdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdMicroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdMicroGamma) & False\n";
                        return 1;
                    }
                }
            }

            if ( i >= 52 ){

                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroGradientGamma" ][ 3 * j + i - 52 ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dPlasticGradientMicroVelocityGradientdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
            }

            /*!======================================
            | Expected plastic deformation measures |
            =======================================*/

            if ( i == 50 ){

                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                    DEBUG[ "dExpectedPlasticDeformationGradientdMacroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroVelocityGradientdMacroGamma & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdMacroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMacroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdMacroGamma) & False\n";
                        return 1;
                    }
                }
            }

            else if ( i == 51 ){

                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                    DEBUG[ "dExpectedPlasticDeformationGradientdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroVelocityGradientdMicroGamma & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                    DEBUG[ "dExpectedPlasticMicroDeformationdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdMicroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
                        std::cout << "derp\n";
                        std::cout << "numeric:  "; vectorTools::print( numericGradients[ "expectedPlasticGradientMicroDeformation" ] );
                        std::cout << "analytic: "; vectorTools::print( DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGamma" ] );
                        results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdMicroGamma) & False\n";
                        return 1;
                    }
                }
            }

            else if ( i >= 52 ){

                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroVelocityGradientdMicroGradientGamma & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                    0.,
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                    DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroGradientGamma" ][ 3 * j + i - 52 ],
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
            }

            /*!====================
            | The yield equations |
            =====================*/

            if ( i == 50 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
                                                variableVector( 1, 0 ),
                                                1e-5, 1e-8 ) ){
                    std::cout << "numeric:  "; vectorTools::print( numericGradients[ "macroYieldFunction" ] );
                    std::cout << "analytic: "; vectorTools::print( DEBUG[ "dMacroYielddMacroGamma" ] );
                    results << "test_computePlasticDeformationResidual2 (dMacroYielddMacroGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
                                                variableVector( 1, 0 ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroYielddMacroMacroGamma) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual2 (dMacroGradientYielddMacroGamma) & False\n";
                        return 1;
                    }
                }
            }
            else if ( i == 51 ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
                                                variableVector( 1, 0. ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroYielddMicroGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
                                                variableVector( 1, 0. ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroYielddMicroGamma) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual2 (dMacroGradientYielddMicroStrainISV) & False\n";
                        return 1;
                    }
                }
            }
            else if ( i >= 52 ){

                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ],
                                                variableVector( 1, 0. ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroYielddMicroGradientGamma) & False\n";
                    return 1;
                }
    
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ],
                                                variableVector( 1, 0. ),
                                                1e-5, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroYielddMicroGradientGamma) & False\n";
                    return 1;
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                    0.,
                                                    1e-5, 1e-8 ) ){
                        results << "test_computePlasticDeformationResidual2 (dMacroGradientYielddMicroGradientGamma) & False\n";
                        return 1;
                    }
                }
            } 

        }

        if ( i >= 55 ){
            std::cout << "ERROR in the test!\n";
            assert( 1 == 0 );
        }

#endif


        //Jacobian test

        solverTools::floatVector gradCol = ( residual_P - residual_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ], 1e-5, 1e-7 ) ){
                std::cout << "i, j: " << i << ", " << j << "\n";
                std::cout << "num: " << gradCol[ j ] << "\n";
                std::cout << "ana: " << jacobian[ j ][ i ] << "\n";
                std::cout << "error: " << gradCol[ j ] - jacobian[ j ][ i ] << "\n";
                std::cout << "gradCol:\n"; vectorTools::print( gradCol );
                std::cout << "jacobian:\n"; vectorTools::print( jacobian );
                results << "test_computePlasticDeformationResidual2 (test 2) & False\n";
                return 1;
            }
        }

        //Test the partial derivative of the stress measures w.r.t. the plastic fundamental deformation measures
        gradCol = vectorTools::appendVectors( { fO_P[ 0 ] - fO_M[ 0 ], fO_P[ 1 ] - fO_M[ 1 ], fO_P[ 2 ] - fO_M[ 2 ] } ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 3 ][ 55 * j + i ], 1e-5 ) ){
                std::cout << "i, j: " << i << ", " << j << "\n";
                std::cout << "num:   " << gradCol[ j ] << "\n";
                std::cout << "ana:   " << floatOuts[ 3 ][ 55 * j + i ] << "\n";
                std::cout << "error: " << gradCol[ j ] - floatOuts[ 3 ][ 55 * j + i ] << "\n";
                results << "test_computePlasticDeformationResidual2 (test 3) & False\n";
                return 1;
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

        if ( error ){
            error->print();
            results << "test_computePlasticDeformationResidual2 & False\n";
            return 1;
        }
        
        error = micromorphicElastoPlasticity::computePlasticDeformationResidual( x, floatArgs_M, intArgs, residual_M, _J,
                                                                                 fO_M, iO_M
#ifdef DEBUG_MODE
                                                                                 , DEBUG_M
#endif
                                                                               );

        if ( error ){
            error->print();
            results << "test_computePlasticDeformationResidual2 & False\n";
            return 1;
        }

        //Debug each of the sub-jacobians if required. This can be very slow so it isn't done all the time

#ifdef DEBUG_MODE

        //Assemble the numeric Jacobians w.r.t. the deformation parameters
        solverTools::debugMap numericGradients;
        
        for ( auto it_P = DEBUG_P.begin(); it_P != DEBUG_P.end(); it_P++ ){
            
            auto it_M = DEBUG_M.find( it_P->first );
            if ( it_M == DEBUG_M.end() ){
                std::cerr << "ERROR: A KEY EXISTS IN DEBUG_P THAT DOESNT EXIST IN DEBUG_M\n";
                results << "test_computePlasticDeformationResidual2 & False\n";
                return 1;
            }

            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );
        }

        //Debug the Jacobians w.r.t. the current deformation gradient
        if ( i < 9 ) {

            /*==========================
            | The Elastic Deformations |
            ==========================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                DEBUG[ "dElasticDeformationGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticDeformationGradientdDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroDeformationdDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGradientMicroDeformationdDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!=========================================
            | The Elastic Derived Deformation Measures |
            ==========================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticRightCauchyGreendDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticRightCauchyGreendDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroRightCauchyGreendDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticPsidDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGammadDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPK2StressdDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceMicroStressdDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceHigherOrderStressdDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!==============
            | The Cohesions |
            ===============*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dMacroCohesiondDeformationGradient) & False\n";
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dMicroCohesiondDeformationGradient) & False\n";
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            solverTools::floatVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dMicroGradientCohesiondDeformationGradient) & False\n";
            }

            /*!====================
            | The Flow Directions |
            =====================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroFlowDirectiondDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroFlowDirectiondDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondDeformationGradient" ][ 9 * j + i ],
                                                1e-8, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientFlowDirectiondDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!=====================
            | The Strain-Like ISVs |
            ======================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdDeformationGradient) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            solverTools::floatVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdDeformationGradient) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic Velocity Gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroGradientVelocityGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroGradientVelocityGradientdDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!==============================
            | Plastic Deformation Gradients |
            ===============================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticDeformationGradientdDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdDeformationGradient" ][ 9 * j + i ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdDeformationGradient) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The Yield Equations |
            =====================*/

            for ( unsigned int j = 0; j < numericGradients[ "macroYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ j ],
                                                DEBUG[ "dMacroYielddDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroYielddDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "microYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ j ],
                                                DEBUG[ "dMicroYielddDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroYielddDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                DEBUG[ "dMicroGradientYielddDeformationGradient" ][ 9 * j + i ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientYielddDeformationGradient) & False\n";
                    return 1;
                }
            }
        }

        if ( ( i >= 9 ) && ( i < 18 ) ){

            /*==========================
            | The Elastic Deformations |
            ==========================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticDeformationGradientdMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticMicroDeformationdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroDeformationdMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGradientMicroDeformationdMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!=========================================
            | The Elastic Derived Deformation Measures |
            ==========================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticRightCauchyGreendDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                DEBUG[ "dElasticMicroRightCauchyGreendMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroRightCauchyGreendMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                DEBUG[ "dElasticPsidMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticPsidMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGammadMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-9 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    std::cout << "num:   " << numericGradients[ "currentPK2Stress" ][ j ] << "\n";
                    std::cout << "ana:   " << DEBUG[ "dPK2StressdMicroDeformation" ][ 9 * j + i - 9 ] << "\n";
                    std::cout << "error: " << numericGradients[ "currentPK2Stress" ][ j ] - DEBUG[ "dPK2StressdMicroDeformation" ][ 9 * j + i - 9 ] << "\n";
                    results << "test_computePlasticDeformationResidual2 (dPK2StressdMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-9 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    std::cout << "num:   " << numericGradients[ "currentReferenceMicroStress" ][ j ] << "\n";
                    std::cout << "ana:   " << DEBUG[ "dReferenceMicroStressdMicroDeformation" ][ 9 * j + i - 9 ] << "\n";
                    std::cout << "error: " << numericGradients[ "currentReferenceMicroStress" ][ j ] - DEBUG[ "dReferenceMicroStressdMicroDeformation" ][ 9 * j + i - 9 ] << "\n";
                    results << "test_computePlasticDeformationResidual2 (dReferenceMicroStressdMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-9 ) ){
                    std::cout << "i, j: " << i << ", " << j << "\n";
                    std::cout << "num:   " << numericGradients[ "currentReferenceHigherOrderStress" ][ j ] << "\n";
                    std::cout << "ana:   " << DEBUG[ "dReferenceHigherOrderStressdMicroDeformation" ][ 9 * j + i - 9 ] << "\n";
                    std::cout << "error: " << numericGradients[ "currentReferenceHigherOrderStress" ][ j ] - DEBUG[ "dReferenceHigherOrderStressdMicroDeformation" ][ 9 * j + i - 9 ] << "\n";
                    results << "test_computePlasticDeformationResidual2 (dReferenceHigherOrderStressdMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!==============
            | The Cohesions |
            ===============*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dMacroCohesiondMicroDeformation) & False\n";
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dMicroCohesiondMicroDeformation) & False\n";
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            solverTools::floatVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dMicroGradientCohesiondMicroDeformation) & False\n";
            }

            /*!====================
            | The Flow Directions |
            =====================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroFlowDirectiondMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroFlowDirectiondMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientFlowDirectiondMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!=====================
            | The Strain-Like ISVs |
            ======================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            solverTools::floatVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdMicroDeformation) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic Velocity Gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroGradientVelocityGradientdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroGradientVelocityGradientdMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!==============================
            | Plastic Deformation Gradients |
            ===============================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticDeformationGradientdMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The Yield Equations |
            =====================*/

            for ( unsigned int j = 0; j < numericGradients[ "macroYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ j ],
                                                DEBUG[ "dMacroYielddMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroYielddMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "microYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ j ],
                                                DEBUG[ "dMicroYielddMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-8 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroYielddMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                DEBUG[ "dMicroGradientYielddMicroDeformation" ][ 9 * j + i - 9 ],
                                                1e-6, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientYielddMicroDeformation) & False\n";
                    return 1;
                }
            }
        }

        if ( ( i >= 18 ) && ( i < 45 ) ){

            /*==========================
            | The Elastic Deformations |
            ==========================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticDeformationGradientdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroDeformationdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dElasticGradientMicroDeformationdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGradientMicroDeformationdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!=========================================
            | The Elastic Derived Deformation Measures |
            ==========================================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticRightCauchyGreendGradientDeformationGradient) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroRightCauchyGreen" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroRightCauchyGreen" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticMicroRightCauchyGreendGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticPsi" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticPsi" ][ j ],
                                                0.,
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticPsidGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentElasticGamma" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticGamma" ][ j ],
                                                DEBUG[ "dElasticGammadGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dElasticGammadGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!================
            | Stress Measures |
            =================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPK2Stress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPK2Stress" ][ j ],
                                                DEBUG[ "dPK2StressdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPK2StressdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceMicroStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceMicroStress" ][ j ],
                                                DEBUG[ "dReferenceMicroStressdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceMicroStressdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentReferenceHigherOrderStress" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentReferenceHigherOrderStress" ][ j ],
                                                DEBUG[ "dReferenceHigherOrderStressdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-5, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dReferenceHigherOrderStressdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!==============
            | The Cohesions |
            ===============*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dMacroCohesiondGradientMicroDeformation) & False\n";
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroCohesion" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dMicroCohesiondGradientMicroDeformation) & False\n";
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientCohesion" ],
                                            solverTools::floatVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dMicroGradientCohesiondGradientMicroDeformation) & False\n";
            }

            /*!====================
            | The Flow Directions |
            =====================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
                                                DEBUG[ "dMacroFlowDirectiondGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroFlowDirectiondGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
                                                DEBUG[ "dMicroFlowDirectiondGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroFlowDirectiondGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
                                                DEBUG[ "dMicroGradientFlowDirectiondGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-8, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientFlowDirectiondGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!=====================
            | The Strain-Like ISVs |
            ======================*/

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMacroStrainISV" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMacroStrainISVdGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroStrainISV" ],
                                            solverTools::floatVector( 1, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroStrainISVdGradientMicroDeformation) & False\n";
                return 1;
            }

            if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedMicroGradientStrainISV" ],
                                            solverTools::floatVector( 3, 0 ),
                                            1e-9, 1e-9 ) ){
                results << "test_computePlasticDeformationResidual2 (dExpectedMicroGradientStrainISVdGradientMicroDeformation) & False\n";
                return 1;
            }

            /*!===========================
            | Plastic Velocity Gradients |
            ============================*/

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMacroVelocityGradientdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMacroVelocityGradientdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroVelocityGradientdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroVelocityGradientdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
                                                DEBUG[ "dPlasticMicroGradientVelocityGradientdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dPlasticMicroGradientVelocityGradientdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!==============================
            | Plastic Deformation Gradients |
            ===============================*/

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticDeformationGradient" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticDeformationGradient" ][ j ],
                                                DEBUG[ "dExpectedPlasticDeformationGradientdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticDeformationGradientdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticMicroDeformationdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticMicroDeformationdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "expectedPlasticGradientMicroDeformation" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "expectedPlasticGradientMicroDeformation" ][ j ],
                                                DEBUG[ "dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-9, 1e-9 ) ){
                    results << "test_computePlasticDeformationResidual2 (dExpectedPlasticGradientMicroDeformationdGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            /*!====================
            | The Yield Equations |
            =====================*/

            for ( unsigned int j = 0; j < numericGradients[ "macroYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "macroYieldFunction" ][ j ],
                                                DEBUG[ "dMacroYielddGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMacroYielddGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "microYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microYieldFunction" ][ j ],
                                                DEBUG[ "dMicroYielddGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroYielddGradientMicroDeformation) & False\n";
                    return 1;
                }
            }

            for ( unsigned int j = 0; j < numericGradients[ "microGradientYieldFunction" ].size(); j++ ){
                if ( !vectorTools::fuzzyEquals( numericGradients[ "microGradientYieldFunction" ][ j ],
                                                DEBUG[ "dMicroGradientYielddGradientMicroDeformation" ][ 27 * j + i - 18 ],
                                                1e-6, 1e-7 ) ){
                    results << "test_computePlasticDeformationResidual2 (dMicroGradientYielddGradientMicroDeformation) & False\n";
                    return 1;
                }
            }
        }

        if ( i >= 45 ){
            assert( 1 == 0 );
        }

#endif

        //Test the partial derivative of the residual w.r.t. the fundamental deformation measures 
        solverTools::floatVector gradCol = ( residual_P - residual_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 5 ][ 45 * j + i ] ) ){
                results << "test_computePlasticDeformationResidual2 (test 4) & False\n";
                return 1;
            }
        }

        //Test the partial derivative of the stress measures w.r.t. the fundamental deformation measures
        gradCol = vectorTools::appendVectors( { fO_P[ 0 ] - fO_M[ 0 ], fO_P[ 1 ] - fO_M[ 1 ], fO_P[ 2 ] - fO_M[ 2 ] } ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], floatOuts[ 4 ][ 45 * j + i ], 1e-5 ) ){
                results << "test_computePlasticDeformationResidual2 (test 5) & False\n";
                return 1;
            }
        }
    }

    results << "test_computePlasticDeformationResidual2 & True\n";
    return 0;
}

int test_materialLibraryInterface( std::ofstream &results ){
    /*!
     * Test the interface to the linear elastic model
     * via the material library.
     *
     * :param std::ofstream &results: The output file.
     */

    //Initialize the model
    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    //Set up the inputs
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
    solverTools::homotopyMap DEBUG;
#endif

    solverTools::floatVector PK2_answer = { 177.067  , 13.287 ,  -0.577489,
                                             10.7222 , 152.621,  -0.288668,
                                             -1.38898, 1.61632, 150.85 };
    solverTools::floatVector SIGMA_answer = { 178.961  ,  13.5623 ,  -2.43027,
                                               13.5623 , 151.785  ,   1.62465,
                                               -2.43027,   1.62465, 149.31 };
    solverTools::floatVector M_answer = { 0.541081, -0.533639,  0.640843,  2.92886 ,  1.1505  ,
                                          1.14946 ,  0.605546, -2.62366 ,  1.54942 , -2.40585 ,
                                         -0.670181, -0.848562,  0.678723,  0.433934, -0.206886,
                                         -2.7026  ,  0.964407,  1.68227 , -0.480138,  2.68016 ,
                                         -0.628373,  1.14616 , -0.10665 , -2.2393  , -0.765349,
                                          0.746722,  0.994415 };

    solverTools::floatVector SDVS_answer = { -0.0394628  ,  0.102019   ,  0.0453858  ,  0.0262009 ,  0.0232167 ,
                                             -0.0196885  ,  0.0357449  ,  0.0453858  ,  0.0262009 ,  0.0232167 ,
                                              0.00813011 ,  0.00956498 , -0.000691731,  0.00885806, -0.0103998 ,
                                              0.000451969, -0.00115188 ,  0.00088041 , -0.00882951,  0.0622934 ,
                                              0.0274584  , -0.00113356 ,  0.0403186  , -0.00945181,  0.00119733,
                                             -0.00227071 ,  0.000880704, -0.00792983 ,  0.0554737 , -0.0202861 ,
                                             -0.0208186  ,  0.0122486  ,  0.0119279  ,  0.037049  ,  0.0168925 ,
                                             -0.0430597  , -0.0242213  ,  0.0244271  , -0.00880588,  0.0383104 ,
                                              0.00423276 ,  0.0111449  ,  0.0155982  ,  0.0136987 , -0.00813611,
                                             -0.0305008  ,  0.0171903  , -0.0126259  ,  0.0500443 , -0.0276203 ,
                                              0.0193062  ,  0.0310543  ,  0.0157931  ,  0.0168109 ,  0.0196658 };

    std::vector< double > SDVS = SDVSDefault;

    std::vector< double > PK2_result, SIGMA_result, M_result;

    //Evaluate the model
    int errorCode = material->evaluate_model( time, fparams,
                                              current_grad_u, current_phi, current_grad_phi,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              SDVS,
                                              current_ADD_DOF, current_ADD_grad_DOF,
                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                              PK2_result, SIGMA_result, M_result,
                                              ADD_TERMS,
                                              output_message
#ifdef DEBUG_MODE
                                              , DEBUG
#endif
                                            );

    if ( errorCode > 0 ){
        std::cout << output_message << "\n";
        results << "test_materialLibraryInterface & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) ){
        std::cout << "PK2_result:\n"; vectorTools::print( PK2_result );
        std::cout << "PK2_answer:\n"; vectorTools::print( PK2_answer );
        results << "test_materialLibraryInterface (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) ){
        results << "test_materialLibraryInterface (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) ){
        results << "test_materialLibraryInterface (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SDVS, SDVS_answer ) ){
        results << "test_materialLibraryInterface (test 4) & False\n";
        return 1;
    }

    //Check the Jacobian using the previously tested jacobian
    std::vector< std::vector< double > > DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer,
                                         DSIGMADgrad_u_answer, DSIGMADphi_answer, DSIGMADgrad_phi_answer,
                                         DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer;

    std::vector< std::vector< std::vector< double > > > ADD_JACOBIANS;

    SDVS = SDVSDefault;

#ifdef DEBUG_MODE
    DEBUG.clear();
#endif

    errorCode = micromorphicElastoPlasticity::evaluate_model(
                                time, fparams,
                                current_grad_u, current_phi, current_grad_phi,
                                previous_grad_u, previous_phi, previous_grad_phi,
                                SDVS,
                                current_ADD_DOF, current_ADD_grad_DOF,
                                previous_ADD_DOF, previous_ADD_grad_DOF,
                                PK2_result, SIGMA_result, M_result,
                                DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer,
                                DSIGMADgrad_u_answer, DSIGMADphi_answer, DSIGMADgrad_phi_answer,
                                DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer,
                                ADD_TERMS, ADD_JACOBIANS,
                                output_message
#ifdef DEBUG_MODE
                                , DEBUG
#endif
                              );

    if ( errorCode > 0 ){
        std::cout << output_message << "\n";
        results << "test_materialLibraryInterface & False\n";
        return 1;
    }

    PK2_result.clear();
    SIGMA_result.clear();
    M_result.clear();

    SDVS = SDVSDefault;

    std::vector< std::vector< double > > DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
                                         DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
                                         DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result;

    errorCode = material->evaluate_model( time, fparams,
                                          current_grad_u, current_phi, current_grad_phi,
                                          previous_grad_u, previous_phi, previous_grad_phi,
                                          SDVS,
                                          current_ADD_DOF, current_ADD_grad_DOF,
                                          previous_ADD_DOF, previous_ADD_grad_DOF,
                                          PK2_result, SIGMA_result, M_result,
                                          DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
                                          DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
                                          DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result,
                                          ADD_TERMS, ADD_JACOBIANS,
                                          output_message
#ifdef DEBUG_MODE
                                          , DEBUG
#endif
                                        );

    if ( !vectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) ){
        results << "test_materialLibraryInterface (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) ){
        results << "test_materialLibraryInterface (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) ){
        results << "test_materialLibraryInterface (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SDVS_answer, SDVS ) ){
        results << "test_materialLibraryInterface (test 8) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DPK2Dgrad_u_result, DPK2Dgrad_u_answer ) ){
        results << "test_materialLibraryInterface (test 9) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DPK2Dphi_result, DPK2Dphi_answer ) ){
        results << "test_materialLibraryInterface (test 10) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer ) ){
        results << "test_materialLibraryInterface (test 11) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DSIGMADgrad_u_result, DSIGMADgrad_u_answer ) ){
        results << "test_materialLibraryInterface (test 12) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DSIGMADphi_result, DSIGMADphi_answer ) ){
        results << "test_materialLibraryInterface (test 13) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer ) ){
        results << "test_materialLibraryInterface (test 14) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DMDgrad_u_result, DMDgrad_u_answer ) ){
        results << "test_materialLibraryInterface (test 15) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DMDphi_result, DMDphi_answer ) ){
        results << "test_materialLibraryInterface (test 16) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DMDgrad_phi_result, DMDgrad_phi_answer ) ){
        results << "test_materialLibraryInterface (test 17) & False\n";
        return 1;
    }

#ifdef DEBUG_MODE
    DEBUG.clear();
#endif

    //Test the computed numeric Jacobian values
    SDVS = SDVSDefault;
    errorCode = material->evaluate_model_numeric_gradients( time, fparams,
                                                            current_grad_u, current_phi, current_grad_phi,
                                                            previous_grad_u, previous_phi, previous_grad_phi,
                                                            SDVS,
                                                            current_ADD_DOF, current_ADD_grad_DOF,
                                                            previous_ADD_DOF, previous_ADD_grad_DOF,
                                                            PK2_result, SIGMA_result, M_result,
                                                            DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
                                                            DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
                                                            DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result,
                                                            ADD_TERMS, ADD_JACOBIANS,
                                                            output_message,
#ifdef DEBUG_MODE
                                                            DEBUG,
#endif
                                                            1e-6 );

    if ( errorCode > 0 ){
        results << "test_materialLibraryInterface & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) ){
        results << "test_materialLibraryInterface (test 18) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) ){
        results << "test_materialLibraryInterface (test 19) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) ){
        results << "test_materialLibraryInterface (test 20) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SDVS, SDVS_answer ) ){
        std::cout << "error: "; vectorTools::print( SDVS - SDVS_answer );
        results << "test_materialLibraryInterface (test 21) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DPK2Dgrad_u_result, DPK2Dgrad_u_answer ) ){
        std::cout << "num:\n"; vectorTools::print( DPK2Dgrad_u_result );
        std::cout << "ana:\n"; vectorTools::print( DPK2Dgrad_u_answer );
        results << "test_materialLibraryInterface (test 22) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DPK2Dphi_result, DPK2Dphi_answer ) ){
        results << "test_materialLibraryInterface (test 23) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer ) ){
        results << "test_materialLibraryInterface (test 24) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DSIGMADgrad_u_result, DSIGMADgrad_u_answer ) ){
        results << "test_materialLibraryInterface (test 25) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DSIGMADphi_result, DSIGMADphi_answer ) ){
        results << "test_materialLibraryInterface (test 26) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer ) ){
        results << "test_materialLibraryInterface (test 27) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DMDgrad_u_result, DMDgrad_u_answer ) ){
        results << "test_materialLibraryInterface (test 28) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DMDphi_result, DMDphi_answer ) ){
        results << "test_materialLibraryInterface (test 29) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( DMDgrad_phi_result, DMDgrad_phi_answer ) ){
        results << "test_materialLibraryInterface (test 30) & False\n";
        return 1;
    }

    results << "test_materialLibraryInterface & True\n";
    return 1;
}

int test_materialLibraryInterface2( std::ofstream &results ){
    /*!
     * Test the interface to the linear elastic model
     * via the material library.
     *
     * NOTE: This function mostly exists to perform debugging
     *       on the implementation of the function into a 
     *       larger solver code.
     *
     * :param std::ofstream &results: The output file.
     */

    //Initialize the model
    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    //Set up the inputs
    //Initialize the time
    std::vector< double > time = { 0.045, 0.01 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 170, 15, 2, 140, 20, 2, 2, 27, 2, 0.56, 0.2, 2, 0.15, 0.3, 2, 0.82, 0.1, 2, 0.42, 0.3, 2, 0.05, 0.2, 2, 0.52, 0.4, 2, 29480, 25480, 5, 1000, 400, -1500, -1400, -3000, 11, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0, 2, 400, -3000, 0.5, 0.5, 0.5, 1e-09, 1e-09 };

    //Initialize the gradient of the macro displacement
    double current_grad_u[ 3 ][ 3 ] =
    {
        { -0.00124343, -6.55319e-14, 3.99657e-13},
        { 0, 0.0045, 0},
        { -1.75135e-13, -1.35481e-13, -0.00124343 },
    };

    double previous_grad_u[ 3 ][ 3 ] =
    {
        { -0.00123858, -1.22379e-17, 5.04154e-18},
        { 0, 0.004, 0},
        { -1.47723e-18, 4.44523e-18, -0.00123858 },
    };

    //Initialize the micro displacement
    double current_phi[ 9 ] = { -0.00153489, -3.04626e-13, 5.16537e-13, 1.58771e-13, 0.00303407, 4.29828e-14, -4.38368e-13, -1.80694e-13, -0.00153489 };

    double previous_phi[ 9 ] = { -0.00164749, -2.63663e-17, 1.35603e-17, 8.65138e-19, 0.00325613, -2.13082e-20, -1.17433e-17, 2.24626e-18, -0.00164749 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };

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
    solverTools::homotopyMap DEBUG;
#endif

    std::vector< double > SDVS = SDVSDefault;

    std::vector< double > PK2_result, SIGMA_result, M_result;

    //Evaluate the model
    int errorCode = material->evaluate_model( time, fparams,
                                              current_grad_u, current_phi, current_grad_phi,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              SDVS,
                                              current_ADD_DOF, current_ADD_grad_DOF,
                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                              PK2_result, SIGMA_result, M_result,
                                              ADD_TERMS,
                                              output_message
#ifdef DEBUG_MODE
                                              , DEBUG
#endif
                                            );

    if ( errorCode > 0 ){
        std::cout << output_message << "\n";
        results << "test_materialLibraryInterface & False\n";
        return 1;
    }


    std::cout << "SDVS:\n"; vectorTools::print( SDVS );
    std::cout << "PK2:\n"; vectorTools::print( PK2_result );
    std::cout << "SIGMA:\n"; vectorTools::print( SIGMA_result );
    std::cout << "M:\n"; vectorTools::print( M_result );

#ifdef DEBUG_MODE
    for ( auto step = DEBUG.begin(); step != DEBUG.end(); step++ ){
        std::cout << step->first << "\n";
        for ( auto iteration = step->second.begin(); iteration != step->second.end(); iteration++ ){
            std::cout << "    " << iteration->first << "\n";
            for ( auto value = iteration->second.begin(); value != iteration->second.end(); value++ ){
                if ( value->second.size() <= 27 ) {
                    std::cout << "        " << value->first << "\n";
                    std::cout << "            "; vectorTools::print( value->second );
                }
            }
        }
    }
#endif

    return 0;
}

int test_evaluate_model_continuation( std::ofstream &results){
    /*!
     * Test the evaluation of the constitutive model when being
     * continued from a previous plastic deformation
     *
     * :param std::ofstream &results: The output file.
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

    solverTools::floatVector PK2Answer0 = { 177.067  , 13.287 ,  -0.577489,
                                             10.7222 , 152.621,  -0.288668,
                                             -1.38898, 1.61632, 150.85 };
    solverTools::floatVector SigmaAnswer0 = { 178.961  ,  13.5623 ,  -2.43027,
                                              13.5623 , 151.785  ,   1.62465,
                                              -2.43027,   1.62465, 149.31 };
    solverTools::floatVector MAnswer0 = { 0.541081, -0.533639,  0.640843,  2.92886 ,  1.1505  ,
                                          1.14946 ,  0.605546, -2.62366 ,  1.54942 , -2.40585 ,
                                         -0.670181, -0.848562,  0.678723,  0.433934, -0.206886,
                                         -2.7026  ,  0.964407,  1.68227 , -0.480138,  2.68016 ,
                                         -0.628373,  1.14616 , -0.10665 , -2.2393  , -0.765349,
                                          0.746722,  0.994415 };

    solverTools::floatVector SDVSAnswer0 = { -0.0394628  ,  0.102019   ,  0.0453858  ,  0.0262009 ,  0.0232167 ,
                                             -0.0196885  ,  0.0357449  ,  0.0453858  ,  0.0262009 ,  0.0232167 ,
                                              0.00813011 ,  0.00956498 , -0.000691731,  0.00885806, -0.0103998 ,
                                              0.000451969, -0.00115188 ,  0.00088041 , -0.00882951,  0.0622934 ,
                                              0.0274584  , -0.00113356 ,  0.0403186  , -0.00945181,  0.00119733,
                                             -0.00227071 ,  0.000880704, -0.00792983 ,  0.0554737 , -0.0202861 ,
                                             -0.0208186  ,  0.0122486  ,  0.0119279  ,  0.037049  ,  0.0168925 ,
                                             -0.0430597  , -0.0242213  ,  0.0244271  , -0.00880588,  0.0383104 ,
                                              0.00423276 ,  0.0111449  ,  0.0155982  ,  0.0136987 , -0.00813611,
                                             -0.0305008  ,  0.0171903  , -0.0126259  ,  0.0500443 , -0.0276203 ,
                                              0.0193062  ,  0.0310543  ,  0.0157931  ,  0.0168109 ,  0.0196658 };

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

    if ( errorCode != 0 ){
        std::cout << output_message;
        results << "test_evaluate_model_continuation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SDVS, SDVSAnswer0 ) ){
        results << "test_evaluate_model_continuation (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( PK2Answer0, current_PK20, 1e-5, 1e-5 ) ){
        std::cout << "error: "; vectorTools::print( PK2Answer0 - current_PK20 );
        results << "test_evaluate_model_continuation (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( SigmaAnswer0, current_SIGMA0, 1e-5, 1e-5 ) ){
        results << "test_evaluate_model_continuation (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( MAnswer0, current_M0, 1e-5 ) ){
        results << "test_evaluate_model_continuation (test 4) & False\n";
        return 1;
    }

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

    std::vector< double  > SDVSAnswer1 = { -0.0657714, 0.145741, 0, 0, 0,
                                            0, 0, 0, 0, 0,
                                            0.0135502 , 0.0159416  , -0.00115289 ,
                                            0.0147634 , -0.017333  ,  0.000753281,
                                           -0.0019198 ,  0.00146735, -0.0147159  ,
                                            0.0889905 ,  0.0392263 , -0.00161937 ,
                                            0.057598  , -0.0135026 ,  0.00171046 ,
                                           -0.00324386, 0.00125815 , -0.0113283  ,
                                            0.0853441 , -0.0312094 , -0.0320287  ,
                                            0.018844  , 0.0183506  ,  0.0569985  ,
                                            0.0259884 , -0.0662456 , -0.0372635  ,
                                            0.0375801 , -0.0135475 ,  0.0589391  ,
                                            0.00651194, 0.0171461  ,  0.0239973  ,
                                            0.0210749 , -0.0125171 , -0.0469244  ,
                                            0.0264467 , -0.0194244 ,  0.0769913  ,
                                           -0.0424927 , 0.0297019  ,  0.0477759  ,
                                            0.0242971 , 0.0258629  ,  0.030255 };

    std::vector< double > PK2Answer1 = { 1.74594528e+02,  9.90357506e+00, -1.16456128e-01,
                                         7.35100585e+00,  1.57060509e+02, -4.68051721e-01,
                                        -8.40302786e-01,  1.28400640e+00,  1.55106125e+02 };

    std::vector< double > SIGMAAnswer1 = { 175.50856348,  10.16163404,  -1.73547912,
                                            10.16163404, 155.38685738,   1.3379937 ,
                                            -1.73547912,   1.3379937 , 152.87515001 };

    std::vector< double > MAnswer1 = { 0.37220031, -0.51491634,  0.60588705,  2.5117993 ,  1.15120688,
                                       1.01073726,  0.54881061, -2.49187357,  1.42408148, -2.15103266,
                                      -0.81552364, -0.72425736,  0.6621303 ,  0.29584242, -0.22199834,
                                      -2.4997514 ,  1.02467573,  1.50571566, -0.43582018,  2.57530701,
                                      -0.56574275,  0.94510996, -0.17195943, -2.08860752, -0.83703857,
                                       0.6563958 ,  0.94868265 };

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

    if ( !vectorTools::fuzzyEquals( SDVS, SDVSAnswer1 ) ){
        results << "test_evaluate_model_continuation (test 9) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( current_PK21, PK2Answer1 ) ){
        std::cout << "error: "; vectorTools::print( current_PK21 - PK2Answer1 );
        results << "test_evaluate_model_continuation (test 10) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( current_SIGMA1, SIGMAAnswer1 ) ){
        std::cout << "error: "; vectorTools::print( current_SIGMA1 - SIGMAAnswer1 );
        results << "test_evaluate_model_continuation (test 11) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( current_M1, MAnswer1 ) ){
        std::cout << "error: "; vectorTools::print( current_M1 - MAnswer1 );
        results << "test_evaluate_model_continuation (test 12) & False\n";
        return 1;
    }

#ifdef DEBUG_MODE

    if ( !vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPK2Stress" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "intermediatePK2Stress" ] ) ){
        results << "test_evaluate_model_continuation (intermediatePK2Stress) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousReferenceMicroStress" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "intermediateReferenceMicroStress" ] ) ){
        results << "test_evaluate_model_continuation (intermediateReferenceMicroStress) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousReferenceHigherOrderStress" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "intermediateReferenceHigherOrderStress" ] ) ){
        results << "test_evaluate_model_continuation (intermediateReferenceHigherOrderStress) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPlasticDeformationGradient" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "convergedPlasticDeformationGradient" ] ) ){
        results << "test_evaluate_model_continuation (previousElasticDeformationGradient) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPlasticMicroDeformation" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "convergedPlasticMicroDeformation" ] ) ){
        results << "test_evaluate_model_continuation (previousPlasticMicroDeformation) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( homotopyDEBUG1[ "pre_iteration_values" ][ "pre_iteration_values" ][ "previousPlasticGradientMicroDeformation" ],
                                    homotopyDEBUG0[ "converged_values" ][ "converged_values" ][ "convergedPlasticGradientMicroDeformation" ] ) ){
        results << "test_evaluate_model_continuation (previousPlasticGradientMicroDeformation) & False\n";
        return 1;
    }

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

    if ( !vectorTools::fuzzyEquals( SDVS, SDVS1 ) ){
        results << "test_evaluate_model_continuation (test 9) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( current_PK21, current_PK22 ) ){
        results << "test_evaluate_model_continuation (test 10) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( current_SIGMA1, current_SIGMA2 ) ){
        results << "test_evaluate_model_continuation (test 11) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( current_M1, current_M2 ) ){
        results << "test_evaluate_model_continuation (test 12) & False\n";
        return 1;
    }
     

    results << "test_evaluate_model_continuation & True\n";
    return 0;
}

int test_evaluate_model_history( std::ofstream &results ){
    /*!
     * Test the material model undergoing a time history.
     *
     * :param std::ofstream &results: The output file.
     */

#ifdef DEBUG_MODE
    std::ofstream output_file;
    output_file.open( "output_file.txt" );
#endif

    //Initialize the model
    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    std::vector< std::vector< double > > grad_u_0 = { { 0, 0, 0 },
                                                      { 0, 0, 0 },
                                                      { 0, 0, 0 } };

    std::vector< std::vector< double > > grad_u_f = { { 0.5, 0, 0 },
                                                      { 0.0, 0, 0 },
                                                      { 0.0, 0, 0 } };

    std::vector< double > phi_0 = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector< double > phi_f = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    std::vector< std::vector< double > > grad_phi_0 = { { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0} };

    std::vector< std::vector< double > > grad_phi_f = { { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0} };

    double dt = 0.01;
    double t0 = 0.;
    double tf = 1.0;

    double t = t0;

    //Set up the model parameters
    std::vector< double > fparams = { 2, 1e3, 1e2,
                                      2, 1e3, 1e4,
                                      2, 1e3, 1e4,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 29480, 25480,
                                      5, 1000, 400, -1500, -1400, -3000,
                                      11, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0,
                                      2, 400, -3000,
                                      0.5, 0.5, 0.5, 1e-09, 1e-09 };
//                                      0.0, 0.0, 0.0, 1e-09, 1e-09 };

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
    solverTools::homotopyMap DEBUG;
#endif

    std::vector< double > SDVS = SDVSDefault;

    std::vector< double > PK2_result, SIGMA_result, M_result;

    std::vector< std::vector< double > > grad_u_prev   = grad_u_0;
    std::vector< double > phi_prev                     = phi_0;
    std::vector< std::vector< double > > grad_phi_prev = grad_phi_0;

    std::vector< std::vector< double > > grad_u_curr;
    std::vector< double > phi_curr;
    std::vector< std::vector< double > > grad_phi_curr;

    std::vector< double > time;

    double current_grad_u[ 3 ][ 3 ], current_phi[ 9 ], current_grad_phi[ 9 ][ 3 ];
    double previous_grad_u[ 3 ][ 3 ], previous_phi[ 9 ], previous_grad_phi[ 9 ][ 3 ];

    //Initial state
    //Update the arrays
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            previous_grad_u[ i ][ j ] = grad_u_prev[ i ][ j ];
        }
    }

    for ( unsigned int i = 0; i < 9; i++ ){
        previous_phi[ i ] = phi_prev[ i ];

        for ( unsigned int j = 0; j < 3; j++ ){
            previous_grad_phi[ i ][ j ] = grad_phi_prev[ i ][ j ];
        }
    }
    //Evaluate the model
    time = { 0., 0. };
    int errorCode = material->evaluate_model( time, fparams,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              SDVS,
                                              current_ADD_DOF, current_ADD_grad_DOF,
                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                              PK2_result, SIGMA_result, M_result,
                                              ADD_TERMS,
                                              output_message
#ifdef DEBUG_MODE
                                              , DEBUG
#endif
                                            );
    
    if ( errorCode > 0 ){
        std::cout << output_message << "\n";
        results << "test_evaluate_model_history & False\n";
        output_file.close();
        return 1;
    }

#ifdef DEBUG_MODE
    output_file << "NEW_INCREMENT\n";

    //Output the time
    output_file << time[ 0 ] << ", " << time[ 1 ] << "\n";

    //Output the current gradient of u
    output_file << current_grad_u[ 0 ][ 0 ] << ", " << current_grad_u[ 0 ][ 1 ] << ", " << current_grad_u[ 0 ][ 2 ] << ", "
                << current_grad_u[ 1 ][ 0 ] << ", " << current_grad_u[ 1 ][ 1 ] << ", " << current_grad_u[ 1 ][ 2 ] << ", "
                << current_grad_u[ 2 ][ 0 ] << ", " << current_grad_u[ 2 ][ 1 ] << ", " << current_grad_u[ 2 ][ 2 ] << "\n";

    //Output the current micro displacement
    output_file << current_phi[ 0 ] << ", " << current_phi[ 1 ] << ", " << current_phi[ 2 ] << ", "
                << current_phi[ 3 ] << ", " << current_phi[ 4 ] << ", " << current_phi[ 5 ] << ", "
                << current_phi[ 6 ] << ", " << current_phi[ 7 ] << ", " << current_phi[ 8 ] << "\n";

    //Output the current gradient of the micro displacement
    output_file << current_grad_phi[ 0 ][ 0 ] << ", " << current_grad_phi[ 0 ][ 1 ] << ", " << current_grad_phi[ 0 ][ 2 ] << ", "
                << current_grad_phi[ 1 ][ 0 ] << ", " << current_grad_phi[ 1 ][ 1 ] << ", " << current_grad_phi[ 1 ][ 2 ] << ", "
                << current_grad_phi[ 2 ][ 0 ] << ", " << current_grad_phi[ 2 ][ 1 ] << ", " << current_grad_phi[ 2 ][ 2 ] << ", "
                << current_grad_phi[ 3 ][ 0 ] << ", " << current_grad_phi[ 3 ][ 1 ] << ", " << current_grad_phi[ 3 ][ 2 ] << ", "
                << current_grad_phi[ 4 ][ 0 ] << ", " << current_grad_phi[ 4 ][ 1 ] << ", " << current_grad_phi[ 4 ][ 2 ] << ", "
                << current_grad_phi[ 5 ][ 0 ] << ", " << current_grad_phi[ 5 ][ 1 ] << ", " << current_grad_phi[ 5 ][ 2 ] << ", "
                << current_grad_phi[ 6 ][ 0 ] << ", " << current_grad_phi[ 6 ][ 1 ] << ", " << current_grad_phi[ 6 ][ 2 ] << ", "
                << current_grad_phi[ 7 ][ 0 ] << ", " << current_grad_phi[ 7 ][ 1 ] << ", " << current_grad_phi[ 7 ][ 2 ] << ", "
                << current_grad_phi[ 8 ][ 0 ] << ", " << current_grad_phi[ 8 ][ 1 ] << ", " << current_grad_phi[ 8 ][ 2 ] << "\n";

    //Output the PK2 stress
    output_file << PK2_result[ 0 ] << ", " << PK2_result[ 1 ] << ", " << PK2_result[ 2 ] << ", "
                << PK2_result[ 3 ] << ", " << PK2_result[ 4 ] << ", " << PK2_result[ 5 ] << ", "
                << PK2_result[ 6 ] << ", " << PK2_result[ 7 ] << ", " << PK2_result[ 8 ] << "\n";

    //Output the SIGMA stress
    output_file << SIGMA_result[ 0 ] << ", " << SIGMA_result[ 1 ] << ", " << SIGMA_result[ 2 ] << ", "
                << SIGMA_result[ 3 ] << ", " << SIGMA_result[ 4 ] << ", " << SIGMA_result[ 5 ] << ", "
                << SIGMA_result[ 6 ] << ", " << SIGMA_result[ 7 ] << ", " << SIGMA_result[ 8 ] << "\n";

    //Output the M stress
    output_file << M_result[  0 ] << ", " << M_result[  1 ] << ", " << M_result[  2 ] << ", "
                << M_result[  3 ] << ", " << M_result[  4 ] << ", " << M_result[  5 ] << ", "
                << M_result[  6 ] << ", " << M_result[  7 ] << ", " << M_result[  8 ] << ", "
                << M_result[  9 ] << ", " << M_result[ 10 ] << ", " << M_result[ 11 ] << ", "
                << M_result[ 12 ] << ", " << M_result[ 13 ] << ", " << M_result[ 14 ] << ", "
                << M_result[ 15 ] << ", " << M_result[ 16 ] << ", " << M_result[ 17 ] << ", "
                << M_result[ 18 ] << ", " << M_result[ 19 ] << ", " << M_result[ 20 ] << ", "
                << M_result[ 21 ] << ", " << M_result[ 22 ] << ", " << M_result[ 23 ] << ", "
                << M_result[ 24 ] << ", " << M_result[ 25 ] << ", " << M_result[ 26 ] << "\n";

    std::vector< double > PK2_intermediate, SIGMA_intermediate, M_intermediate;

    //Determine if there were non-linear iterations
    auto inc = DEBUG.find( "converged_values" );
    if ( inc != DEBUG.end() ){
        PK2_intermediate   = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediatePK2Stress" ];
        SIGMA_intermediate = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceMicroStress" ];
        M_intermediate     = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceHigherOrderStress" ];
    }
    else{
        PK2_intermediate   = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediatePK2Stress" ];
        SIGMA_intermediate = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceMicroStress" ];
        M_intermediate     = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceHigherOrderStress" ];
    }

    //Output the intermediate PK2 stress
    output_file << PK2_intermediate[ 0 ] << ", " << PK2_intermediate[ 1 ] << ", " << PK2_intermediate[ 2 ] << ", "
                << PK2_intermediate[ 3 ] << ", " << PK2_intermediate[ 4 ] << ", " << PK2_intermediate[ 5 ] << ", "
                << PK2_intermediate[ 6 ] << ", " << PK2_intermediate[ 7 ] << ", " << PK2_intermediate[ 8 ] << "\n";

    //Output the intermediate SIGMA stress
    output_file << SIGMA_intermediate[ 0 ] << ", " << SIGMA_intermediate[ 1 ] << ", " << SIGMA_intermediate[ 2 ] << ", "
                << SIGMA_intermediate[ 3 ] << ", " << SIGMA_intermediate[ 4 ] << ", " << SIGMA_intermediate[ 5 ] << ", "
                << SIGMA_intermediate[ 6 ] << ", " << SIGMA_intermediate[ 7 ] << ", " << SIGMA_intermediate[ 8 ] << "\n";

    //Output the intermediate M stress
    output_file << M_intermediate[  0 ] << ", " << M_intermediate[  1 ] << ", " << M_intermediate[  2 ] << ", "
                << M_intermediate[  3 ] << ", " << M_intermediate[  4 ] << ", " << M_intermediate[  5 ] << ", "
                << M_intermediate[  6 ] << ", " << M_intermediate[  7 ] << ", " << M_intermediate[  8 ] << ", "
                << M_intermediate[  9 ] << ", " << M_intermediate[ 10 ] << ", " << M_intermediate[ 11 ] << ", "
                << M_intermediate[ 12 ] << ", " << M_intermediate[ 13 ] << ", " << M_intermediate[ 14 ] << ", "
                << M_intermediate[ 15 ] << ", " << M_intermediate[ 16 ] << ", " << M_intermediate[ 17 ] << ", "
                << M_intermediate[ 18 ] << ", " << M_intermediate[ 19 ] << ", " << M_intermediate[ 20 ] << ", "
                << M_intermediate[ 21 ] << ", " << M_intermediate[ 22 ] << ", " << M_intermediate[ 23 ] << ", "
                << M_intermediate[ 24 ] << ", " << M_intermediate[ 25 ] << ", " << M_intermediate[ 26 ] << "\n";

    //Output the state variables
    for ( unsigned int i = 0; i < SDVS.size()-1; i++ ){
        output_file << SDVS[ i ] << ", ";
    }
    output_file << SDVS[ SDVS.size() - 1 ] << "\n";

#endif


    //Begin iteration
    while ( t + dt < tf ){

        std::cout << "t: " << t + dt << "\n";
        time = { t + dt, dt };

        //Increment the displacements
        grad_u_curr   = grad_u_prev   + dt * ( grad_u_f - grad_u_0 );
        phi_curr      = phi_prev      + dt * ( phi_f - phi_0 );
        grad_phi_curr = grad_phi_prev + dt * ( grad_phi_f - grad_phi_0 );

        //Update the arrays
        for ( unsigned int i = 0; i < 3; i++ ){
            for ( unsigned int j = 0; j < 3; j++ ){
                current_grad_u[ i ][ j ]  = grad_u_curr[ i ][ j ];
                previous_grad_u[ i ][ j ] = grad_u_prev[ i ][ j ];
            }
        }

        for ( unsigned int i = 0; i < 9; i++ ){
            current_phi[ i ] = phi_curr[ i ];
            previous_phi[ i ] = phi_prev[ i ];

            for ( unsigned int j = 0; j < 3; j++ ){
                current_grad_phi[ i ][ j ] = grad_phi_curr[ i ][ j ];
                previous_grad_phi[ i ][ j ] = grad_phi_prev[ i ][ j ];
            }
        }

        //Evaluate the model
#ifdef DEBUG_MODE
        DEBUG.clear();
#endif
        int errorCode = material->evaluate_model( time, fparams,
                                                  current_grad_u, current_phi, current_grad_phi,
                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                  SDVS,
                                                  current_ADD_DOF, current_ADD_grad_DOF,
                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                  PK2_result, SIGMA_result, M_result,
                                                  ADD_TERMS,
                                                  output_message
#ifdef DEBUG_MODE
                                                  , DEBUG
#endif
                                                );

        std::cout << "SDVS:\n"; vectorTools::print( SDVS );

#ifdef DEBUG_MODE

        if ( ( fabs( SDVS[ 0 ] ) > 1e-8 ) || ( errorCode != 0 ) ){
            for ( auto inc = DEBUG.begin(); inc != DEBUG.end(); inc++ ){
                std::cout << inc->first << "\n";
                for ( auto itr = inc->second.begin(); itr != inc->second.end(); itr++ ){
                    if ( itr->first.compare( "converged_values" ) != 0 ){
                        std::cout << "    " << itr->first << "\n";
                        std::cout << "        currentMacroGamma:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroGamma" ] );
                        std::cout << "        currentMicroGamma:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGamma" ] );
                        std::cout << "        currentMicroGradientGamma:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGradientGamma" ] );
                        std::cout << "        currentDeformationGradient:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentDeformationGradient" ] );
                        std::cout << "        currentElasticDeformationGradient:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentElasticDeformationGradient" ] );
                        std::cout << "        currentPlasticDeformationGradient:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentPlasticDeformationGradient" ] );
                        std::cout << "        currentPK2Stress:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentPK2Stress" ] );
                        std::cout << "        currentMacroCohesion:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroCohesion" ] );
                        std::cout << "        dMacroCohesiondMacroStrainISV:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "dMacroCohesiondMacroStrainISV" ] );
                        std::cout << "        previousMacroFlowDirection:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "previousMacroFlowDirection" ] );
                        std::cout << "        previousMicroFlowDirection:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "previousMicroFlowDirection" ] );
                        std::cout << "        previousMicroGradientFlowDirection:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "previousMicroGradientFlowDirection" ] );
                        std::cout << "        currentMacroFlowDirection:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroFlowDirection" ] );
                        std::cout << "        currentMicroFlowDirection:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroFlowDirection" ] );
                        std::cout << "        currentMicroGradientFlowDirection:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGradientFlowDirection" ] );
                        std::cout << "        expectedMacroStrainISV\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "expectedMacroStrainISV" ] );
                        std::cout << "        expectedMicroStrainISV\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "expectedMicroStrainISV" ] );
                        std::cout << "        expectedMicroGradientStrainISV\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "expectedMicroGradientStrainISV" ] );
                        std::cout << "        macroYieldFunction:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "macroYieldFunction" ] );
                        std::cout << "        microYieldFunction:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "microYieldFunction" ] );
                        std::cout << "        microGradientYieldFunction:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "microGradientYieldFunction" ] );
                    }
                    else{
                        std::cout << "        convergedPlasticDeformationGradient:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "convergedPlasticDeformationGradient" ] );
                        std::cout << "        convergedPlasticMicroDeformation:\n";
                        std::cout << "        "; vectorTools::print( itr->second[ "convergedPlasticMicroDeformation" ] );
                    }
                }
            }
        }

        output_file << "NEW_INCREMENT\n";

        //Output the time
        output_file << time[ 0 ] << ", " << time[ 1 ] << "\n";

        //Output the current gradient of u
        output_file << current_grad_u[ 0 ][ 0 ] << ", " << current_grad_u[ 0 ][ 1 ] << ", " << current_grad_u[ 0 ][ 2 ] << ", "
                    << current_grad_u[ 1 ][ 0 ] << ", " << current_grad_u[ 1 ][ 1 ] << ", " << current_grad_u[ 1 ][ 2 ] << ", "
                    << current_grad_u[ 2 ][ 0 ] << ", " << current_grad_u[ 2 ][ 1 ] << ", " << current_grad_u[ 2 ][ 2 ] << "\n";

        //Output the current micro displacement
        output_file << current_phi[ 0 ] << ", " << current_phi[ 1 ] << ", " << current_phi[ 2 ] << ", "
                    << current_phi[ 3 ] << ", " << current_phi[ 4 ] << ", " << current_phi[ 5 ] << ", "
                    << current_phi[ 6 ] << ", " << current_phi[ 7 ] << ", " << current_phi[ 8 ] << "\n";

        //Output the current gradient of the micro displacement
        output_file << current_grad_phi[ 0 ][ 0 ] << ", " << current_grad_phi[ 0 ][ 1 ] << ", " << current_grad_phi[ 0 ][ 2 ] << ", "
                    << current_grad_phi[ 1 ][ 0 ] << ", " << current_grad_phi[ 1 ][ 1 ] << ", " << current_grad_phi[ 1 ][ 2 ] << ", "
                    << current_grad_phi[ 2 ][ 0 ] << ", " << current_grad_phi[ 2 ][ 1 ] << ", " << current_grad_phi[ 2 ][ 2 ] << ", "
                    << current_grad_phi[ 3 ][ 0 ] << ", " << current_grad_phi[ 3 ][ 1 ] << ", " << current_grad_phi[ 3 ][ 2 ] << ", "
                    << current_grad_phi[ 4 ][ 0 ] << ", " << current_grad_phi[ 4 ][ 1 ] << ", " << current_grad_phi[ 4 ][ 2 ] << ", "
                    << current_grad_phi[ 5 ][ 0 ] << ", " << current_grad_phi[ 5 ][ 1 ] << ", " << current_grad_phi[ 5 ][ 2 ] << ", "
                    << current_grad_phi[ 6 ][ 0 ] << ", " << current_grad_phi[ 6 ][ 1 ] << ", " << current_grad_phi[ 6 ][ 2 ] << ", "
                    << current_grad_phi[ 7 ][ 0 ] << ", " << current_grad_phi[ 7 ][ 1 ] << ", " << current_grad_phi[ 7 ][ 2 ] << ", "
                    << current_grad_phi[ 8 ][ 0 ] << ", " << current_grad_phi[ 8 ][ 1 ] << ", " << current_grad_phi[ 8 ][ 2 ] << "\n";

        //Output the PK2 stress
        output_file << PK2_result[ 0 ] << ", " << PK2_result[ 1 ] << ", " << PK2_result[ 2 ] << ", "
                    << PK2_result[ 3 ] << ", " << PK2_result[ 4 ] << ", " << PK2_result[ 5 ] << ", "
                    << PK2_result[ 6 ] << ", " << PK2_result[ 7 ] << ", " << PK2_result[ 8 ] << "\n";

        //Output the SIGMA stress
        output_file << SIGMA_result[ 0 ] << ", " << SIGMA_result[ 1 ] << ", " << SIGMA_result[ 2 ] << ", "
                    << SIGMA_result[ 3 ] << ", " << SIGMA_result[ 4 ] << ", " << SIGMA_result[ 5 ] << ", "
                    << SIGMA_result[ 6 ] << ", " << SIGMA_result[ 7 ] << ", " << SIGMA_result[ 8 ] << "\n";

        //Output the M stress
        output_file << M_result[  0 ] << ", " << M_result[  1 ] << ", " << M_result[  2 ] << ", "
                    << M_result[  3 ] << ", " << M_result[  4 ] << ", " << M_result[  5 ] << ", "
                    << M_result[  6 ] << ", " << M_result[  7 ] << ", " << M_result[  8 ] << ", "
                    << M_result[  9 ] << ", " << M_result[ 10 ] << ", " << M_result[ 11 ] << ", "
                    << M_result[ 12 ] << ", " << M_result[ 13 ] << ", " << M_result[ 14 ] << ", "
                    << M_result[ 15 ] << ", " << M_result[ 16 ] << ", " << M_result[ 17 ] << ", "
                    << M_result[ 18 ] << ", " << M_result[ 19 ] << ", " << M_result[ 20 ] << ", "
                    << M_result[ 21 ] << ", " << M_result[ 22 ] << ", " << M_result[ 23 ] << ", "
                    << M_result[ 24 ] << ", " << M_result[ 25 ] << ", " << M_result[ 26 ] << "\n";

        //Determine if there were non-linear iterations
        auto inc = DEBUG.find( "converged_values" );
        if ( inc != DEBUG.end() ){
            PK2_intermediate   = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediatePK2Stress" ];
            SIGMA_intermediate = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceMicroStress" ];
            M_intermediate     = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceHigherOrderStress" ];
        }
        else{
            PK2_intermediate   = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediatePK2Stress" ];
            SIGMA_intermediate = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceMicroStress" ];
            M_intermediate     = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceHigherOrderStress" ];
        }

        //Output the intermediate PK2 stress
        output_file << PK2_intermediate[ 0 ] << ", " << PK2_intermediate[ 1 ] << ", " << PK2_intermediate[ 2 ] << ", "
                    << PK2_intermediate[ 3 ] << ", " << PK2_intermediate[ 4 ] << ", " << PK2_intermediate[ 5 ] << ", "
                    << PK2_intermediate[ 6 ] << ", " << PK2_intermediate[ 7 ] << ", " << PK2_intermediate[ 8 ] << "\n";

        //Output the intermediate SIGMA stress
        output_file << SIGMA_intermediate[ 0 ] << ", " << SIGMA_intermediate[ 1 ] << ", " << SIGMA_intermediate[ 2 ] << ", "
                    << SIGMA_intermediate[ 3 ] << ", " << SIGMA_intermediate[ 4 ] << ", " << SIGMA_intermediate[ 5 ] << ", "
                    << SIGMA_intermediate[ 6 ] << ", " << SIGMA_intermediate[ 7 ] << ", " << SIGMA_intermediate[ 8 ] << "\n";

        //Output the intermediate M stress
        output_file << M_intermediate[  0 ] << ", " << M_intermediate[  1 ] << ", " << M_intermediate[  2 ] << ", "
                    << M_intermediate[  3 ] << ", " << M_intermediate[  4 ] << ", " << M_intermediate[  5 ] << ", "
                    << M_intermediate[  6 ] << ", " << M_intermediate[  7 ] << ", " << M_intermediate[  8 ] << ", "
                    << M_intermediate[  9 ] << ", " << M_intermediate[ 10 ] << ", " << M_intermediate[ 11 ] << ", "
                    << M_intermediate[ 12 ] << ", " << M_intermediate[ 13 ] << ", " << M_intermediate[ 14 ] << ", "
                    << M_intermediate[ 15 ] << ", " << M_intermediate[ 16 ] << ", " << M_intermediate[ 17 ] << ", "
                    << M_intermediate[ 18 ] << ", " << M_intermediate[ 19 ] << ", " << M_intermediate[ 20 ] << ", "
                    << M_intermediate[ 21 ] << ", " << M_intermediate[ 22 ] << ", " << M_intermediate[ 23 ] << ", "
                    << M_intermediate[ 24 ] << ", " << M_intermediate[ 25 ] << ", " << M_intermediate[ 26 ] << "\n";
        
        //Output the state variables
        for ( unsigned int i = 0; i < SDVS.size()-1; i++ ){
            output_file << SDVS[ i ] << ", ";
        }
        output_file << SDVS[ SDVS.size() - 1 ] << "\n";

#endif

        if ( errorCode > 0 ){
            std::cout << output_message << "\n";
            results << "test_evaluate_model_history & False\n";
            output_file.close();
            return 1;
        }

        t += dt;

        grad_u_prev   = grad_u_curr;
        phi_prev      = phi_curr;
        grad_phi_prev = grad_phi_curr;
    
        if ( t > 0.12 ){ output_file.close(); return 1; }
    }
    
#ifdef DEBUG_MODE
    output_file.close();
#endif
    results << "test_evaluate_model_history & True\n";
    return 0;
}

int main(){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or False 
    if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    //Run the tests
    test_computeSecondOrderDruckerPragerYieldEquation( results );
    test_computeHigherOrderDruckerPragerYieldEquation( results );
    test_computeElasticPartOfDeformation( results );
    test_computeElasticDeformationMeasures( results );
    test_computePlasticMacroVelocityGradient( results );
    test_computePlasticMicroVelocityGradient( results );
    test_computePlasticMicroGradientVelocityGradient( results );
    test_computePlasticVelocityGradients( results );
    test_evolvePlasticMicroGradChi( results );
    test_evolvePlasticDeformation( results );
    test_evolveStrainStateVariables( results );
    test_computeFlowDirections( results );
    test_computePlasticDeformationResidual( results );
    test_computePlasticDeformationResidual2( results );
    test_extractMaterialParameters( results );
    test_extractStateVariables( results );
    test_assembleFundamentalDeformationMeasures( results );
    test_evaluateYieldFunctions( results );
    test_computeCohesion( results );
    test_cout_redirect( results );
    test_cerr_redirect( results );

    test_evaluate_model( results );
    test_evaluate_model_continuation( results );

    test_materialLibraryInterface( results );

    test_evaluate_model_history( results );

    //Close the results file
    results.close();

    return 0;
}
