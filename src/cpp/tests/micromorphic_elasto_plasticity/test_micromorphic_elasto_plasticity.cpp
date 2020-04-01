//Tests for constitutive_tools

#include<micromorphic_elasto_plasticity.h>
#include<sstream>
#include<fstream>
#include<iostream>

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
        delta[ ( int )( i / 3) ][ i % 3 ] = eps * fabs( currentMicroGradientGamma[ i ] ) + eps;

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

    parameterVector macroFlowParameters = { 0.69418171, 0.48920844, 0.37059555 };
    parameterVector microFlowParameters = { 0.92593374, 0.07788052, 0.70255635 };
    parameterVector microGradientFlowParameters = { 0.75106525, 0.45320165, 0.02355991, 0.68213594, 0.26385375 };

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

    //Close the results file
    results.close();

    return 0;
}
