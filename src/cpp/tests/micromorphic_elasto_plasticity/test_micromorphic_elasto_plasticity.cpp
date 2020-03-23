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
                            std::cout << gradMat[ j ][ 9 * k + 3 * l + m ] << "\n";
                            std::cout << d2FdStress2J2[ j ][ 243 * k + 81 * l + 27 * m + 9 * n + 3 * o + p ] << "\n";
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
                            std::cout << gradMat[ j ][ 9 * k + 3 * l + m ] << "\n";
                            std::cout << d2FdStress2J2[ j ][ 81 * k + 27 * l + 9 * m + 3 * n + o ] << "\n";
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

    variableVector answerGradChie = { 3.43994318,  0.07533594,  1.9900309 ,  2.74890173, -4.07603536,
                                     -1.57164754, -4.68155599,  3.48965518, -8.58540185, -2.61980577,
                                      0.13221771, -3.39243253, -0.96976908,  2.90708208,  1.79870388,
                                     -2.68941213, -4.92137838,  3.02297046, -2.84013306, -4.04337426,
                                     -3.34281866,  2.72965036,  0.752594  ,  2.49143428, -6.444848  ,
                                     -1.27783452, -5.49363904 };

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

    results << "test_computeElasticPartOfDeformation & True\n";
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

    //Close the results file
    results.close();

    return 0;
}
