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

/*struct cout_redirect{
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
*/

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

int test_computeResidual( std::ofstream &results ){
    /*!
     * Test the computation of the residual equation
     *
     * TODO: Update with better values when they are known.
     *
     * :param std::ofstream &results: The output file.
     */

    solverTools::floatVector gammas = { 8.79028423e-06, 5.52226146e-06, 5.12978173e-06, 5.67171873e-06, 6.46993099e-06 };

    solverTools::floatType Dt = 6.867410424824635;

    solverTools::floatVector currentDeformationGradient = { -0.59712837, -0.74175714,  0.12227677,
                                                            -0.73225042,  0.6025886 , -0.25093911,
                                                             0.14896011, -0.21148619, -0.99899327 };

    solverTools::floatVector currentMicroDeformation = { -0.82792037, -0.32233057,  0.44528056,
                                                          0.37287375,  0.16864785,  0.95050342,
                                                         -0.33873036,  0.90420205, -0.02666455 };

    solverTools::floatVector currentGradientMicroDeformation = { -0.0149665 ,  0.06095584,  0.13338065, -0.43638215,  0.44766035,
                                                                  0.03721981, -0.28652256, -0.2099459 ,  0.21949816,  0.33977289,
                                                                  0.12631343, -0.00438821,  0.17011993,  0.30361245,  0.12610632,
                                                                 -0.37983491, -0.43842411, -0.45938569, -0.43937068,  0.4324248 ,
                                                                 -0.47758355, -0.41184988, -0.37158834, -0.22519797, -0.06346757,
                                                                  0.32424538, -0.09800669 };

    solverTools::floatVector previousPlasticDeformationGradient = { 0.71462444,  0.66839242,  0.26307041,
                                                                    0.29434242,  0.03695804, -0.97809092,
                                                                    0.70705832, -0.77772875,  0.21091525 };

    solverTools::floatVector previousPlasticMicroDeformation = { -0.2443346 ,  0.55472099, -0.7664812 ,
                                                                  0.91741992,  0.36762017,  0.00387774,
                                                                  0.18390445, -0.71529299, -0.66325411 };

    solverTools::floatVector previousPlasticMicroGradient = { -0.00180827, -0.00269321, -0.00372204,  0.00240335, -0.00461926,
                                                               0.00117064,  0.00250045,  0.00061191,  0.00038711, -0.00446384,
                                                               0.00460516, -0.00457497, -0.00492235,  0.00109124,  0.00028529,
                                                               0.00229132,  0.00223341,  0.00179924,  0.00388522,  0.00268035,
                                                              -0.00295786, -0.000225  ,  0.00479787, -0.00461512, -0.00373827,
                                                              -0.00025868, -0.00375377 };

    solverTools::floatVector previousPlasticMacroVelocityGradient = {  2.19422772e-03, -4.84673232e-03,  3.63970379e-03,
                                                                      -1.01913921e-04, -3.44264545e-04, -3.57929500e-03,
                                                                      -4.73359799e-03,  2.05961538e-05,  4.38795320e-03 };

    solverTools::floatVector previousPlasticMicroVelocityGradient = { 0.00158273,  0.0001165 , -0.00275621,
                                                                      0.00230001,  0.00302578, -0.00085188,
                                                                      0.00273759, -0.00303639, -0.00138923 };

    solverTools::floatVector previousPlasticMicroGradientVelocityGradient = { 0.00333778,  0.00288747,  0.00303614,  0.00380015,
                                                                             -0.00176945,  0.00427005,  0.00267057,  0.00243994,
                                                                              0.0014864 , -0.00365863,  0.00095241, -0.00135106,
                                                                              0.00149007, -0.00012389,  0.00280992, -0.00289623,
                                                                             -0.00418976, -0.00169763,  0.0004705 ,  0.00325416,
                                                                             -0.00250118, -0.00423587,  0.00464114,  0.00348601,
                                                                              0.00166716, -0.00284342, -0.00142908 };

    solverTools::floatType previousMacroStrainISV = 0.35364362717230047;

    solverTools::floatType previousMicroStrainISV = 0.7193661140797262;

    solverTools::floatVector previousMicroGradientStrainISV = { 0.86263736, 0.18943229, 0.17855118 };

    solverTools::floatType previousMacroGamma = 8.485955888055519e-07;

    solverTools::floatType previousMicroGamma = 5.234944962734852e-06;

    solverTools::floatVector previousMicroGradientGamma = { 5.30991705e-06, 2.19644682e-06, 2.89037284e-06 };

    solverTools::floatType previousdMacroGdMacroCohesion = 0.8269304350280192;

    solverTools::floatType previousdMicroGdMicroCohesion = 0.3314131216353998;

    solverTools::floatVector previousdMicroGradientGdMicroGradientCohesion = { 0.21254439, 0.        , 0.        ,
                                                                               0.        , 0.21254439, 0.        ,
                                                                               0.        , 0.        , 0.21254439 };

    solverTools::floatVector macroHardeningParameters = { 0.53895133, 0.37172145 };

    solverTools::floatVector microHardeningParameters = { 0.37773052, 0.92739145 };

    solverTools::floatVector microGradientHardeningParameters = { 0.53186824, 0.75454313 };

    solverTools::floatVector macroFlowParameters = { 0.95338442, 0.74042148 };

    solverTools::floatVector microFlowParameters = { 0.38093104, 0.49241325 };

    solverTools::floatVector microGradientFlowParameters = { 0.82121039, 0.90566759 };

    solverTools::floatVector macroYieldParameters = { 0.01166325, 0.05331896 };

    solverTools::floatVector microYieldParameters = { 0.32982199, 0.60161431 };

    solverTools::floatVector microGradientYieldParameters = { 0.58881096, 0.11473813 };

    solverTools::floatVector Amatrix = { 0.5101772 , 0.74331815, 0.54170405, 0.66738607, 0.35637472,
                                         0.64655396, 0.91069292, 0.95300466, 0.82022059, 0.55783418,
                                         0.54157517, 0.83207366, 0.46177772, 0.69092231, 0.76725898,
                                         0.47029883, 0.82405298, 0.3661099 , 0.07293262, 0.51684264,
                                         0.89762081, 0.76841927, 0.39590146, 0.39869885, 0.03044499,
                                         0.70246531, 0.79241141, 0.79776407, 0.40123693, 0.95625319,
                                         0.72248093, 0.08186357, 0.65237553, 0.5281702 , 0.65666359,
                                         0.01927321, 0.06968033, 0.78995691, 0.91746247, 0.10008994,
                                         0.10948784, 0.12832417, 0.78842437, 0.88991891, 0.12318031,
                                         0.21100658, 0.09401477, 0.08875122, 0.43593466, 0.53847729,
                                         0.80731536, 0.44121323, 0.35066743, 0.04239952, 0.86934317,
                                         0.01533633, 0.4818819 , 0.65919776, 0.3688348 , 0.46028593,
                                         0.80390149, 0.40320712, 0.14903849, 0.73070543, 0.19201557,
                                         0.20539311, 0.2651657 , 0.22318962, 0.86158435, 0.05044948,
                                         0.37777148, 0.00819649, 0.60030552, 0.02528907, 0.20399632,
                                         0.91727648, 0.18553634, 0.53538532, 0.38834813, 0.57283409,
                                         0.43195479 };

    solverTools::floatVector Bmatrix = { 0.91158305, 0.790986  , 0.96083758, 0.75514658, 0.92363221,
                                         0.00784569, 0.01689886, 0.8013245 , 0.71237869, 0.90812232,
                                         0.99769622, 0.86699072, 0.71595408, 0.39973814, 0.77767249,
                                         0.50122314, 0.78247379, 0.3746412 , 0.32386691, 0.1931714 ,
                                         0.49819057, 0.59903694, 0.44070301, 0.77242856, 0.77285442,
                                         0.28782501, 0.42986148, 0.45877532, 0.72972672, 0.54542404,
                                         0.2379147 , 0.59603716, 0.37736195, 0.58636566, 0.65543309,
                                         0.01035109, 0.00417787, 0.33885496, 0.57758378, 0.35769777,
                                         0.46650109, 0.99414964, 0.63916312, 0.16170839, 0.0765414 ,
                                         0.88174832, 0.8138165 , 0.00787116, 0.79992681, 0.2065459 ,
                                         0.91834793, 0.88084718, 0.75418913, 0.85539333, 0.05843479,
                                         0.45901681, 0.89307376, 0.31313206, 0.5488607 , 0.24296201,
                                         0.49804496, 0.71486487, 0.75971963, 0.90344445, 0.2070579 ,
                                         0.22442539, 0.82219968, 0.79129232, 0.11145554, 0.20873864,
                                         0.49238459, 0.91907367, 0.854275  , 0.27648138, 0.58150749,
                                         0.64596543, 0.71265391, 0.47846603, 0.43617981, 0.34361584,
                                         0.91137032 };

    solverTools::floatVector Cmatrix = { 7.16078743e-01, 6.54533811e-01, 4.92669701e-01, 7.25673077e-01,
                                         5.09967328e-01, 4.08433385e-01, 8.05655384e-02, 9.80063222e-02,
                                         1.04725354e-01, 7.67136349e-01, 4.97956019e-01, 9.39075819e-01,
                                         6.45106748e-01, 9.21772866e-03, 4.29997409e-01, 4.92448118e-01,
                                         7.02019184e-01, 7.51262806e-02, 9.48261502e-01, 2.57804962e-01,
                                         2.43861442e-01, 5.27705754e-01, 1.93193191e-01, 6.48355536e-01,
                                         8.59922030e-01, 8.84050247e-01, 7.67549820e-01, 7.27032606e-01,
                                         1.72152114e-01, 9.07676842e-01, 2.83243945e-02, 4.30996151e-02,
                                         4.32416941e-01, 8.87628792e-01, 3.17053278e-01, 1.60209715e-01,
                                         9.43791722e-01, 8.49305323e-01, 3.73774568e-01, 8.32175015e-05,
                                         7.38307751e-01, 4.58147026e-01, 6.21058805e-01, 5.11025211e-01,
                                         4.04920559e-01, 6.42759445e-01, 1.72623550e-01, 7.60821653e-02,
                                         8.40239264e-01, 6.93527508e-01, 5.42838621e-01, 9.53003954e-01,
                                         1.05165991e-01, 3.37779673e-01, 1.84495445e-01, 3.70970070e-01,
                                         2.63348719e-01, 8.29303819e-02, 8.28102887e-01, 7.80840856e-01,
                                         9.35823609e-01, 5.20580069e-01, 2.00452456e-01, 6.69073952e-01,
                                         4.87496902e-01, 2.93260668e-01, 8.42782274e-01, 2.94389923e-01,
                                         1.75229486e-01, 7.39695235e-03, 5.89040329e-01, 5.41508718e-01,
                                         9.16831636e-01, 6.29564930e-01, 6.26539192e-01, 5.26590717e-01,
                                         1.04756149e-01, 8.55277101e-01, 3.81779583e-01, 9.85132298e-01,
                                         7.99357251e-01, 8.03072863e-01, 3.15128657e-01, 3.23095171e-01,
                                         9.69129278e-01, 1.27739434e-01, 5.57504096e-01, 9.41152173e-01,
                                         3.25674812e-01, 7.04721543e-01, 2.76820997e-01, 9.87641643e-02,
                                         3.11907575e-01, 9.99604863e-01, 1.26744137e-01, 7.79399580e-01,
                                         8.84852633e-01, 7.92254617e-01, 9.16479277e-01, 8.95196074e-01,
                                         4.43163663e-04, 1.49451416e-01, 7.70233929e-01, 5.94364233e-02,
                                         9.94426969e-01, 1.40389935e-01, 9.90592516e-01, 9.41412729e-01,
                                         5.55187783e-01, 3.11140814e-01, 3.45352885e-01, 2.55041373e-01,
                                         3.82300590e-01, 8.51117530e-01, 6.61397369e-01, 9.41826576e-01,
                                         8.16461015e-01, 8.73566304e-02, 3.23235441e-01, 1.56518549e-02,
                                         4.58382813e-01, 3.17901745e-01, 8.26699741e-01, 5.13889039e-01,
                                         5.30930799e-01, 5.32342370e-01, 8.15407325e-01, 4.03920971e-01,
                                         6.79756593e-01, 2.30255927e-01, 6.26806121e-01, 7.56863745e-01,
                                         6.36493873e-01, 1.31681778e-02, 6.18293024e-01, 9.04757575e-01,
                                         9.57495949e-02, 9.70878995e-01, 9.49321066e-01, 6.91391732e-01,
                                         8.51163894e-01, 2.73203234e-01, 3.48536107e-01, 4.90327008e-01,
                                         3.18865266e-01, 7.33615216e-01, 2.34577072e-01, 1.69258183e-01,
                                         3.47338289e-01, 1.29391192e-01, 5.46278676e-01, 4.32157813e-01,
                                         2.65938801e-01, 7.50066270e-01, 5.10554203e-01, 8.46223050e-01,
                                         2.76502068e-01, 3.55686070e-01, 8.77019523e-01, 3.51840317e-01,
                                         8.72264343e-01, 1.02515259e-01, 5.09074539e-01, 3.78694356e-01,
                                         9.51219452e-01, 1.51467203e-01, 7.85874058e-01, 5.80272582e-01,
                                         8.53099230e-01, 9.92419866e-01, 1.82989325e-01, 6.54718446e-02,
                                         6.16550145e-03, 6.77396427e-01, 3.61848771e-01, 2.31992505e-01,
                                         5.99014216e-02, 5.35632588e-01, 5.58533837e-02, 6.36049056e-01,
                                         2.30503282e-01, 3.03627768e-01, 6.31391927e-01, 1.15726788e-01,
                                         6.34934249e-01, 4.83115809e-01, 5.27998970e-01, 3.28403576e-01,
                                         8.67139032e-01, 9.14908917e-01, 4.81872812e-01, 4.77973647e-01,
                                         9.55639567e-01, 5.53133317e-01, 2.06760015e-01, 1.79061402e-01,
                                         1.50070548e-01, 7.33490866e-01, 7.66185762e-01, 1.13711904e-01,
                                         6.07474359e-01, 4.03380883e-01, 3.53065537e-01, 1.17554817e-01,
                                         2.43909991e-01, 6.76077559e-01, 2.28715506e-01, 8.96349931e-01,
                                         8.63635703e-01, 9.30318072e-01, 9.96620231e-01, 9.42321117e-02,
                                         7.60086060e-01, 4.13949029e-01, 1.69930969e-01, 8.42530581e-02,
                                         3.42661301e-01, 9.91690290e-01, 4.11689438e-01, 7.25256039e-02,
                                         7.94499171e-01, 9.83895530e-01, 6.25047792e-01, 7.94248333e-02,
                                         4.90012260e-01, 5.39185737e-01, 8.52760591e-01, 3.66079202e-01,
                                         5.34651441e-01, 3.80298217e-01, 6.52424789e-02, 3.55121746e-01,
                                         5.93206762e-01, 4.03711507e-01, 4.47993509e-01, 5.41601272e-01,
                                         3.88063795e-01, 3.76234814e-01, 5.10785604e-01, 9.34202309e-01,
                                         7.12643971e-01, 5.41916712e-01, 3.28864944e-01, 9.44482342e-01,
                                         8.64160416e-02, 2.83681250e-01, 1.49416499e-01, 1.59088165e-01,
                                         2.82881897e-02, 6.88674441e-01, 1.34654099e-01, 1.60333272e-01,
                                         5.60124547e-01, 4.78954441e-01, 1.53585867e-01, 2.58217861e-01,
                                         2.78852452e-01, 9.44994445e-01, 5.53084034e-01, 7.92335796e-01,
                                         7.99437385e-01, 5.55374160e-01, 8.13887461e-01, 9.32636355e-01,
                                         3.77691050e-01, 9.08391239e-01, 4.14666777e-01, 6.42549746e-01,
                                         8.79638950e-01, 8.04972696e-01, 1.75427646e-01, 3.65554006e-01,
                                         8.50257955e-01, 1.33044746e-01, 5.76853126e-01, 6.67415597e-01,
                                         1.56866719e-01, 7.00024296e-02, 2.27591053e-02, 5.78150392e-01,
                                         4.68473246e-01, 4.71944456e-01, 4.74897209e-01, 2.99880291e-01,
                                         1.29559241e-01, 9.94036735e-01, 1.96607405e-01, 8.02240205e-01,
                                         6.27781519e-01, 2.44546650e-01, 3.66289594e-01, 9.77416050e-01,
                                         2.39645018e-01, 8.68776402e-01, 8.12732710e-01, 7.64896368e-01,
                                         9.73892955e-01, 2.47565761e-01, 5.42437960e-02, 3.20308787e-01,
                                         8.56517070e-01, 4.19437602e-01, 8.73441598e-01, 5.25224292e-01,
                                         5.22616775e-01, 5.94893952e-01, 1.48649549e-02, 1.62793657e-01,
                                         2.19288138e-01, 8.32614228e-01, 2.45835319e-01, 8.81125304e-01,
                                         1.77871770e-01, 3.01472011e-01, 4.31969236e-01, 6.10635257e-01,
                                         1.33685743e-01, 7.17185802e-01, 3.28887234e-01, 4.54045146e-01,
                                         2.08409820e-01, 9.44744104e-01, 9.10848818e-01, 1.40092949e-01,
                                         7.66015177e-01, 2.55800800e-01, 8.61776934e-01, 8.53936283e-01,
                                         6.14173626e-02, 9.57839377e-01, 1.95874483e-01, 2.61138874e-02,
                                         7.01198729e-01, 9.97427068e-01, 6.03577905e-01, 4.01746154e-01,
                                         9.28278545e-01, 1.36915489e-01, 4.18846935e-01, 8.04181406e-01,
                                         4.45399193e-01, 3.25017398e-01, 7.92906666e-01, 5.61459819e-01,
                                         7.80642957e-01, 4.92125229e-02, 6.42394423e-01, 3.58836920e-01,
                                         6.83028167e-01, 1.28339286e-01, 1.57937966e-01, 1.40092631e-01,
                                         2.85247031e-01, 6.46776464e-01, 6.45188430e-01, 4.11177904e-01,
                                         3.38103070e-01, 9.71812928e-01, 4.42317727e-01, 9.56051389e-01,
                                         7.82289735e-01, 3.33754760e-01, 6.65635978e-01, 3.63705263e-01,
                                         4.17828070e-01, 5.54785953e-01, 4.33981404e-01, 7.99359863e-01,
                                         3.49437115e-01, 9.23036455e-01, 2.19194083e-01, 2.53512378e-01,
                                         6.79483410e-01, 8.79939474e-01, 4.44267000e-02, 2.65684387e-02,
                                         7.74382849e-01, 4.32992472e-01, 7.26891004e-01, 5.70827596e-01,
                                         9.67562572e-01, 7.91047221e-01, 8.49044838e-02, 7.35956104e-01,
                                         2.96818190e-01, 4.69558623e-01, 9.02274540e-01, 7.86645085e-02,
                                         1.52180809e-02, 5.22215007e-01, 9.26919938e-01, 1.30264811e-01,
                                         2.62248809e-01, 5.82398910e-01, 6.86495444e-01, 8.20748499e-01,
                                         6.53023356e-01, 3.16136336e-01, 4.11731766e-01, 5.44675219e-01,
                                         2.55941947e-01, 9.65917791e-01, 7.29965550e-01, 3.81581675e-01,
                                         7.89879700e-01, 7.29741604e-01, 2.15842950e-01, 5.37783289e-01,
                                         4.71351654e-01, 4.37684960e-01, 4.73363495e-01, 8.05173284e-01,
                                         4.10787654e-01, 1.17309902e-01, 2.76157061e-01, 8.91503764e-01,
                                         4.82771924e-01, 4.01176234e-01, 2.06848064e-01, 3.79619446e-01,
                                         4.81267693e-01, 4.08404369e-01, 6.46060775e-02, 6.81302771e-01,
                                         2.06807524e-02, 5.50721334e-01, 1.69866402e-01, 6.62254241e-01,
                                         7.70789136e-01, 6.04696481e-01, 7.17381947e-01, 2.70863752e-02,
                                         2.26605583e-01, 6.30257855e-01, 6.16782503e-01, 5.52706041e-02,
                                         4.59653589e-01, 1.76282762e-01, 2.00411296e-01, 8.00672009e-01,
                                         9.23030810e-01, 7.02394183e-01, 2.98750166e-01, 8.29638680e-01,
                                         5.29039000e-02, 2.18847005e-01, 6.24943327e-02, 7.18858275e-01,
                                         5.85714519e-01, 8.95663866e-01, 9.19433850e-01, 2.23537217e-01,
                                         6.76981065e-01, 7.43375257e-01, 4.73045926e-01, 9.30597866e-01,
                                         7.37347076e-01, 2.98563979e-01, 2.94107804e-01, 2.34378719e-01,
                                         6.66092458e-01, 7.18504051e-02, 2.72176678e-01, 1.23176295e-01,
                                         9.62623312e-01, 2.26542131e-01, 2.25031528e-01, 5.42584180e-01,
                                         1.61734977e-01, 6.27909394e-02, 7.28977994e-01, 8.70259798e-01,
                                         3.59251040e-01, 5.51354983e-01, 9.44685814e-01, 2.97385419e-01,
                                         9.88507574e-01, 3.13744010e-02, 1.91261490e-01, 2.37405622e-01,
                                         1.85690266e-02, 2.00404211e-01, 2.82350295e-01, 6.89846968e-01,
                                         7.44868955e-01, 8.24432237e-01, 2.72873391e-01, 5.08124865e-01,
                                         5.83268108e-01, 1.65433055e-01, 5.22205397e-01, 6.72169054e-01,
                                         4.50928773e-02, 4.54758045e-01, 6.15369929e-01, 1.25110124e-04,
                                         4.30662105e-01, 9.92292082e-01, 6.64003697e-01, 8.54979328e-01,
                                         1.03476619e-01, 6.71162330e-02, 3.48812529e-02, 5.15496293e-01,
                                         5.91575020e-01, 2.83710234e-01, 2.76232369e-01, 4.11397519e-01,
                                         4.28684891e-01, 5.78403169e-01, 6.10435092e-01, 7.68102921e-01,
                                         2.49521172e-01, 8.39248251e-01, 7.92288211e-01, 6.46321925e-01,
                                         5.94022940e-01, 4.71058683e-01, 7.24139826e-01, 8.69411089e-01,
                                         7.83576224e-01, 7.24885753e-01, 5.97115693e-01, 2.56620616e-01,
                                         3.61294906e-03, 3.88906282e-01, 9.45270638e-02, 6.85227130e-01,
                                         1.84505023e-01, 9.34186533e-02, 5.68676781e-01, 4.91224969e-01,
                                         6.97561358e-01, 5.64455150e-01, 6.24639548e-01, 3.89256854e-01,
                                         5.79385338e-01, 9.23713336e-01, 1.78437908e-01, 2.28857497e-01,
                                         8.94885502e-01, 7.19346072e-01, 9.98160778e-01, 4.99114584e-01,
                                         5.11366740e-01, 5.17338403e-01, 1.66754579e-02, 9.14587670e-01,
                                         6.59415885e-01, 5.01295027e-01, 5.32036128e-01, 2.67602006e-01,
                                         2.60363071e-01, 9.07221597e-01, 4.07318593e-01, 7.90591813e-01,
                                         4.92818994e-01, 4.35826271e-01, 2.77098652e-01, 2.31756933e-01,
                                         6.30535318e-01, 2.15473200e-01, 7.42179662e-01, 2.47085990e-01,
                                         5.72589586e-01, 8.82438759e-01, 7.85982404e-01, 3.93700485e-01,
                                         3.04424039e-02, 8.73282135e-01, 3.85698590e-01, 8.80338031e-01,
                                         2.68986810e-01, 8.94467611e-01, 3.77862813e-02, 3.71869345e-01,
                                         1.19102480e-01, 8.33983266e-01, 3.78398773e-01, 1.73555596e-03,
                                         6.24559091e-01, 7.43047324e-01, 8.74409629e-02, 8.73068385e-01,
                                         5.50129965e-02, 5.69450247e-01, 5.36854585e-01, 8.89658769e-01,
                                         6.19142454e-01, 8.04241113e-01, 7.07266413e-01, 3.58608521e-01,
                                         8.20396714e-01, 7.01044055e-01, 3.40214884e-01, 4.59002084e-01,
                                         9.48081040e-01, 8.77344146e-01, 7.06098407e-02, 4.87234971e-01,
                                         7.94733498e-01, 4.59914612e-01, 9.39412345e-02, 5.70547607e-01,
                                         3.85085304e-02, 5.03391297e-01, 4.43798124e-01, 3.20811985e-01,
                                         5.41666556e-01, 9.70856806e-03, 2.34964772e-01, 7.50481327e-01,
                                         7.34870712e-01, 6.02018914e-01, 5.36873081e-01, 5.39976952e-01,
                                         4.64915995e-01, 1.67180838e-02, 8.25868796e-01, 9.48148914e-01,
                                         9.56612362e-01, 4.07612458e-01, 1.28214084e-02, 6.71480905e-01,
                                         4.21558361e-01, 8.51884414e-01, 7.07666712e-01, 1.52372511e-01,
                                         3.03788570e-01, 1.27832625e-02, 4.77840084e-01, 7.71984459e-01,
                                         9.11538510e-01, 9.73664889e-01, 3.71854116e-01, 9.80347165e-01,
                                         6.79372521e-02, 3.00107491e-01, 5.65747566e-01, 8.58660189e-01,
                                         7.30488455e-02, 8.11344503e-01, 4.71219598e-01, 2.89617532e-01,
                                         3.47833011e-01, 9.21814829e-01, 1.25974365e-01, 7.91273030e-01,
                                         4.51719101e-02, 9.33997980e-01, 9.41043477e-01, 2.11900514e-01,
                                         3.21121691e-01, 6.25578075e-01, 8.43845386e-02, 6.04355535e-01,
                                         2.47409502e-01, 4.26046895e-01, 5.42151693e-01, 5.33000823e-01,
                                         1.74906323e-01, 8.57930647e-01, 1.93766797e-01, 6.77369961e-01,
                                         2.89441161e-01, 9.93021531e-01, 3.05527054e-01, 4.57102794e-01,
                                         3.34989433e-01, 5.28267132e-01, 9.93494488e-01, 8.06245534e-01,
                                         9.51331431e-01, 8.94310022e-02, 9.69686311e-02, 4.19514603e-01,
                                         5.03048246e-01, 1.33955968e-01, 2.04824171e-01, 8.24997080e-01,
                                         9.32769271e-01, 6.22541392e-01, 9.86213514e-01, 6.47077549e-01,
                                         2.49569382e-01, 9.53760535e-01, 2.84954064e-01, 4.81259346e-01,
                                         9.25669608e-01, 2.68679507e-01, 2.73417305e-01, 6.95263925e-01,
                                         1.93444813e-01, 3.32621390e-01, 9.43018250e-02, 4.71452005e-02,
                                         1.94488707e-01, 6.28288011e-01, 1.84308812e-01, 5.24017794e-01,
                                         3.80316482e-01, 9.51808512e-01, 2.66610432e-01, 3.05295850e-01,
                                         8.29350794e-01, 4.09690920e-01, 7.24218487e-01, 6.53769477e-01,
                                         2.26116827e-01, 5.26310597e-01, 1.07856020e-01, 2.68010784e-01,
                                         5.81579835e-02, 6.05645745e-01, 2.08305036e-02, 8.42081064e-01,
                                         5.12197761e-01, 6.08640087e-01, 7.08011230e-01, 1.35377840e-01,
                                         4.03393722e-01, 5.67416129e-01, 8.67877545e-01, 6.16676090e-01,
                                         7.03160778e-01, 3.64559987e-01, 7.47423341e-01, 6.23061646e-01,
                                         3.57459133e-01 };

    solverTools::floatVector Dmatrix = { 0.72247484, 0.96837388, 0.16378992, 0.56396088, 0.49388827,
                                         0.29232953, 0.55390765, 0.97434329, 0.79686194, 0.70115803,
                                         0.69670239, 0.31835329, 0.98967071, 0.0710957 , 0.64394145,
                                         0.75447061, 0.62735492, 0.2988142 , 0.39556958, 0.69287339,
                                         0.97507727, 0.63330884, 0.32744789, 0.54114372, 0.63408464,
                                         0.75717759, 0.49859976, 0.93064155, 0.37533856, 0.10242752,
                                         0.2585685 , 0.03115378, 0.15650373, 0.94232915, 0.02525425,
                                         0.51954038, 0.06884208, 0.10213294, 0.2079744 , 0.68143151,
                                         0.625201  , 0.32727099, 0.06991378, 0.6359506 , 0.91870599,
                                         0.39520538, 0.42548738, 0.71426428, 0.39849126, 0.45900502,
                                         0.2899069 , 0.6585042 , 0.36694765, 0.34926767, 0.76545296,
                                         0.06051368, 0.14079775, 0.3398395 , 0.39089776, 0.42926743,
                                         0.15551115, 0.80774478, 0.30957378, 0.59528034, 0.0634203 ,
                                         0.61653922, 0.67869773, 0.27402482, 0.48442433, 0.51676883,
                                         0.44382337, 0.76163941, 0.39436827, 0.64539714, 0.67287057,
                                         0.47043048, 0.62600939, 0.44564747, 0.47883474, 0.2762934 ,
                                         0.44582656 };

    solverTools::floatType alphaMacro = 0.018021292949062184;

    solverTools::floatType alphaMicro = 0.06262948150422443;

    solverTools::floatType alphaMicroGradient = 0.8003116264732453;

    solverTools::floatVector currentElasticDeformationGradient = { 1.00689861,  0.09224087,  0.21654602,
                                                                   0.11113114,  0.75383361, -0.7106512 ,
                                                                   0.13831688, -0.64454366, -0.74094417 };

    solverTools::floatVector currentElasticMicroDeformation = { 0.82618105, -0.01336797,  0.52899063,
                                                                0.42643081, -0.61595002, -0.6270148 ,
                                                               -0.39621571, -0.74087789,  0.53303634 };

    solverTools::floatVector currentElasticMicroGradient = { 2.88193613e-03, -3.52459190e-03, -4.49466929e-03,  2.55546031e-03,
                                                            -1.67398392e-03, -3.39944700e-03, -3.60362996e-03,  3.20305328e-03,
                                                             4.30320373e-04,  2.45914590e-03,  4.98278548e-03, -4.71569558e-03,
                                                             2.55113385e-03,  8.58188313e-04, -2.04621902e-05,  1.53872604e-03,
                                                             3.39648683e-03,  3.75594952e-03,  1.00871697e-03, -1.65978469e-03,
                                                             2.74453064e-03, -1.01917751e-04,  6.86011563e-04, -3.34467284e-03,
                                                            -6.80464731e-04, -2.30576428e-03,  2.91842751e-03 };

    solverTools::floatVector currentPK2Stress = { 0.52380888, 0.21677263, 0.00070288,
                                                  0.23692701, 0.62366402, 0.61310665,
                                                  0.49080997, 0.42909226, 0.5573792 };

    solverTools::floatVector currentReferenceMicroStress = { 0.40006576, 0.75448479, 1.37814846,
                                                             0.75448479, 0.22195633, 0.97292321,
                                                             1.37814846, 0.97292321, 1.26451356 };

    solverTools::floatVector currentReferenceHigherOrderStress = { -0.59756606,  0.81747598, -0.73440323,  0.56693003,  0.68597663,
                                                                   -0.75839129, -0.2002921 ,  0.91988327, -0.54828376,  0.57095722,
                                                                    0.35562692,  0.95820166, -0.11285787,  0.35232342,  0.40325674,
                                                                    0.5408038 ,  0.60630815, -0.53225929, -0.62218873, -0.83222143,
                                                                    0.70922318,  0.22261098,  0.81772968,  0.30289168, -0.19538147,
                                                                    0.02573448, -0.7437062 };

    solverTools::floatType currentMacroStrainISV = 0.5951074462908348; 

    solverTools::floatType currentMicroStrainISV = 0.13847634395281005;

    solverTools::floatVector currentMicroGradientStrainISV = { 0.81024298, 0.06618223, 0.2580753 };

    solverTools::floatVector currentPlasticDeformationGradient = { -0.64232754,  0.76045526,  0.11740517,
                                                                   -0.55963441, -0.37452247, -0.71437547,
                                                                   -0.50872934, -0.50367488,  0.64971862 };

    solverTools::floatVector currentPlasticMicroDeformation = { -0.73611818, -0.34637929,  0.66070056,
                                                                 0.59650727,  0.21262375,  0.75061071,
                                                                -0.37480155,  0.88808075, -0.00981704 };

    solverTools::floatVector currentPlasticMicroGradient = { 9.09221091e-05, -4.97556298e-03,  3.75364875e-03, -6.95765698e-04,
                                                            -2.32785285e-03, -1.73842209e-03, -1.28249451e-03, -3.95809840e-03,
                                                            -4.77385177e-03,  4.00567574e-03,  2.60598145e-03,  4.33336843e-03,
                                                            -4.81096527e-03, -3.75750650e-03, -5.87322930e-04,  1.76290780e-03,
                                                            -3.60791005e-03,  4.95079033e-03, -3.15106866e-04, -1.86616625e-03,
                                                             3.84291870e-03, -2.22612924e-03,  2.47510119e-03,  1.09440750e-03,
                                                             2.43493719e-03, -2.12909924e-03, -3.79114796e-03 };

//    //DEBUGGING
//    std::vector< double > fparams = { 2, 1e2, 1.5e1,               //Macro hardening parameters 
//                                      2, 2e2, 2.0e1,               //Micro hardening parameters
//                                      2, 2.5e2, 2.7e1,             //Micro gradient hardening parameters
//                                      2, 0.56, .79,                //Macro flow parameters
//                                      2, 0.15, -0.40,              //Micro flow parameters
//                                      2, 0.82, 0.90,               //Micro gradient flow parameters
//                                      2, .7, -0.3,                 //Macro yield parameters
//                                      2, .4, .76,                  //Micro yield parameters
//                                      2, 0.52, -.70,                //Micro gradient yield parameters
//                                      2, 696.47, 65.84,            //A stiffness tensor parameters
//                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
//                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
//                                      2, -51.92, 5.13,             //D stiffness tensor parameters
//                                      0.4, 0.3, 0.35, 1e-9, 1e-9   //Integration parameters
//                                    };
//    constantType absoluteTolerance, relativeTolerance;
//    micromorphicElastoPlasticity::extractMaterialParameters( fparams, macroHardeningParameters, microHardeningParameters,
//                                                             microGradientHardeningParameters, macroFlowParameters, microFlowParameters,
//                                                             microGradientFlowParameters, macroYieldParameters, microYieldParameters,
//                                                             microGradientYieldParameters, Amatrix, Bmatrix, Cmatrix, Dmatrix,
//                                                             alphaMacro, alphaMicro, alphaMicroGradient, relativeTolerance, absoluteTolerance );
//
//    Dt = 2.5;
//    currentDeformationGradient = { 1.01, 0, 0, 0, 1, 0, 0, 0, 1 };
//    currentMicroDeformation = { 1.0, 0, 0, 0, 0.98, 0, 0, 0, 1 };
//    currentGradientMicroDeformation = variableVector( 27, 0 );
//    previousPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
//    previousPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
//    previousPlasticMicroGradient = variableVector( 27, 0 );
//    previousPlasticMacroVelocityGradient = variableVector( 9, 0 );
//    previousPlasticMicroVelocityGradient = variableVector( 9, 0 );
//    previousPlasticMicroGradientVelocityGradient = variableVector( 27, 0 );
//    previousMacroStrainISV = 0.;
//    previousMicroStrainISV = 0.;
//    previousMicroGradientStrainISV = variableVector( 3, 0 );
//    previousMacroGamma = 0;
//    previousMicroGamma = 0;
//    previousMicroGradientGamma = variableVector( 3, 0 );
//
//    currentElasticDeformationGradient = currentDeformationGradient;
//    currentElasticMicroDeformation = currentMicroDeformation;
//    currentElasticMicroGradient = currentGradientMicroDeformation;
//    currentPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
//    currentPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
//    currentPlasticMicroGradient = variableVector( 27, 0 );
//
//    currentPK2Stress = { 434.22, 386.082, -531.382, 657.265, -626.479, 438.899, 54.008, 580.606, 255.122 };
//    currentReferenceMicroStress = { 1314.38, 1161.46, -513.317, 1161.46, -1074.57, 1121.81, -513.317, 1121.81, 899.495 };
//    currentReferenceHigherOrderStress  = { -3.13881, 21.1205, -32.2236, 61.3638, 54.4966, 80.7927, -18.3752, 10.5808, 65.8143, 46.5659, 126.087, 19.542, 125.007, 88.992, -34.8056, 10.9744, -93.0911, 15.1465, -44.2434, 71.769, 9.05663, 77.7283, -69.6531, -0.0611322, 52.8567, 32.069, -236.034 };
//
//    gammas = { 0., 0., 0., 0., 0. };
//
//    //DEBUGGING

    solverTools::intVector activePlasticity = { 1, 1, 0, 0, 0 };

    solverTools::floatMatrix floatArgsDefault =
        {
            { Dt },
            currentDeformationGradient,
            currentMicroDeformation,
            currentGradientMicroDeformation,
            previousPlasticDeformationGradient,
            previousPlasticMicroDeformation,
            previousPlasticMicroGradient,
            previousPlasticMacroVelocityGradient,
            previousPlasticMicroVelocityGradient,
            previousPlasticMicroGradientVelocityGradient,
            { previousMacroStrainISV },
            { previousMicroStrainISV },
            previousMicroGradientStrainISV,
            { previousMacroGamma },
            { previousMicroGamma },
            previousMicroGradientGamma,
            { previousdMacroGdMacroCohesion },
            { previousdMicroGdMicroCohesion },
            previousdMicroGradientGdMicroGradientCohesion,
            macroHardeningParameters,
            microHardeningParameters,
            microGradientHardeningParameters,
            macroFlowParameters,
            microFlowParameters,
            microGradientFlowParameters,
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

    solverTools::floatMatrix floatOutsDefault = 
        {
            currentElasticDeformationGradient,
            currentElasticMicroDeformation,
            currentElasticMicroGradient,
            currentPlasticDeformationGradient,
            currentPlasticMicroDeformation,
            currentPlasticMicroGradient,
            currentPK2Stress,
            currentReferenceMicroStress,
            currentReferenceHigherOrderStress,
            { currentMacroStrainISV },
            { currentMicroStrainISV },
            currentMicroGradientStrainISV,
        };

    solverTools::intMatrix intArgs;
    
    solverTools::intMatrix intOutsDefault = { activePlasticity };

    variableVector x = gammas;

    variableVector answerResidual = { 2.77092, 4.0997, 1.31455, 0.730796, 1.20269 };

    variableVector answerElasticMicroGradient = { 0.19740483,  0.17269619, -0.3678423 ,  0.08953459, -0.1405664 ,
                                                 -0.18935258,  0.18035407,  0.24419328,  0.39392356,  0.67019271,
                                                 -0.26558889, -0.08769537,  0.46701901,  0.07696557,  0.09980104,
                                                  0.2731155 , -0.10287654,  0.1019402 , -0.47220379, -0.09295431,
                                                  0.33946221, -0.3351329 ,  0.35476669, -0.6010457 ,  0.30065071,
                                                 -0.03147447, -0.00606308 };

    variableVector resultResidual;
    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intOuts = intOutsDefault;
#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG;
    errorOut error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, resultResidual, floatOuts, intOuts, DEBUG );
#else
    errorOut error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, resultResidual, floatOuts, intOuts );
#endif

    if ( error ){
        error->print();
        results << "test_computeResidual & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals< solverTools::floatType >( resultResidual, answerResidual, 1e-5 ) ){
        results << "test_computeResidual (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( currentDeformationGradient, vectorTools::matrixMultiply( floatOuts[ 0 ], floatOuts[ 3 ], 3, 3, 3, 3 ) ) ){
        results << "test_computeResidual (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( currentMicroDeformation, vectorTools::matrixMultiply( floatOuts[ 1 ], floatOuts[ 4 ], 3, 3, 3, 3 ) ) ){
        results << "test_computeResidual (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerElasticMicroGradient, floatOuts[ 2 ], 1e-5, 1e-4 ) ){
        results << "test_computeResidual (test 4) & False\n";
        return 1;
    }

    //Test the Jacobians
    std::cout << "\nCOMPUTE RESIDUAL WITH JACOBIAN\n";
    solverTools::floatVector resultResidualJ;
    floatArgs = floatArgsDefault;
    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    solverTools::floatMatrix jacobian;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG_JACOBIAN;
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, resultResidualJ,
                                                           jacobian, floatOuts, intOuts, DEBUG_JACOBIAN );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, resultResidualJ,
                                                           jacobian, floatOuts, intOuts );
#endif

    if ( error ){
        error->print();
        results << "test_computeResidual & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals< solverTools::floatType >( resultResidualJ, answerResidual, 1e-5 ) ){
        results << "test_computeResidual (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( currentDeformationGradient, vectorTools::matrixMultiply( floatOuts[ 0 ], floatOuts[ 3 ], 3, 3, 3, 3 ) ) ){
        results << "test_computeResidual (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( currentMicroDeformation, vectorTools::matrixMultiply( floatOuts[ 1 ], floatOuts[ 4 ], 3, 3, 3, 3 ) ) ){
        results << "test_computeResidual (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerElasticMicroGradient, floatOuts[ 2 ], 1e-5, 1e-4 ) ){
        results << "test_computeResidual (test 8) & False\n";
        return 1;
    }

    constantType eps = 1e-6;
#ifdef DEBUG_MODE
    std::cout << "DEBUGGING macro Gamma derivatives\n";
    for ( unsigned int i = 0; i < 1; i++ ){
        solverTools::floatVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector P, M;

        floatArgs = floatArgsDefault;
        floatOuts = floatOutsDefault;
        std::map< std::string, solverTools::floatVector > DEBUG_P, DEBUG_M;
        error = micromorphicElastoPlasticity::computeResidual( x + delta, floatArgs, intArgs, resultResidualJ,
                                                               floatOuts, intOuts, DEBUG_P );

        if ( error ){
            error->print();
            results << "test_computeResidual & False\n";
            return 1;
        }

        floatArgs = floatArgsDefault;
        floatOuts = floatOutsDefault;
        error = micromorphicElastoPlasticity::computeResidual( x - delta, floatArgs, intArgs, resultResidualJ,
                                                               floatOuts, intOuts, DEBUG_M );
    
        if ( error ){
            error->print();
            results << "test_computeResidual & False\n";
            return 1;
        }

        P = DEBUG_P.find( "currentMacroStrainISV" )->second;
        M = DEBUG_M.find( "currentMacroStrainISV" )->second;
        solverTools::floatVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        if ( !vectorTools::fuzzyEquals( gradCol, DEBUG_JACOBIAN[ "dCurrentMacroISVdCurrentMacroGamma" ] ) ){
            std::cerr << "ERROR: dCurrentMacroISVdCurrentMacroGamma is not consistent!\n";
            std::cerr << "result: "; vectorTools::print( gradCol );
            std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dCurrentMacroISVdCurrentMacroGamma" ] );
        }

        P = DEBUG_P.find( "currentMicroStrainISV" )->second;
        M = DEBUG_M.find( "currentMicroStrainISV" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dCurrentMicroISVdCurrentMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: " << 0 << "\n";
            }
        }

        P = DEBUG_P.find( "currentMicroGradientStrainISV" )->second;
        M = DEBUG_M.find( "currentMicroGradientStrainISV" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dCurrentMicroISVdCurrentMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: " << 0 << "\n";
            }
        }

        P = DEBUG_P.find( "currentMacroCohesion_2" )->second;
        M = DEBUG_M.find( "currentMacroCohesion_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMacroCdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMacroCdMacroGamma_2 is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMacroCdMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentMicroCohesion_2" )->second;
        M = DEBUG_M.find( "currentMicroCohesion_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dMicroCdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: " << 0. << "\n";
            }
        }

        P = DEBUG_P.find( "currentMicroGradientCohesion_2" )->second;
        M = DEBUG_M.find( "currentMicroGradientCohesion_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dMicroGradientCdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: " << 0. << "\n";
            }
        }

        P = DEBUG_P.find( "currentPlasticMacroVelocityGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMacroVelocityGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMacroLpdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMacroLpdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMacroLpdMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroVelocityGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMicroVelocityGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dMicroLpdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroGradientVelocityGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMicroGradientVelocityGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dPlasticMicroGradientVelocityGradientdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentPlasticDeformationGradient" )->second;
        M = DEBUG_M.find( "currentPlasticDeformationGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dPlasticFpdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dPlasticFpdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dPlasticFpdMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroDeformation" )->second;
        M = DEBUG_M.find( "currentPlasticMicroDeformation" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dPlasticMicroDeformationdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMicroGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dPlasticMicroGradientdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dPlasticMicroGradientdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dPlasticMicroGradientdMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentElasticDeformationGradient" )->second;
        M = DEBUG_M.find( "currentElasticDeformationGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dElasticFdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dElasticFdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dElasticFdMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentElasticMicroDeformation" )->second;
        M = DEBUG_M.find( "currentElasticMicroDeformation" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dElasticChidMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentElasticMicroGradient" )->second;
        M = DEBUG_M.find( "currentElasticMicroGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dElasticGradChidMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dElasitcGradChidMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dElasticGradChidMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentElasticRightCauchyGreen_2" )->second;
        M = DEBUG_M.find( "currentElasticRightCauchyGreen_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dElasticRCGdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dElasticRCGdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dElasticRCGdMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPK2Stress" )->second;
        M = DEBUG_M.find( "currentPK2Stress" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dPK2dMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dPK2dMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dPK2dMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentReferenceMicroStress" )->second;
        M = DEBUG_M.find( "currentReferenceMicroStress" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dSigmadMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dSigmadMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dSigmadMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentReferenceHigherOrderStress" )->second;
        M = DEBUG_M.find( "currentReferenceHigherOrderStress" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMdMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "macroYieldFunction" )->second;
        M = DEBUG_M.find( "macroYieldFunction" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMacroFdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMacroFdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMacroFdMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "microYieldFunction" )->second;
        M = DEBUG_M.find( "microYieldFunction" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroFdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMicroFdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroFdMacroGamma" ] );
            }
        }

        P = DEBUG_P.find( "microGradientYieldFunction" )->second;
        M = DEBUG_M.find( "microGradientYieldFunction" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroGradientFdMacroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMicroGradientFdMacroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroGradientFdMacroGamma" ] );
            }
        }
    }

    std::cout << "DEBUGGING micro Gamma derivatives\n";
    for ( unsigned int i = 1; i < 2; i++ ){
        solverTools::floatVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector P, M;

        floatArgs = floatArgsDefault;
        floatOuts = floatOutsDefault;
        std::map< std::string, solverTools::floatVector > DEBUG_P, DEBUG_M;
        error = micromorphicElastoPlasticity::computeResidual( x + delta, floatArgs, intArgs, resultResidualJ,
                                                               floatOuts, intOuts, DEBUG_P );

        if ( error ){
            error->print();
            results << "test_computeResidual & False\n";
            return 1;
        }

        floatArgs = floatArgsDefault;
        floatOuts = floatOutsDefault;
        error = micromorphicElastoPlasticity::computeResidual( x - delta, floatArgs, intArgs, resultResidualJ,
                                                               floatOuts, intOuts, DEBUG_M );
    
        if ( error ){
            error->print();
            results << "test_computeResidual & False\n";
            return 1;
        }

        P = DEBUG_P.find( "currentMacroStrainISV" )->second;
        M = DEBUG_M.find( "currentMacroStrainISV" )->second;
        solverTools::floatVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            std::cerr << "ERROR: dCurrentMacroISVdCurrentMicroGamma is not consistent!\n";
            std::cerr << "result: "; vectorTools::print( gradCol );
            std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
        }

        P = DEBUG_P.find( "currentMicroStrainISV" )->second;
        M = DEBUG_M.find( "currentMicroStrainISV" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dCurrentMicroISVdCurrentMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dCurrentMicroISVdCurrentMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: " << 0 << "\n";
            }
        }

        P = DEBUG_P.find( "currentMicroGradientStrainISV" )->second;
        M = DEBUG_M.find( "currentMicroGradientStrainISV" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dCurrentMicroISVdCurrentMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: " << 0 << "\n";
            }
        }

        P = DEBUG_P.find( "currentMacroCohesion_2" )->second;
        M = DEBUG_M.find( "currentMacroCohesion_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dMacroCdMicroGamma_2 is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMacroCdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentMicroCohesion_2" )->second;
        M = DEBUG_M.find( "currentMicroCohesion_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroCdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMicroCdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroCdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentMicroGradientCohesion_2" )->second;
        M = DEBUG_M.find( "currentMicroGradientCohesion_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dMicroGradientCdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentPlasticMacroVelocityGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMacroVelocityGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMacroLpdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMacroLpdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMacroLpdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroVelocityGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMicroVelocityGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroLpdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMicroLpdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroLpdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroGradientVelocityGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMicroGradientVelocityGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroGradientLpdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dPlasticMicroGradientVelocityGradientdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroGradientLpdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPlasticDeformationGradient" )->second;
        M = DEBUG_M.find( "currentPlasticDeformationGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dPlasticFpdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dPlasticFpdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dPlasticFpdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroDeformation" )->second;
        M = DEBUG_M.find( "currentPlasticMicroDeformation" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dPlasticMicroDeformationdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dPlasticMicroDeformationdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dPlasticMicroDeformationdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMicroGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dPlasticMicroGradientdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dPlasticMicroGradientdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dPlasticMicroGradientdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentElasticDeformationGradient" )->second;
        M = DEBUG_M.find( "currentElasticDeformationGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dElasticFdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dElasticFdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dElasticFdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentElasticMicroDeformation" )->second;
        M = DEBUG_M.find( "currentElasticMicroDeformation" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dElasticChidMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dElasticChidMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dElasticChidMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentElasticMicroGradient" )->second;
        M = DEBUG_M.find( "currentElasticMicroGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dElasticGradChidMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dElasticGradChidMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dElasticGradChidMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentElasticRightCauchyGreen_2" )->second;
        M = DEBUG_M.find( "currentElasticRightCauchyGreen_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dElasticRCGdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dElasticRCGdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dElasticRCGdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPK2Stress" )->second;
        M = DEBUG_M.find( "currentPK2Stress" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dPK2dMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dPK2dMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dPK2dMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentReferenceMicroStress" )->second;
        M = DEBUG_M.find( "currentReferenceMicroStress" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dSigmadMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dSigmadMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dSigmadMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentReferenceHigherOrderStress" )->second;
        M = DEBUG_M.find( "currentReferenceHigherOrderStress" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "macroYieldFunction" )->second;
        M = DEBUG_M.find( "macroYieldFunction" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMacroFdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMacroFdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMacroFdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "microYieldFunction" )->second;
        M = DEBUG_M.find( "microYieldFunction" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroFdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMicroFdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroFdMicroGamma" ] );
            }
        }

        P = DEBUG_P.find( "microGradientYieldFunction" )->second;
        M = DEBUG_M.find( "microGradientYieldFunction" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroGradientFdMicroGamma" ][ j ] ) ){
                std::cerr << "ERROR: dMicroGradientFdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroGradientFdMicroGamma" ] );
            }
        }
    }

    std::cout << "DEBUGGING micro gradient gamma derivatives\n";
    for ( unsigned int i = 2; i < 5; i++ ){
        solverTools::floatVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector P, M;

        floatArgs = floatArgsDefault;
        floatOuts = floatOutsDefault;
        std::map< std::string, solverTools::floatVector > DEBUG_P, DEBUG_M;
        error = micromorphicElastoPlasticity::computeResidual( x + delta, floatArgs, intArgs, resultResidualJ,
                                                               floatOuts, intOuts, DEBUG_P );

        if ( error ){
            error->print();
            results << "test_computeResidual & False\n";
            return 1;
        }

        floatArgs = floatArgsDefault;
        floatOuts = floatOutsDefault;
        error = micromorphicElastoPlasticity::computeResidual( x - delta, floatArgs, intArgs, resultResidualJ,
                                                               floatOuts, intOuts, DEBUG_M );
    
        if ( error ){
            error->print();
            results << "test_computeResidual & False\n";
            return 1;
        }

        P = DEBUG_P.find( "currentMacroStrainISV" )->second;
        M = DEBUG_M.find( "currentMacroStrainISV" )->second;
        solverTools::floatVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        if ( !vectorTools::fuzzyEquals( gradCol, variableVector( gradCol.size(), 0 ) ) ){
            std::cerr << "ERROR: dCurrentMacroISVdCurrentMicroGradientGamma is not consistent!\n";
            std::cerr << "result: "; vectorTools::print( gradCol );
            std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
        }

        P = DEBUG_P.find( "currentMicroStrainISV" )->second;
        M = DEBUG_M.find( "currentMicroStrainISV" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dCurrentMicroISVdCurrentMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentMicroGradientStrainISV" )->second;
        M = DEBUG_M.find( "currentMicroGradientStrainISV" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dCurrentMicroGradISVdCurrentMicroGradGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dCurrentMicroISVdCurrentMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer:\n"; vectorTools::print( DEBUG_JACOBIAN[ "dCurrentMicroGradISVdCurrentMicroGradGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentMacroCohesion_2" )->second;
        M = DEBUG_M.find( "currentMacroCohesion_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dMacroCdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentMicroCohesion_2" )->second;
        M = DEBUG_M.find( "currentMicroCohesion_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dMicroGradientCdMicroGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentMicroGradientCohesion_2" )->second;
        M = DEBUG_M.find( "currentMicroGradientCohesion_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroGradientCdMicroGadientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dMicroCdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroGradientCdMicroGadientGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPlasticMacroVelocityGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMacroVelocityGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dMacroLpdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroVelocityGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMicroVelocityGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dMicroLpdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroGradientVelocityGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMicroGradientVelocityGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroGradientLpdMicroGradientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dPlasticMicroGradientVelocityGradientdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroGradientLpdMicroGradientGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentPlasticDeformationGradient" )->second;
        M = DEBUG_M.find( "currentPlasticDeformationGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dPlasticFpdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroDeformation" )->second;
        M = DEBUG_M.find( "currentPlasticMicroDeformation" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dPlasticMicroDeformationdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentPlasticMicroGradient" )->second;
        M = DEBUG_M.find( "currentPlasticMicroGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dPlasticMicroGradientdMicroGradientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dPlasticMicroGradientdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dPlasticMicroGradientdMicroGradientGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentElasticDeformationGradient" )->second;
        M = DEBUG_M.find( "currentElasticDeformationGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dElasticFdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentElasticMicroDeformation" )->second;
        M = DEBUG_M.find( "currentElasticMicroDeformation" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dElasticChidMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentElasticMicroGradient" )->second;
        M = DEBUG_M.find( "currentElasticMicroGradient" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dElasticGradChidMicroGradientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dElasticGradChidMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dElasticGradChidMicroGradientGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentElasticRightCauchyGreen_2" )->second;
        M = DEBUG_M.find( "currentElasticRightCauchyGreen_2" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], 0. ) ){
                std::cerr << "ERROR: dElasticRCGdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( variableVector( gradCol.size(), 0 ) );
            }
        }

        P = DEBUG_P.find( "currentPK2Stress" )->second;
        M = DEBUG_M.find( "currentPK2Stress" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dPK2dMicroGradientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dPK2dMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dPK2dMicroGradientGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentReferenceMicroStress" )->second;
        M = DEBUG_M.find( "currentReferenceMicroStress" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dSigmadMicroGradientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dSigmadMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dSigmadMicroGradientGamma" ] );
            }
        }

        P = DEBUG_P.find( "currentReferenceHigherOrderStress" )->second;
        M = DEBUG_M.find( "currentReferenceHigherOrderStress" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMdMicroGradientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dMdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMdMicroGradientGamma" ] );
            }
        }

        P = DEBUG_P.find( "macroYieldFunction" )->second;
        M = DEBUG_M.find( "macroYieldFunction" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMacroFdMicroGradientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dMacroFdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMacroFdMicroGradientGamma" ] );
            }
        }

        P = DEBUG_P.find( "microYieldFunction" )->second;
        M = DEBUG_M.find( "microYieldFunction" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroFdMicroGradientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dMicroFdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroFdMicroGradientGamma" ] );
            }
        }

        P = DEBUG_P.find( "microGradientYieldFunction" )->second;
        M = DEBUG_M.find( "microGradientYieldFunction" )->second;
        gradCol = ( P - M ) / ( 2 * delta[ i ] );

        if ( gradCol.size() == 0 ){
            std::cout << "Map error!\n";
            assert( 1 == 0 );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], DEBUG_JACOBIAN[ "dMicroGradientFdMicroGradientGamma" ][ 3 * j + i - 2 ] ) ){
                std::cerr << "ERROR: dMicroGradientFdMicroGradientGamma is not consistent!\n";
                std::cerr << "result: "; vectorTools::print( gradCol );
                std::cerr << "answer: "; vectorTools::print( DEBUG_JACOBIAN[ "dMicroGradientFdMicroGradientGamma" ] );
            }
        }
    }
#endif

    for ( unsigned int i = 0; i < x.size(); i++ ){

        constantVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector P, M;

        floatArgs = floatArgsDefault;
        floatOuts = floatOutsDefault;

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector > DEBUG_P, DEBUG_M;
        error = micromorphicElastoPlasticity::computeResidual( x + delta, floatArgs, intArgs, P,
                                                               floatOuts, intOuts, DEBUG_P );
#else
        error = micromorphicElastoPlasticity::computeResidual( x + delta, floatArgs, intArgs, P,
                                                               floatOuts, intOuts );
#endif

        if ( error ){
            error->print();
            results << "test_computeResidual & False\n";
            return 1;
        }

        floatArgs = floatArgsDefault;
        floatOuts = floatOutsDefault;
#ifdef DEBUG_MODE
        error = micromorphicElastoPlasticity::computeResidual( x - delta, floatArgs, intArgs, M,
                                                               floatOuts, intOuts, DEBUG_M );
#else
        error = micromorphicElastoPlasticity::computeResidual( x - delta, floatArgs, intArgs, M,
                                                               floatOuts, intOuts );
#endif
    
        if ( error ){
            error->print();
            results << "test_computeResidual & False\n";
            return 1;
        }

        solverTools::floatVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ] ) ){
                std::cout << "i, j: " << i << ", " << j << "\n";
                std::cout << "numeric:\n"; vectorTools::print( gradCol );
                std::cout << "jacobian:\n"; vectorTools::print( jacobian );
                results << "test_computeResidual (test 9) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeResidual & True\n";
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
    std::vector< double > fparams = { 2, 1e2, 1.5e1,               //Macro hardening parameters 
                                      2, 2e2, 2.0e1,               //Micro hardening parameters
                                      2, 2.5e2, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.,                 //Macro flow parameters
                                      2, 0.15, 0.,                 //Micro flow parameters
                                      2, 0.82, 0.,                 //Micro gradient flow parameters
                                      2, 0.70, 0.,                 //Macro yield parameters
                                      2, 0.40, 0.,                 //Micro yield parameters
                                      2, 0.52, 0.,                 //Micro gradient yield parameters
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

    double current_grad_u[ 3 ][ 3 ] = { {0.200, 0.000, 0.000 },
                                        {0.000, 0.000, 0.000 },
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

    double current_phi[ 9 ] = { 0.600, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000 };

    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    //Initialize the gradient of the micro displacement
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
//    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
//                                           {  0.00294417,  0.34480654, -0.34450988 },
//                                           {  0.21056511, -0.28113967, -0.45726839 },
//                                           { -0.26431286, -0.09985721,  0.47322301 },
//                                           { -0.18156887, -0.32226199, -0.37295847 },
//                                           {  0.15062371,  0.09439471,  0.09167948 },
//                                           { -0.46869859,  0.018301  ,  0.45013866 },
//                                           { -0.15455446,  0.40552715, -0.4216042  },
//                                           { -0.38930237,  0.10974753, -0.31188239 } };

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
    std::vector< double > SDVS( 55, 0 );

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
    std::map< std::string, solverTools::floatVector > DEBUG;
#endif

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
                                                                  , DEBUG
#endif
                                                                  );

    if ( errorCode != 0 ){
        std::cout << output_message;
        results << "test_evaluate_model & False\n";
        return 1;
    }

    std::cout << "SDVS:\n"; vectorTools::print( SDVS );

    results << "test_evaluate_model & True\n";
    return 1;
}

int test_computeStrainISVResidual( std::ofstream &results ){
    /*!
     * Test the computation of the strain ISV residual.
     *
     * :param std::ofstream &results: The output file.
     */

    solverTools::floatVector x = { 1, 1, 1, 1, 1 };

    solverTools::floatMatrix floatArgsDefault =
    {
        { 2.5 },
        { 1.7 },
        { 2.7 },
        { 0.31, 0.162, 0.991 },
        { 9.92430424e-01,  3.70411571e-03,  2.00684100e-04,
          3.70411571e-03,  9.96861203e-01, -1.45757049e-04,
          2.00684100e-04, -1.45757049e-04,  1.00752325e+00 },
        {  1.66749433,  0.35416765, -0.939809  ,  
           1.13282152,  1.67200425,  1.30942005,
           0.91180593, -0.9580783 ,  1.6470536 },
        { -0.95697378,  1.06550409, -0.95387546,
          -1.51083453, -0.45597786,  1.36032498,
          -0.88728256, -1.72034523,  0.53243859 },
        {  0.16953206,  0.16246387,  0.35441079,  0.3743102 , -0.12214526,
          -0.35814699,  0.14110444, -0.04030276,  0.72214869,  0.66960948,
          -0.78846779, -0.87918393,  0.19376438,  0.58478994, -0.54728853,
           0.0704016 , -0.72786767, -0.2555111 , -0.69604652, -0.1403563 ,
           0.58541244, -0.18608698, -0.64430013,  0.81850409,  0.090662  ,
          -0.79900641,  0.43744119 },
        { .28 },
        { .927 },
        { .162, .01782, .772 },
        { .272 },
        { .8761 },
        { .281, .982, 1.271 },
        { .171 },
        { 1.72 },
        { 1.272, 0, 0, 0, 1.272, 0, 0, 0, 1.272 },
        { .272,  0.80 },
        { .272,  0.15 },
        { .272,  0.60 },
        { .272, -0.20 },
        { .272, -0.30 },
        { .272,  0.92 },
        { .262 },
        { .620 },
        { .480 }
    };

    solverTools::floatMatrix floatOutsDefault =
    {             
        { 9.271 },
        { 4.176 },
        { .271, .981, 1.32 },
        variableVector(  9, 0 ),
        variableVector(  9, 0 ),
        variableVector( 27, 0 ),
        variableVector( 81, 0 ),
        variableVector( 81, 0 ),
        variableVector( 81, 0 ),
        variableVector( 81, 0 ),
        variableVector( 27 * 27, 0 ),
        variableVector( 27 * 9, 0 )
    };

    solverTools::floatVector residual( 5, 0 );
    solverTools::floatMatrix jacobian( 5, solverTools::floatVector( 5, 0 ) );

    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intArgs, intOuts;

    solverTools::floatVector answer = { 4.26419, 1.55074, -0.380626, 0.260849, 0.964812 };

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG;
    errorOut error = micromorphicElastoPlasticity::computeStrainISVResidual( x, floatArgs, intArgs, residual,
                                                                             jacobian, floatOuts, intOuts, DEBUG );
#else
    errorOut error = micromorphicElastoPlasticity::computeStrainISVResidual( x, floatArgs, intArgs, residual,
                                                                             jacobian, floatOuts, intOuts );
#endif

    if ( error ){
        error->print();
        results << "test_computeStrainISVResidual & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, residual ) ){
        results << "test_computeStrainISVResidual (test 1) & False\n";
        return 1;
    }

    //Test the jacobians
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < x.size(); i++ ){
        constantVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector P, M;
        solverTools::floatMatrix J;
        
        solverTools::floatMatrix fO = floatOutsDefault;
        solverTools::intMatrix iO;

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector> DEBUG_P, DEBUG_M;
        error = micromorphicElastoPlasticity::computeStrainISVResidual( x + delta, floatArgs, intArgs, P,
                                                                        J, fO, iO, DEBUG_P );
#else
        error = micromorphicElastoPlasticity::computeStrainISVResidual( x + delta, floatArgs, intArgs, P,
                                                                        J, fO, iO );
#endif

        if ( error ){
            error->print();
            results << "test_computeStrainISVResidual & False\n";
            return 1;
        }

        fO = floatOutsDefault;
        iO = solverTools::intMatrix( 0 );

#ifdef DEBUG_MODE
        error = micromorphicElastoPlasticity::computeStrainISVResidual( x - delta, floatArgs, intArgs, M,
                                                                        J, fO, iO, DEBUG_M );
#else
        error = micromorphicElastoPlasticity::computeStrainISVResidual( x - delta, floatArgs, intArgs, M,
                                                                        J, fO, iO );
#endif

        if ( error ){
            error->print();
            results << "test_computeStrainISVResidual & False\n";
            return 1;
        }

        solverTools::floatVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ] ) ){
                results << "test_computeStrainISVResidual (test 2 ) & False\n";
                return 1;
            }
        }

    }

    results << "test_computeStrainISVResidual & True\n";
    return 0;
}

int test_solveForStrainISV( std::ofstream &results ){
    /*!
     * Test the solve for the strain-like ISV
     *
     * :param std::ofstream &results: The output file.
     */

    constantType Dt = 2.5;

    variableType currentMacroGamma = 1.;
    variableType currentMicroGamma = 1.;
    variableVector currentMicroGradientGamma = { 1., 1., 1. };

    variableVector currentElasticRightCauchyGreen = { 9.92430424e-01,  3.70411571e-03,  2.00684100e-04,
                                                      3.70411571e-03,  9.96861203e-01, -1.45757049e-04,
                                                      2.00684100e-04, -1.45757049e-04,  1.00752325e+00 };

    variableVector currentPK2Stress = { 1.66749433,  0.35416765, -0.939809  ,  
                                        1.13282152,  1.67200425,  1.30942005,
                                        0.91180593, -0.9580783 ,  1.6470536 };

    variableVector currentReferenceMicroStress = { -0.95697378,  1.06550409, -0.95387546,
                                                   -1.51083453, -0.45597786,  1.36032498,
                                                   -0.88728256, -1.72034523,  0.53243859 };

    variableVector currentReferenceHigherOrderStress = {  0.16953206,  0.16246387,  0.35441079,  0.3743102 , -0.12214526,
                                                         -0.35814699,  0.14110444, -0.04030276,  0.72214869,  0.66960948,
                                                         -0.78846779, -0.87918393,  0.19376438,  0.58478994, -0.54728853, 
                                                          0.0704016 , -0.72786767, -0.2555111 , -0.69604652, -0.1403563 ,
                                                          0.58541244, -0.18608698, -0.64430013,  0.81850409,  0.090662  ,
                                                         -0.79900641,  0.43744119 };

    variableType   previousMacroGamma = 0.280;
    variableType   previousMicroGamma = 0.927;
    variableVector previousMicroGradientGamma = { .162, .01782, .772 };
    variableType   previousMacroStrainISV = 0.272;
    variableType   previousMicroStrainISV = 0.8761;
    variableVector previousMicroGradientStrainISV = { 0.281, 0.982, 1.271 };
    variableType   previousdMacroGdMacroCohesion = 0.171;
    variableType   previousdMicroGdMicroCohesion = 1.720;
    variableMatrix previousdMicroGradientGdMicroGradientCohesion = { { 1.272, 0.000, 0.000 },
                                                                     { 0.000, 1.272, 0.000 },
                                                                     { 0.000, 0.000, 1.272 } };
    variableType   currentMacroStrainISV = 9.271;
    variableType   currentMicroStrainISV = 4.176;
    variableVector currentMicroGradientStrainISV = { .271, .981, 1.32 };

    variableMatrix dMacroFlowDirectiondPK2Stress, dMacroFlowDirectiondElasticRCG;
    variableMatrix dMicroFlowDirectiondReferenceMicroStress, dMicroFlowDirectiondElasticRCG;
    variableMatrix dMicroGradientFlowDirectiondReferenceHigherOrderStress, dMicroGradientFlowDirectiondElasticRCG;

    bool convergeFlag = false;
    bool fatalErrorFlag = false;

    parameterVector macroHardeningParameters = { 0.272, 0.80 };
    parameterVector microHardeningParameters = { 0.372, 0.15 };
    parameterVector microGradientHardeningParameters = { 0.072, 0.60 };
    
    parameterVector macroFlowParameters = { 0.272, 0.80 };
    parameterVector microFlowParameters = { 0.372, 0.15 };
    parameterVector microGradientFlowParameters = { 0.072, 0.60 };

    parameterType alphaMacro = 0.262;
    parameterType alphaMicro = 0.620;
    parameterType alphaMicroGradient = 0.480;

    variableType currentMacroCohesion, currentMicroCohesion;
    variableVector currentMicroGradientCohesion;

    variableVector currentMacroFlowDirection, currentMicroFlowDirection, currentMicroGradientFlowDirection;

    variableType answerMacroISV = 2.94873;
    variableType answerMicroISV = -0.175844;
    variableVector answerMicroGradientISV = { 2.12108, 3.04216, 2.17998 };

    std::map< std::string, solverTools::floatVector > DEBUG;
    errorOut error = micromorphicElastoPlasticity::solveForStrainISV(
                                        Dt, currentMacroGamma, currentMicroGamma, currentMicroGradientGamma,
                                        currentElasticRightCauchyGreen,
                                        currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress,
                                        previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
                                        previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
                                        previousdMacroGdMacroCohesion, previousdMicroGdMicroCohesion,
                                        previousdMicroGradientGdMicroGradientCohesion,
                                        currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                        currentMacroCohesion, currentMicroCohesion, currentMicroGradientCohesion,
                                        currentMacroFlowDirection, currentMicroFlowDirection, currentMicroGradientFlowDirection,
                                        dMacroFlowDirectiondPK2Stress, dMacroFlowDirectiondElasticRCG,
                                        dMicroFlowDirectiondReferenceMicroStress, dMicroFlowDirectiondElasticRCG,
                                        dMicroGradientFlowDirectiondReferenceHigherOrderStress,
                                        dMicroGradientFlowDirectiondElasticRCG,
                                        convergeFlag, fatalErrorFlag,
                                        macroHardeningParameters, microHardeningParameters,
                                        microGradientHardeningParameters,
                                        macroFlowParameters, microFlowParameters,
                                        microGradientFlowParameters,
                                        alphaMacro, alphaMicro, alphaMicroGradient
#ifdef DEBUG_MODE
                                        , DEBUG
#endif
                                        );

    if ( error ){
        error->print();
        results << "test_solveForStrainISV & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMacroISV, currentMacroStrainISV ) ){
        results << "test_solveForStrainISV (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroISV, currentMicroStrainISV ) ){
        results << "test_solveForStrainISV (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMicroGradientISV, currentMicroGradientStrainISV ) ){
        results << "test_solveForStrainISV (test 3) & False\n";
        return 1;
    }

    const solverTools::floatMatrix strainISVResidualFloatArgs =
        {
            { Dt },
            { currentMacroGamma },
            { currentMicroGamma },
            currentMicroGradientGamma,
            currentElasticRightCauchyGreen,
            currentPK2Stress,
            currentReferenceMicroStress,
            currentReferenceHigherOrderStress,
            { previousMacroGamma },
            { previousMicroGamma },
            previousMicroGradientGamma,
            { previousMacroStrainISV },
            { previousMicroStrainISV },
            previousMicroGradientStrainISV,
            { previousdMacroGdMacroCohesion },
            { previousdMicroGdMicroCohesion },
            vectorTools::appendVectors( previousdMicroGradientGdMicroGradientCohesion ),
            macroHardeningParameters,
            microHardeningParameters,
            microGradientHardeningParameters,
            macroFlowParameters,
            microFlowParameters,
            microGradientFlowParameters,
            { alphaMacro },
            { alphaMicro },
            { alphaMicroGradient }
        };

    solverTools::floatVector vec_dMacroFlowDirectiondPK2Stress( 81, 0 );
    solverTools::floatVector vec_dMacroFlowDirectiondElasticRCG( 81, 0 );
    solverTools::floatVector vec_dMicroFlowDirectiondReferenceMicroStress( 81, 0 );
    solverTools::floatVector vec_dMicroFlowDirectiondElasticRCG( 81, 0 );
    solverTools::floatVector vec_dMicroGradientFlowDirectiondReferenceHigherOrderStress( 81 * 27, 0 );
    solverTools::floatVector vec_dMicroGradientFlowDirectiondElasticRCG( 81 * 9, 0 );

    solverTools::floatMatrix strainISVResidualFloatOuts =
    {
        { currentMacroCohesion },
        { currentMicroCohesion },
        currentMicroGradientCohesion,
        currentMacroFlowDirection,
        currentMicroFlowDirection,
        currentMicroGradientFlowDirection,
        vec_dMacroFlowDirectiondPK2Stress,
        vec_dMacroFlowDirectiondElasticRCG,
        vec_dMicroFlowDirectiondReferenceMicroStress,
        vec_dMicroFlowDirectiondElasticRCG,
        vec_dMicroGradientFlowDirectiondReferenceHigherOrderStress,
        vec_dMicroGradientFlowDirectiondElasticRCG
    };

    solverTools::floatVector x = { currentMacroStrainISV, currentMicroStrainISV,
                                   currentMicroGradientStrainISV[ 0 ],
                                   currentMicroGradientStrainISV[ 1 ],
                                   currentMicroGradientStrainISV[ 2 ] };

    solverTools::floatMatrix floatOuts = strainISVResidualFloatOuts;

    solverTools::floatVector residual( 5, 0 );
    solverTools::floatMatrix jacobian( 5, solverTools::floatVector( 5, 0 ) );

    solverTools::intMatrix intOuts, intArgs;

#ifdef DEBUG_MODE
    DEBUG.clear();
    error = micromorphicElastoPlasticity::computeStrainISVResidual( x, strainISVResidualFloatArgs, intArgs, residual, jacobian,
                                                                    floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeStrainISVResidual( x, strainISVResidualFloatArgs, intArgs, residual, jacobian,
                                                                    floatOuts, intOuts );
#endif

    if ( error ){
        error->print();
        results << "test_solveForStrainISV & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( residual, solverTools::floatVector( x.size(), 0 ) ) ){
        results << "test_solveForStrainISV (test 4) & False\n";
        return 1;
    }


    results << "test_solveForStrainISV & True\n";
    return 0;
}

int test_convergence_of_computeResidual( std::ofstream &results ){
    /*!
     * Test the convergence of computeResidual. This test may be removed at a later date.
     *
     * :param std::ofstream &results: The output file.
     */

    std::vector< double > fparams = { 2, 1e2, 1.5e1,               //Macro hardening parameters
                                      2, 2e2, 2.0e1,               //Micro hardening parameters
                                      2, 2.5e2, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.,                 //Macro flow parameters
                                      2, 0.15, 0.,                 //Micro flow parameters
                                      2, 0.82, 0.,                 //Micro gradient flow parameters
                                      2, 0.70, 0.,                 //Macro yield parameters
                                      2, 0.40, 0.,                 //Micro yield parameters
                                      2, 0.52, 0.,                 //Micro gradient yield parameters
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
        results << "test_convergence_of_computeResidual & False\n";
        return 1;
    }

    constantType Dt = 2.5;
    variableType   currentMacroStrainISV = 0.;
    variableType   currentMicroStrainISV = 0.;
    variableVector currentMicroGradientStrainISV( 3., 0. );
    variableVector currentDeformationGradient = { 1.2, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentMicroDeformation = { 1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0 };
    variableVector currentGradientMicroDeformation = variableVector( 27, 0 );
    variableVector previousPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticMicroGradient = variableVector( 27, 0 );
    variableVector previousPlasticMacroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroGradientVelocityGradient = variableVector( 27, 0 );
    variableType   previousdMacroGdMacroCohesion = 0.;
    variableType   previousdMicroGdMicroCohesion = 0.;
    variableVector previousdMicroGradientGdMicroGradientCohesion( 9., 0. );
    variableType   previousMacroStrainISV = 0.;
    variableType   previousMicroStrainISV = 0.;
    variableVector previousMicroGradientStrainISV = variableVector( 3, 0 );
    variableType   previousMacroGamma = 0;
    variableType   previousMicroGamma = 0;
    variableVector previousMicroGradientGamma = variableVector( 3, 0 );

    variableVector currentElasticDeformationGradient = currentDeformationGradient;
    variableVector currentElasticMicroDeformation = currentMicroDeformation;
    variableVector currentElasticMicroGradient = currentGradientMicroDeformation;
    variableVector currentPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticMicroGradient = variableVector( 27, 0 );

    variableVector currentPK2Stress = { 434.22, 386.082, -531.382, 657.265, -626.479, 438.899, 54.008, 580.606, 255.122 };
    variableVector currentReferenceMicroStress = { 1314.38, 1161.46, -513.317, 1161.46, -1074.57, 1121.81, -513.317, 1121.81, 899.495 };
    variableVector currentReferenceHigherOrderStress  = { -3.13881, 21.1205, -32.2236, 61.3638, 54.4966, 80.7927, -18.3752, 10.5808, 65.8143, 46.5659, 126.087, 19.542, 125.007, 88.992, -34.8056, 10.9744, -93.0911, 15.1465, -44.2434, 71.769, 9.05663, 77.7283, -69.6531, -0.0611322, 52.8567, 32.069, -236.034 };

    variableVector gammas = { 0., 0., 0., 0., 0. };

    solverTools::intVector activePlasticity = { 1, 0, 0, 0, 0 };

    solverTools::floatMatrix floatArgsDefault =
        {
            { Dt },
            currentDeformationGradient,
            currentMicroDeformation,
            currentGradientMicroDeformation,
            previousPlasticDeformationGradient,
            previousPlasticMicroDeformation,
            previousPlasticMicroGradient,
            previousPlasticMacroVelocityGradient,
            previousPlasticMicroVelocityGradient,
            previousPlasticMicroGradientVelocityGradient,
            { previousMacroStrainISV },
            { previousMicroStrainISV },
            previousMicroGradientStrainISV,
            { previousMacroGamma },
            { previousMicroGamma },
            previousMicroGradientGamma,
            { previousdMacroGdMacroCohesion },
            { previousdMicroGdMicroCohesion },
            previousdMicroGradientGdMicroGradientCohesion,
            macroHardeningParameters,
            microHardeningParameters,
            microGradientHardeningParameters,
            macroFlowParameters,
            microFlowParameters,
            microGradientFlowParameters,
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

    solverTools::floatMatrix floatOutsDefault = 
        {
            currentElasticDeformationGradient,
            currentElasticMicroDeformation,
            currentElasticMicroGradient,
            currentPlasticDeformationGradient,
            currentPlasticMicroDeformation,
            currentPlasticMicroGradient,
            currentPK2Stress,
            currentReferenceMicroStress,
            currentReferenceHigherOrderStress,
            { currentMacroStrainISV },
            { currentMicroStrainISV },
            currentMicroGradientStrainISV,
        };

    solverTools::intMatrix intArgs;
    
    solverTools::intMatrix intOutsDefault = { activePlasticity };

    variableVector x = gammas;

    variableVector answerResidual = { 2.77092, 4.0997, 1.31455, 0.730796, 1.20269 };

    variableVector answerElasticMicroGradient = { 0.19740483,  0.17269619, -0.3678423 ,  0.08953459, -0.1405664 ,
                                                 -0.18935258,  0.18035407,  0.24419328,  0.39392356,  0.67019271,
                                                 -0.26558889, -0.08769537,  0.46701901,  0.07696557,  0.09980104,
                                                  0.2731155 , -0.10287654,  0.1019402 , -0.47220379, -0.09295431,
                                                  0.33946221, -0.3351329 ,  0.35476669, -0.6010457 ,  0.30065071,
                                                 -0.03147447, -0.00606308 };

    solverTools::floatVector residual, tol;
    solverTools::floatMatrix jacobian;

    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intOuts = intOutsDefault;
#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG;
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, floatOuts, intOuts );
#endif

    std::cout << "residual: "; vectorTools::print( residual );
    tol = relativeTolerance * residual + absoluteTolerance;
    std::cout << "tolerance: "; vectorTools::print( tol );
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );
    std::cout << "intOuts:\n"; vectorTools::print( intOuts );

    std::cout << "\n\nInitial Iteration\n\n";
    std::cout << "gammas: "; vectorTools::print( gammas );
    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;
#ifdef DEBUG_MODE
    DEBUG.clear();
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts );
#endif
    
    std::cout << "residual: "; vectorTools::print( residual );
    std::cout << "jacobian:\n"; vectorTools::print( jacobian );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );
    std::cout << "intOuts:\n"; vectorTools::print( intOuts );

    unsigned int rank;
    solverTools::floatVector dGamma = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "dGamma: "; vectorTools::print( dGamma );

    gammas += dGamma;
    std::cout << "\n\nFirst Iteration\n\n";
    std::cout << "gamma: "; vectorTools::print( gammas );
    x = gammas;
#ifdef DEBUG_MODE
    DEBUG.clear();
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts );
#endif

    std::cout << "residual: "; vectorTools::print( residual );
    std::cout << "jacobian:\n"; vectorTools::print( jacobian );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );
    std::cout << "intOuts:\n"; vectorTools::print( intOuts );

    dGamma = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "dGamma: "; vectorTools::print( dGamma );

    gammas += dGamma;
    std::cout << "\n\nSecond Iteration\n\n";
    std::cout << "gamma: "; vectorTools::print( gammas );
    x = gammas;

#ifdef DEBUG_MODE
    DEBUG.clear();
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts );
#endif

    std::cout << "residual: "; vectorTools::print( residual );
    std::cout << "jacobian:\n"; vectorTools::print( jacobian );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );
    std::cout << "intOuts:\n"; vectorTools::print( intOuts );

    dGamma = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "dGamma: "; vectorTools::print( dGamma );

    gammas += dGamma;
    std::cout << "\n\nThird Iteration\n\n";
    std::cout << "gamma: "; vectorTools::print( gammas );
    x = gammas;

#ifdef DEBUG_MODE
    DEBUG.clear();
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts );
#endif

    std::cout << "residual: "; vectorTools::print( residual );
    std::cout << "jacobian:\n"; vectorTools::print( jacobian );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );
    std::cout << "intOuts:\n"; vectorTools::print( intOuts );

    dGamma = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "dGamma: "; vectorTools::print( dGamma );

    gammas += dGamma;
    std::cout << "\n\nFourth Iteration\n\n";
    std::cout << "gamma: "; vectorTools::print( gammas );
    x = gammas;

#ifdef DEBUG_MODE
    DEBUG.clear();
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts );
#endif

    std::cout << "residual: "; vectorTools::print( residual );
    std::cout << "jacobian:\n"; vectorTools::print( jacobian );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );
    std::cout << "intOuts:\n"; vectorTools::print( intOuts );

    dGamma = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "dGamma: "; vectorTools::print( dGamma );

    gammas += dGamma;
    std::cout << "\n\nFifth Iteration\n\n";
    std::cout << "gamma: "; vectorTools::print( gammas );
    x = gammas;

#ifdef DEBUG_MODE
    DEBUG.clear();
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts );
#endif

    std::cout << "residual: "; vectorTools::print( residual );
    std::cout << "jacobian:\n"; vectorTools::print( jacobian );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );
    std::cout << "intOuts:\n"; vectorTools::print( intOuts );

    dGamma = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "dGamma: "; vectorTools::print( dGamma );

    gammas += dGamma;
    std::cout << "\n\nSixth Iteration\n\n";
    std::cout << "gamma: "; vectorTools::print( gammas );
    x = gammas;

#ifdef DEBUG_MODE
    DEBUG.clear();
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts );
#endif

    std::cout << "residual: "; vectorTools::print( residual );
    std::cout << "jacobian:\n"; vectorTools::print( jacobian );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );
    std::cout << "intOuts:\n"; vectorTools::print( intOuts );

    dGamma = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "dGamma: "; vectorTools::print( dGamma );

    gammas += dGamma;
    std::cout << "\n\nSeventh Iteration\n\n";
    std::cout << "gamma: "; vectorTools::print( gammas );
    x = gammas;

#ifdef DEBUG_MODE
    DEBUG.clear();
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts, DEBUG );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residual, jacobian, floatOuts, intOuts );
#endif

    std::cout << "residual: "; vectorTools::print( residual );
    std::cout << "jacobian:\n"; vectorTools::print( jacobian );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );
    std::cout << "intOuts:\n"; vectorTools::print( intOuts );

    dGamma = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "dGamma: "; vectorTools::print( dGamma );

    //Check for sub variation
    solverTools::floatMatrix f0 = floatOuts;
    solverTools::intMatrix i0 = intOuts;

    solverTools::floatVector gammasA = gammas + dGamma;
    solverTools::floatVector residualA;
    solverTools::floatMatrix jacobianA;
    std::cout << "\n\nEigth Iteration A\n\n";
    std::cout << "gamma: "; vectorTools::print( gammasA );
    x = gammasA;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG_A;
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residualA, f0, i0, DEBUG_A );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residualA, f0, i0 );
#endif

    std::cout << "residual: "; vectorTools::print( residualA );
    std::cout << "jacobian:\n"; vectorTools::print( jacobian );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( f0 );
    std::cout << "intOuts:\n"; vectorTools::print( i0 );

    f0 = floatOuts;
    i0 = intOuts;

    solverTools::floatVector gammasB = gammas + 0.5 * dGamma;
    solverTools::floatVector residualB;
    solverTools::floatMatrix jacobianB;
    std::cout << "\n\nEigth Iteration B\n\n";
    std::cout << "gamma: "; vectorTools::print( gammasB );
    x = gammasB;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG_B;
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residualB, f0, i0, DEBUG_B );
#else
    error = micromorphicElastoPlasticity::computeResidual( x, floatArgs, intArgs, residualB, f0, i0 );
#endif

    std::cout << "residual: "; vectorTools::print( residualB );
    std::cout << "jacobian:\n"; vectorTools::print( jacobianB );
    std::cout << "det( jacobian ): " << vectorTools::determinant( vectorTools::appendVectors( jacobian), 5, 5 ) << "\n";
    std::cout << "floatOuts:\n"; vectorTools::print( f0 );
    std::cout << "intOuts:\n"; vectorTools::print( i0 );

    std::cout << "\nCHANGE IN VALUES\n";

//    for ( auto it_A = DEBUG_A.begin(); it_A != DEBUG_A.end(); it_A++ ){
//
//        auto it_B = DEBUG_B.find( it_A->first );
//
//        if ( it_B == DEBUG_B.end() ){
//            std::cout << "WRONG\n";
//            assert( 1 == 0 );
//        }
//
//        std::cout << it_A->first << ":\n";
//        std::cout.precision( 17 );
//        std::cout << "    A    : "; vectorTools::print( it_A->second );
//        std::cout << "    B    : "; vectorTools::print( it_B->second );
//        std::cout << "    Delta: "; vectorTools::print( it_A->second - it_B->second );
//        std::cout << "\n";
//        
//    }
//
//    std::cout << "change in residual: "; vectorTools::print( residualB - residualA );
//
//    assert( 1 == 0 );
//
//    dGamma = -vectorTools::solveLinearSystem( jacobian, residual, rank );
//    std::cout << "dGamma: "; vectorTools::print( dGamma );
//
//    f0 = floatOuts;
//    i0 = intOuts;
//
//    variableType dg = 0.001;
//    variableType g0 = gammas[ 0 ] - 0.01;
//    variableType gf = gammas[ 0 ] + 0.01;
//    variableType g = g0;
//
//    //Check the jacobian of Gamma
//    solverTools::floatType eps = 1e-6;
//
//    std::cout << "\n\nBEGIN NUMERIC JACOBIAN\n\n";
//    for ( unsigned int i = 0; i < gammas.size(); i++ ){
//        solverTools::floatVector delta( gammas.size(), 0 );
//        delta[ i ] = eps * fabs( gammas[ i ] ) + eps;
//
//        f0 = floatOuts;
//        i0 = intOuts;
//
//        solverTools::floatVector residual_P, residual_M;
//
//        solverTools::debugMap DEBUG_P, DEBUG_M;
//        error = micromorphicElastoPlasticity::computeResidual( x + delta, floatArgs, intArgs, residual_P, f0, i0, DEBUG_P );
//
//        if ( error ){
//            error->print();
//            assert( 1 == 0 );
//        }
//
////        std::cout << "f0:\n"; vectorTools::print( f0 );
////        std::cout << "i0:\n"; vectorTools::print( i0 );
//
//        f0 = floatOuts;
//        i0 = intOuts;
//
//        error = micromorphicElastoPlasticity::computeResidual( x - delta, floatArgs, intArgs, residual_M, f0, i0, DEBUG_M );
//
//        if ( error ){
//            error->print();
//            assert( 1 == 0 );
//        }
//
//        solverTools::floatVector gradCol = ( residual_P - residual_M ) / ( 2 * delta[ i ] );
//
//        std::cout << "gradCol:\n"; vectorTools::print( gradCol );
//
//        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
//            if ( !vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ] ) ){
//                results << "test_convergence_of_computeResidual & False\n";
//                return 1;
//            }
//        }
//
//        gradCol = ( DEBUG_P[ "currentMacroCohesion_2" ] - DEBUG_M[ "currentMacroCohesion_2" ] ) / ( 2 * delta[ i ] );
//
//        std::cout << "numeric dMCdGamma: "; vectorTools::print( gradCol );
//        std::cout << "analytic : "; vectorTools::print( DEBUG[ "dMacroCdMacroGamma" ] );
//
//        gradCol = ( DEBUG_P[ "currentPK2Stress" ] - DEBUG_M[ "currentPK2Stress" ] ) / ( 2 * delta[ i ] );
//
//        std::cout << "numeric dPK2dGamma: "; vectorTools::print( gradCol );
//        std::cout << "analytic : "; vectorTools::print( DEBUG[ "dMacroCdMacroGamma" ] );
//    }

    return 0;
}

int test_computeStressResidual( std::ofstream &results ){
    /*!
     * Test the computation of the stress Jacobian
     *
     * :param std::ofstream &results: The output file.
     */

    std::vector< double > fparams = { 2, 1e2, 1.5e1,               //Macro hardening parameters
                                      2, 2e2, 2.0e1,               //Micro hardening parameters
                                      2, 2.5e2, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.,                 //Macro flow parameters
                                      2, 0.15, 0.,                 //Micro flow parameters
                                      2, 0.82, 0.,                 //Micro gradient flow parameters
                                      2, 0.70, 0.,                 //Macro yield parameters
                                      2, 0.40, 0.,                 //Micro yield parameters
                                      2, 0.52, 0.,                 //Micro gradient yield parameters
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

    constantType Dt = 2.5;
    variableType   macroGamma = 0.1;//0.271;
    variableType   microGamma = 0.;//0.132;
    variableVector microGradientGamma = {0, 0, 0};//{ 0.082, 0.091, 0.021 };
    variableType   currentMacroStrainISV = 0.;
    variableType   currentMicroStrainISV = 0.;
    variableVector currentMicroGradientStrainISV( 3., 0. );
    variableVector currentDeformationGradient = { 1.2, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentMicroDeformation = { 1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0 };
    variableVector currentGradientMicroDeformation = variableVector( 27, 0 );
    variableVector previousPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector previousPlasticMicroGradient = variableVector( 27, 0 );
    variableVector previousPlasticMacroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroVelocityGradient = variableVector( 9, 0 );
    variableVector previousPlasticMicroGradientVelocityGradient = variableVector( 27, 0 );
    variableType   previousdMacroGdMacroCohesion = 0.;
    variableType   previousdMicroGdMicroCohesion = 0.;
    variableVector previousdMicroGradientGdMicroGradientCohesion( 9., 0. );
    variableType   previousMacroStrainISV = 0.;
    variableType   previousMicroStrainISV = 0.;
    variableVector previousMicroGradientStrainISV = variableVector( 3, 0 );
    variableType   previousMacroGamma = 0;
    variableType   previousMicroGamma = 0;
    variableVector previousMicroGradientGamma = variableVector( 3, 0 );

    variableVector currentElasticDeformationGradient = currentDeformationGradient;
    variableVector currentElasticMicroDeformation = currentMicroDeformation;
    variableVector currentElasticMicroGradient = currentGradientMicroDeformation;
    variableVector currentPlasticDeformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticMicroDeformation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    variableVector currentPlasticMicroGradient = variableVector( 27, 0 );

    //Compute the first value of the stresses
    variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;
    error = micromorphicLinearElasticity::linearElasticityReference(  currentElasticDeformationGradient,
                                                                      currentElasticMicroDeformation,
                                                                      currentElasticMicroGradient,
                                                                      Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                      currentPK2Stress, currentReferenceMicroStress,
                                                                      currentReferenceHigherOrderStress );

    if ( error ){
        error->print();
        assert( 1 == 0 );
    }

    solverTools::intVector activePlasticity = { 1, 0, 0, 0, 0 };

    solverTools::floatMatrix floatArgsDefault =
        {
            { Dt },
            { macroGamma },
            { microGamma },
            microGradientGamma,
            currentDeformationGradient,
            currentMicroDeformation,
            currentGradientMicroDeformation,
            previousPlasticDeformationGradient,
            previousPlasticMicroDeformation,
            previousPlasticMicroGradient,
            previousPlasticMacroVelocityGradient,
            previousPlasticMicroVelocityGradient,
            previousPlasticMicroGradientVelocityGradient,
            { previousMacroStrainISV },
            { previousMicroStrainISV },
            previousMicroGradientStrainISV,
            { previousMacroGamma },
            { previousMicroGamma },
            previousMicroGradientGamma,
            { previousdMacroGdMacroCohesion },
            { previousdMicroGdMicroCohesion },
            previousdMicroGradientGdMicroGradientCohesion,
            macroHardeningParameters,
            microHardeningParameters,
            microGradientHardeningParameters,
            macroFlowParameters,
            microFlowParameters,
            microGradientFlowParameters,
            Amatrix,
            Bmatrix,
            Cmatrix,
            Dmatrix,
            { alphaMacro },
            { alphaMicro },
            { alphaMicroGradient }
        };

    solverTools::floatMatrix floatOutsDefault = 
        {
            currentElasticDeformationGradient,
            currentElasticMicroDeformation,
            currentElasticMicroGradient,
            currentPlasticDeformationGradient,
            currentPlasticMicroDeformation,
            currentPlasticMicroGradient,
            { currentMacroStrainISV },
            { currentMicroStrainISV },
            currentMicroGradientStrainISV,
        };

    solverTools::floatMatrix floatArgs = floatArgsDefault;
    solverTools::floatMatrix floatOuts = floatOutsDefault;

    solverTools::intMatrix intArgs, intOuts;

#ifdef DEBUG_MODE
    std::map< std::string, solverTools::floatVector > DEBUG;
#endif

    solverTools::floatVector x( 45, 0 );
    for ( unsigned int i = 0; i < currentPK2Stress.size(); i++ ){
        x[ i + 0 ] = currentPK2Stress[ i ];
        x[ i + 9 ] = currentReferenceMicroStress[ i ];
    }
    for ( unsigned int i = 0; i < currentReferenceHigherOrderStress.size(); i++ ){
        x[ i + 18 ] = currentReferenceHigherOrderStress[ i ];
    }

    solverTools::floatVector residual;
    solverTools::floatMatrix jacobian;

    error = micromorphicElastoPlasticity::computeStressResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                 floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                 , DEBUG
#endif
                                                               );

    if ( error ){
        error->print();
        results << "test_computeStressResidual & False\n";
        return 1;
    }

    //Test the Jacobians
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < x.size(); i++ ){
        constantVector delta( 45, 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        solverTools::floatVector residual_P, residual_M;
        solverTools::floatMatrix jacobian_P, jacobian_M;

        solverTools::floatMatrix fO = floatOutsDefault;
        solverTools::intMatrix   iO = {};

#ifdef DEBUG_MODE
        std::map< std::string, solverTools::floatVector > DEBUG_P, DEBUG_M;
#endif

        error = micromorphicElastoPlasticity::computeStressResidual( x + delta, floatArgs, intArgs, 
                                                                     residual_P, jacobian_P, fO, iO
#ifdef DEBUG_MODE
                                                                     , DEBUG_P
#endif
                                                                   );
        if ( error ){
            error->print();
            results << "test_computeStressResidual & False\n";
            return 1;
        }

        fO = floatOutsDefault;
        iO = {};
        error = micromorphicElastoPlasticity::computeStressResidual( x - delta, floatArgs, intArgs, 
                                                                     residual_M, jacobian_M, fO, iO 
#ifdef DEBUG_MODE
                                                                     , DEBUG_M
#endif
                                                                   );
        if ( error ){
            error->print();
            results << "test_computeStressResidual & False\n";
            return 1;
        }

//#ifdef DEBUG_MODE
//        //Construct the numeric jacobians of the debug terms
//
//        std::map< std::string, solverTools::floatVector > numericGradients;
//
//        for ( auto it_P = DEBUG_P.begin(); it_P != DEBUG_P.end(); it_P++ ){
//            auto it_M = DEBUG_M.find( it_P->first );
//            if ( it_M == DEBUG_M.end() ){
//                std::cout << "test coding error\n";
//                assert ( 1 == 0 );
//            }
//
//            numericGradients.emplace( it_P->first, ( it_P->second - it_M->second ) / ( 2 * delta[ i ] ) );
//        }
//
//        //Check the Jacobians w.r.t. the PK2 stress
//        if ( i < 9 ){
//            std::cout << "\n\nJACOBIAN W.R.T. PK2 " << i << "\n\n";
//
//            //Flow directions
//            std::cout << "numericGrad (macro dir): "; vectorTools::print( numericGradients[ "currentMacroFlowDirection" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ],
//                                                DEBUG[ "dMacroFlowDirectiondPK2Stress" ][ 9 * j + i ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentMacroFlowDirection" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dMacroFlowDirectiondPK2Stress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of the macro flow direction\n";
//                }
//            }
//
//            std::cout << "numericGrad (micro dir): "; vectorTools::print( numericGradients[ "currentMicroFlowDirection" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ], 0. ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric:  "; vectorTools::print( numericGradients[ "currentMicroFlowDirection" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of the micro flow direction\n";
//                }
//            }
//
//            std::cout << "numericGrad (micro grad dir): "; vectorTools::print( numericGradients[ "currentMicroGradientFlowDirection" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ], 0. ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric:  "; vectorTools::print( numericGradients[ "currentMicroGradientFlowDirection" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 27, 0 ) );
//                    std::cout << "ERROR in jacobian of the micro gradient flow direction\n";
//                }
//            }
//
//            //Plastic velocity gradients
//            std::cout << "numericGrad (Lp): "; vectorTools::print( numericGradients[ "currentPlasticMacroVelocityGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
//                                                DEBUG[ "dPlasticMacroLdPK2Stress" ][ 9 * j + i ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMacroVelocityGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticMacroLdPK2Stress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of plastic macro velocity gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (micro Lp): "; vectorTools::print( numericGradients[ "currentPlasticMicroVelocityGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ], 0. ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroVelocityGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of plastic micro velocity gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (grad micro Lp): "; vectorTools::print( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ], 0. ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 27, 0 ) );
//                    std::cout << "ERROR in jacobian of plastic micro gradient velocity gradient\n";
//                }
//            }
//
//            //Plastic deformation measures
//            std::cout << "numericGrad (Fp): "; vectorTools::print( numericGradients[ "currentPlasticDeformationGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticDeformationGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticDeformationGradient" ][ j ],
//                                                DEBUG[ "dPlasticFdPK2Stress" ][ 9 * j + i ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticDeformationGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticFdPK2Stress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of plastic deformation gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (chi p): "; vectorTools::print( numericGradients[ "currentPlasticMicroDeformation" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroDeformation" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroDeformation" ][ j ], 0. ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroDeformation" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of plastic micro deformation\n";
//                }
//            }
//
//            std::cout << "numericGrad (grad chi p): "; vectorTools::print( numericGradients[ "currentPlasticMicroGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradient" ][ j ],
//                                                DEBUG[ "dPlasticMicroGradientdPK2Stress" ][ 9 * j + i ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticMicroGradientdPK2Stress" ], 27, 9 ) );
//                    std::cout << "ERROR in jacobian of plastic micro gradient\n";
//                }
//            }
//
//            //Elastic deformation measures
//            std::cout << "numericGrad (Fe): "; vectorTools::print ( numericGradients[ "currentElasticDeformationGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
//                                                DEBUG[ "dElasticFdPK2Stress" ][ 9 * j + i ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentElasticDeformationGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dElasticFdPK2Stress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of elastic deformation gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (chi e): "; vectorTools::print( numericGradients[ "currentElasticMicroDeformation" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ], 0. ) ){
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentElasticMicroDeformation" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of elastic micro deformation\n";
//                }
//            }
//
//            std::cout << "numericGrad (grad chi e): "; vectorTools::print( numericGradients[ "currentElasticMicroGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroGradient" ][ j ],
//                                                DEBUG[ "dElasticGradChidPK2Stress" ][ 9 * j + i ], 1e-6, 1e-9 ) ){
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentElasticMicroGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dElasticGradChidPK2Stress" ], 27, 9 ) );
//                    std::cout << "ERROR in jacobian of elastic micro gradient\n";
//                }
//            }
//
//            //Jacobian of the stress measures
//            std::cout << "numericGrad (PK2): "; vectorTools::print( numericGradients[ "newPK2Stress" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "newPK2Stress" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "newPK2Stress" ][ j ],
//                                                DEBUG[ "dNewPK2StressdPK2Stress" ][ 9 * j + i ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "newPK2Stress " ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dNewPK2StressdPK2Stress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of new PK2 stress\n";
//                }
//            }
//
//            std::cout << "numericGrad (Sigma): "; vectorTools::print( numericGradients[ "newReferenceMicroStress" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "newReferenceMicroStress" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "newReferenceMicroStress" ][ j ],
//                                                DEBUG[ "dNewReferenceMicroStressdPK2Stress" ][ 9 * j + i ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print(  numericGradients[ "newReferenceMicroStress" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dNewReferenceMicroStressdPK2Stress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of new Reference micro stress\n";
//                }
//            }
//
//            std::cout << "numericGrad (M): "; vectorTools::print( numericGradients[ "newReferenceHigherOrderStress" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "newReferenceHigherOrderStress" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "newReferenceHigherOrderStress" ][ j ],
//                                                DEBUG[ "dNewReferenceHigherOrderStressdPK2Stress" ][ 9 * j + i ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "newReferenceHigherOrderStress" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dNewReferenceHigherOrderStressdPK2Stress" ], 27, 9 ) );
//                    std::cout << "ERROR in jacobian of new Reference higher order stress\n";
//                }
//            }
//
//        }
//
//        if ( ( i >= 9 ) && ( i < 18 ) ){
//            std::cout << "\n\nJACOBIAN W.R.T. SIGMA " << i - 9 << "\n\n";
//
//            //Flow directions
//            std::cout << "numericGrad (macro dir): "; vectorTools::print( numericGradients[ "currentMacroFlowDirection" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ], 0. ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentMacroFlowDirection" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of the macro flow direction\n";
//                }
//            }
//
//            std::cout << "numericGrad (micro dir): "; vectorTools::print( numericGradients[ "currentMicroFlowDirection" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ],
//                                                 DEBUG[ "dMicroFlowDirectiondReferenceMicroStress" ][ 9 * j + i - 9 ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric:  "; vectorTools::print( numericGradients[ "currentMicroFlowDirection" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dMicroFlowDirectiondReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of the micro flow direction\n";
//                }
//            }
//
//            std::cout << "numericGrad (micro grad dir): "; vectorTools::print( numericGradients[ "currentMicroGradientFlowDirection" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ], 0. ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric:  "; vectorTools::print( numericGradients[ "currentMicroGradientFlowDirection" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 27, 0 ) );
//                    std::cout << "ERROR in jacobian of the micro gradient flow direction\n";
//                }
//            }
//
//            //Plastic velocity gradients
//            std::cout << "numericGrad (Lp): "; vectorTools::print( numericGradients[ "currentPlasticMacroVelocityGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ],
//                                                DEBUG[ "dPlasticMacroLdReferenceMicroStress" ][ 9 * j + i - 9 ], 1e-6, 1e-9 ) ) {
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMacroVelocityGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticMacroLdReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of plastic macro velocity gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (micro Lp): "; vectorTools::print( numericGradients[ "currentPlasticMicroVelocityGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
//                                                DEBUG[ "dPlasticMicroLdReferenceMicroStress" ][ 9 * j + i - 9 ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroVelocityGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticMicroLdReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of plastic micro velocity gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (grad micro Lp): "; vectorTools::print( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ], 0. ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 27, 0 ) );
//                    std::cout << "ERROR in jacobian of plastic micro gradient velocity gradient\n";
//                }
//            }
//
//            //Plastic deformation measures
//            std::cout << "numericGrad (Fp): "; vectorTools::print( numericGradients[ "currentPlasticDeformationGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticDeformationGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticDeformationGradient" ][ j ],
//                                                DEBUG[ "dPlasticFdReferenceMicroStress" ][ 9 * j + i - 9 ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticDeformationGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticFdReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of plastic deformation gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (chi p): "; vectorTools::print( numericGradients[ "currentPlasticMicroDeformation" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroDeformation" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroDeformation" ][ j ],
//                                                DEBUG[ "dPlasticMicroDeformationdReferenceMicroStress" ][ 9 * j + i - 9 ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroDeformation" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticMicroDeformationdReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of plastic micro deformation\n";
//                }
//            }
//
//            std::cout << "numericGrad (grad chi p): "; vectorTools::print( numericGradients[ "currentPlasticMicroGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradient" ][ j ],
//                                                DEBUG[ "dPlasticMicroGradientdReferenceMicroStress" ][ 9 * j + i - 9 ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticMicroGradientdReferenceMicroStress" ], 27, 9 ) );
//                    std::cout << "ERROR in jacobian of plastic micro gradient\n";
//                }
//            }
//
//            //Elastic deformation measures
//            std::cout << "numericGrad (Fe): "; vectorTools::print ( numericGradients[ "currentElasticDeformationGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ],
//                                                DEBUG[ "dElasticFdReferenceMicroStress" ][ 9 * j + i - 9 ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentElasticDeformationGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dElasticFdReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of elastic deformation gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (chi e): "; vectorTools::print( numericGradients[ "currentElasticMicroDeformation" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ],
//                                                DEBUG[ "dElasticChidReferenceMicroStress" ][ 9 * j + i - 9 ], 1e-6, 1e-9 ) ){
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentElasticMicroDeformation" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dElasticChidReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of elastic micro deformation\n";
//                }
//            }
//
//            std::cout << "numericGrad (grad chi e): "; vectorTools::print( numericGradients[ "currentElasticMicroGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroGradient" ][ j ],
//                                                DEBUG[ "dElasticGradChidReferenceMicroStress" ][ 9 * j + i - 9 ], 1e-6, 1e-9 ) ){
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentElasticMicroGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dElasticGradChidReferenceMicroStress" ], 27, 9 ) );
//                    std::cout << "ERROR in jacobian of elastic micro gradient\n";
//                }
//            }
//
//            //Jacobian of the stress measures
//            std::cout << "numericGrad (PK2): "; vectorTools::print( numericGradients[ "newPK2Stress" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "newPK2Stress" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "newPK2Stress" ][ j ],
//                                                DEBUG[ "dNewPK2StressdReferenceMicroStress" ][ 9 * j + i - 9 ] ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "newPK2Stress " ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dNewPK2StressdReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of new PK2 stress\n";
//                }
//            }
//
//            std::cout << "numericGrad (Sigma): "; vectorTools::print( numericGradients[ "newReferenceMicroStress" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "newReferenceMicroStress" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "newReferenceMicroStress" ][ j ],
//                                                DEBUG[ "dNewReferenceMicroStressdReferenceMicroStress" ][ 9 * j + i - 9 ] ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print(  numericGradients[ "newReferenceMicroStress" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dNewReferenceMicroStressdReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of new Reference micro stress\n";
//                }
//            }
//
//            std::cout << "numericGrad (M): "; vectorTools::print( numericGradients[ "newReferenceHigherOrderStress" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "newReferenceHigherOrderStress" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "newReferenceHigherOrderStress" ][ j ],
//                                                DEBUG[ "dNewReferenceHigherOrderStressdReferenceMicroStress" ][ 9 * j + i - 9 ] ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "newReferenceHigherOrderStress" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dNewReferenceHigherOrderStressdReferenceMicroStress" ], 27, 9 ) );
//                    std::cout << "ERROR in jacobian of new Reference higher order stress\n";
//                }
//            }
//        }
//
//        if ( i >= 18 ){
//            std::cout << "\n\nJACOBIAN W.R.T. M " << i - 18 << "\n\n";
//
//            //Flow directions
//            std::cout << "numericGrad (macro dir): "; vectorTools::print( numericGradients[ "currentMacroFlowDirection" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentMacroFlowDirection" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMacroFlowDirection" ][ j ], 0. ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentMacroFlowDirection" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of the macro flow direction\n";
//                }
//            }
//
//            std::cout << "numericGrad (micro dir): "; vectorTools::print( numericGradients[ "currentMicroFlowDirection" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroFlowDirection" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroFlowDirection" ][ j ], 0., 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric:  "; vectorTools::print( numericGradients[ "currentMicroFlowDirection" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of the micro flow direction\n";
//                }
//            }
//
//            std::cout << "numericGrad (micro grad dir): "; vectorTools::print( numericGradients[ "currentMicroGradientFlowDirection" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentMicroGradientFlowDirection" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentMicroGradientFlowDirection" ][ j ],
//                                                DEBUG[ "dMicroGradientFlowDirectiondReferenceHigherOrderStress" ][ 27 * j + i - 18 ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric:  "; vectorTools::print( numericGradients[ "currentMicroGradientFlowDirection" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dMicroGradientFlowDirectiondReferenceHigherOrderStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of the micro gradient flow direction\n";
//                }
//            }
//
//            //Plastic velocity gradients
//            std::cout << "numericGrad (Lp): "; vectorTools::print( numericGradients[ "currentPlasticMacroVelocityGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMacroVelocityGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMacroVelocityGradient" ][ j ], 0., 1e-6, 1e-9 ) ) {
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMacroVelocityGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0. ) );
//                    std::cout << "ERROR in jacobian of plastic macro velocity gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (micro Lp): "; vectorTools::print( numericGradients[ "currentPlasticMicroVelocityGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ], 0., 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroVelocityGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticMicroLdReferenceMicroStress" ], 9, 9 ) );
//                    std::cout << "ERROR in jacobian of plastic micro velocity gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (grad micro Lp): "; vectorTools::print( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradientVelocityGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ][ j ],
//                                                DEBUG[ "dPlasticMicroGradientLdReferenceHigherOrderStress" ][ 27 * j + i - 18 ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroGradientVelocityGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticMicroGradientLdReferenceHigherOrderStress" ], 27, 27 ) );
//                    std::cout << "ERROR in jacobian of plastic micro gradient velocity gradient\n";
//                }
//            }
//
//            //Plastic deformation measures
//            std::cout << "numericGrad (Fp): "; vectorTools::print( numericGradients[ "currentPlasticDeformationGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticDeformationGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticDeformationGradient" ][ j ], 0., 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticDeformationGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of plastic deformation gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (chi p): "; vectorTools::print( numericGradients[ "currentPlasticMicroDeformation" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroDeformation" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroDeformation" ][ j ], 0., 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroDeformation" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of plastic micro deformation\n";
//                }
//            }
//
//            std::cout << "numericGrad (grad chi p): "; vectorTools::print( numericGradients[ "currentPlasticMicroGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroGradient" ][ j ],
//                                                DEBUG[ "dPlasticMicroGradientdReferenceHigherOrderStress" ][ 27 * j + i - 18 ], 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentPlasticMicroGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dPlasticMicroGradientdReferenceHigherOrderStress" ], 27, 27 ) );
//                    std::cout << "ERROR in jacobian of plastic micro gradient\n";
//                }
//            }
//
//            //Elastic deformation measures
//            std::cout << "numericGrad (Fe): "; vectorTools::print ( numericGradients[ "currentElasticDeformationGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticDeformationGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticDeformationGradient" ][ j ], 0., 1e-6, 1e-9 ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentElasticDeformationGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of elastic deformation gradient\n";
//                }
//            }
//
//            std::cout << "numericGrad (chi e): "; vectorTools::print( numericGradients[ "currentElasticMicroDeformation" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroDeformation" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroDeformation" ][ j ], 0., 1e-6, 1e-9 ) ){
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentElasticMicroDeformation" ] );
//                    std::cout << "analytic: "; vectorTools::print( variableVector( 9, 0 ) );
//                    std::cout << "ERROR in jacobian of elastic micro deformation\n";
//                }
//            }
//
//            std::cout << "numericGrad (grad chi e): "; vectorTools::print( numericGradients[ "currentElasticMicroGradient" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "currentElasticMicroGradient" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "currentElasticMicroGradient" ][ j ],
//                                                DEBUG[ "dElasticGradChidReferenceHigherOrderStress" ][ 27 * j + i - 18 ], 1e-6, 1e-9 ) ){
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "currentElasticMicroGradient" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dElasticGradChidReferenceHigherOrderStress" ], 27, 27 ) );
//                    std::cout << "ERROR in jacobian of elastic micro gradient\n";
//                }
//            }
//
//            //Jacobian of the stress measures
//            std::cout << "numericGrad (PK2): "; vectorTools::print( numericGradients[ "newPK2Stress" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "newPK2Stress" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "newPK2Stress" ][ j ],
//                                                DEBUG[ "dNewPK2StressdReferenceHigherOrderStress" ][ 27 * j + i - 18 ] ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "newPK2Stress " ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dNewPK2StressdReferenceHigherOrderStress" ], 27, 9 ) );
//                    std::cout << "ERROR in jacobian of new PK2 stress\n";
//                }
//            }
//
//            std::cout << "numericGrad (Sigma): "; vectorTools::print( numericGradients[ "newReferenceMicroStress" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "newReferenceMicroStress" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "newReferenceMicroStress" ][ j ],
//                                                DEBUG[ "dNewReferenceMicroStressdReferenceHigherOrderStress" ][ 27 * j + i - 18 ] ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print(  numericGradients[ "newReferenceMicroStress" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dNewReferenceMicroStressdReferenceHigherOrderStress" ], 27, 9 ) );
//                    std::cout << "ERROR in jacobian of new Reference micro stress\n";
//                }
//            }
//
//            std::cout << "numericGrad (M): "; vectorTools::print( numericGradients[ "newReferenceHigherOrderStress" ] );
//            for ( unsigned int j = 0; j < numericGradients[ "newReferenceHigherOrderStress" ].size(); j++ ){
//                if ( !vectorTools::fuzzyEquals( numericGradients[ "newReferenceHigherOrderStress" ][ j ],
//                                                DEBUG[ "dNewReferenceHigherOrderStressdReferenceHigherOrderStress" ][ 27 * j + i - 18 ] ) ){
//                    std::cout << "i, j: " << i << ", " << j << "\n";
//                    std::cout << "numeric: "; vectorTools::print( numericGradients[ "newReferenceHigherOrderStress" ] );
//                    std::cout << "analytic: "; vectorTools::print( vectorTools::inflate( DEBUG[ "dNewReferenceHigherOrderStressdReferenceHigherOrderStress" ], 27, 27 ) );
//                    std::cout << "ERROR in jacobian of new Reference higher order stress\n";
//                }
//            }
//        }
//#endif

        //Check the residual jacobian
        solverTools::floatVector gradCol = ( residual_P - residual_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ], 1e-6, 1e-7 ) ){
                std::cout << "i, j: " << i << ", " << j << "\n";
                std::cout << "error: " << gradCol[ j ] - jacobian[ j ][ i ] << "\n";
                std::cout << "gradCol:\n"; vectorTools::print( gradCol );
                std::cout << "jacobian:\n"; vectorTools::print( jacobian );
                results << "test_computeStressResidual (test 1) & False\n";
                return 1;
            }
        }
    }

    //Try solving the equation
    
    std::cout << "\nBEGINING SOLVE\n\n";

    std::ofstream residual_out, jacobian_out;
    residual_out.open("first_residual.txt");
    for ( unsigned int i = 0; i < residual.size(); i++ ){
        residual_out << residual[ i ] << ", ";
    }
    residual_out.close();

    jacobian_out.open("first_jacobian.txt");
    jacobian_out << "[\n";
    for ( unsigned int i = 0; i < jacobian.size(); i++ ){
        jacobian_out << "    [ ";
        for ( unsigned int j = 0; j < jacobian[ i ].size(); j++ ){
            jacobian_out << jacobian[ i ][ j ] << ", ";
        }
        jacobian_out << "]\n";
    }
    jacobian_out << "]\n";
    jacobian_out.close();

    std::cout << "residual: " << vectorTools::l2norm( residual ) << "\n";
    std::cout << "det residual: " << vectorTools::determinant( vectorTools::appendVectors( jacobian ), 45, 45 ) << "\n";
    std::vector< double > flat_jacobian = vectorTools::appendVectors( jacobian );
    Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > A1( flat_jacobian.data(), 45, 45 ); 
    Eigen::JacobiSVD<Eigen::MatrixXd> svd1(A1);
    double cond = svd1.singularValues()(0)
        / svd1.singularValues()(svd1.singularValues().size()-1);
    std::cout << "cond: " << cond << "\n";
    unsigned int rank;
    solverTools::floatVector dx = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );

//    std::cout << "dx: "; vectorTools::print( dx );
    x += dx;

    error = micromorphicElastoPlasticity::computeStressResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                 floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                 , DEBUG
#endif
                                                               );

    if ( error ){
        error->print();
        results << "test_computeStressResidual & False\n";
        return 1;
    }

    residual_out.open("second_residual.txt");
    for ( unsigned int i = 0; i < residual.size(); i++ ){
        residual_out << residual[ i ] << ", ";
    }
    residual_out.close();

    jacobian_out.open("second_jacobian.txt");
    jacobian_out << "[\n";
    for ( unsigned int i = 0; i < jacobian.size(); i++ ){
        jacobian_out << "    [ ";
        for ( unsigned int j = 0; j < jacobian[ i ].size(); j++ ){
            jacobian_out << jacobian[ i ][ j ] << ", ";
        }
        jacobian_out << "]\n";
    }
    jacobian_out << "]\n";
    jacobian_out.close();
    std::cout << "residual: " << vectorTools::l2norm( residual ) << "\n";
    std::cout << "det residual: " << vectorTools::determinant( vectorTools::appendVectors( jacobian ), 45, 45 ) << "\n";
    flat_jacobian = vectorTools::appendVectors( jacobian );
    Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > A2( flat_jacobian.data(), 45, 45 ); 
    Eigen::JacobiSVD<Eigen::MatrixXd> svd2(A2);
    double cond2 = svd2.singularValues()(0)
        / svd2.singularValues()(svd2.singularValues().size()-1);
    std::cout << "cond2: " << cond2 << "\n";

    dx = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );

//    std::cout << "dx: "; vectorTools::print( dx );
    x += dx;

    error = micromorphicElastoPlasticity::computeStressResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                 floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                 , DEBUG
#endif
                                                               );

    if ( error ){
        error->print();
        results << "test_computeStressResidual & False\n";
        return 1;
    }

    std::cout << "residual: " << vectorTools::l2norm( residual ) << "\n";
    std::cout << "det residual: " << vectorTools::determinant( vectorTools::appendVectors( jacobian ), 45, 45 ) << "\n";

    dx = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );

//    std::cout << "dx: "; vectorTools::print( dx );
    x += 0.5 * dx;

    error = micromorphicElastoPlasticity::computeStressResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                 floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                 , DEBUG
#endif
                                                               );

    if ( error ){
        error->print();
        results << "test_computeStressResidual & False\n";
        return 1;
    }

    std::cout << "residual: " << vectorTools::l2norm( residual ) << "\n";
    std::cout << "det residual: " << vectorTools::determinant( vectorTools::appendVectors( jacobian ), 45, 45 ) << "\n";

    dx = -vectorTools::solveLinearSystem( jacobian, residual, rank );
    std::cout << "floatOuts:\n"; vectorTools::print( floatOuts );

//    std::cout << "dx: "; vectorTools::print( dx );
    x += 0.5 * dx;

    error = micromorphicElastoPlasticity::computeStressResidual( x, floatArgs, intArgs, residual, jacobian,
                                                                 floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                 , DEBUG
#endif
                                                               );

    if ( error ){
        error->print();
        results << "test_computeStressResidual & False\n";
        return 1;
    }

    std::cout << "residual: " << vectorTools::l2norm( residual ) << "\n";
    std::cout << "det residual: " << vectorTools::determinant( vectorTools::appendVectors( jacobian ), 45, 45 ) << "\n";

    assert( 1 == 0 );

    results << "test_computeStressResidual & True\n";
    return 0;
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

    solverTools::intMatrix intArgs;
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
    solverTools::floatVector residualAnswer = {   -0.018157  ,   0.0150136 , -0.00834171,  0.0143829  ,  0.00306498,
                                                  -0.0253978 ,  -0.0152216 , -0.0199094 , -0.00694793 ,  0.00409966,
                                                   0.00240259,  -0.0214048 ,  0.0189483 ,  0.0076235  , -0.00285891,
                                                   0.0027911 ,  -0.0208777 ,  0.00599161, -0.000942068, -0.0107383 ,
                                                   0.00681374,  -0.0067848 ,  0.0135221 ,  0.00821126 ,  0.0196686 ,
                                                   0.0136814 ,   0.0159186 , -0.00800705, -0.00138683 ,  0.00228973,
                                                  -0.00966137,  -0.0165582 ,  0.00513825, -0.0161062  ,  0.0269113 ,
                                                   0.0188245 ,   0.008651  ,  0.00110232,  0.0116071  , -0.0142204 ,
                                                   0.0244683 ,   0.00267681,  0.0313831 ,  0.0174318  ,  0.00307478,
                                                  -0.0200436 ,  -0.0570817 , -0.0194401 , -0.0371129  , -0.0547857 ,
                                                 381.101     , 464.145     ,  0.011     ,  0.021      ,  0.031 };

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
        iA = { };

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
                                                1e-9, 1e-9 ) ){
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
                                                1e-9, 1e-9 ) ){
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
                                                1e-9, 1e-9 ) ){
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
                                                1e-9, 1e-9 ) ){
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
                                                1e-9, 1e-9 ) ){
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
                    results << "test_computePlasticDeformationResidual (dPlasticGradientMicroVelocityGradientdPlasticGradientMicroDeformation & False\n";
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
                                                1e-9, 1e-9 ) ){
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
//                return 1;
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
                                                    1e-9, 1e-9 ) ){
                        results << "test_computePlasticDeformationResidual (dPlasticMacroVelocityGradientdMicroGamma) & False\n";
                        return 1;
                    }
                }
    
                for ( unsigned int j = 0; j < numericGradients[ "currentPlasticMicroVelocityGradient" ].size(); j++ ){
                    if ( !vectorTools::fuzzyEquals( numericGradients[ "currentPlasticMicroVelocityGradient" ][ j ],
                                                    DEBUG[ "dPlasticMicroVelocityGradientdMicroGamma" ][ j ],
                                                    1e-9, 1e-9 ) ){
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
                        results << "test_computePlasticDeformationResidual (dExpectedPlasticGradientMicroVelocityGradientdMicroGamma & False\n";
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
                        std::cout << "derp\n";
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
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ], 1e-5, 1e-8 ) ){
                std::cout << "i, j: " << i << ", " << j << "\n";
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
//    test_computeResidual( results );
//    test_convergence_of_computeResidual( results );
    test_extractMaterialParameters( results );
    test_extractStateVariables( results );
    test_assembleFundamentalDeformationMeasures( results );
    test_evaluateYieldFunctions( results );
    test_computeCohesion( results );
    test_cout_redirect( results );
    test_cerr_redirect( results );
    test_computeStrainISVResidual( results );
//    test_computeStressResidual( results );
    test_solveForStrainISV( results );

//    test_evaluate_model( results );

    //Close the results file
    results.close();

    return 0;
}
