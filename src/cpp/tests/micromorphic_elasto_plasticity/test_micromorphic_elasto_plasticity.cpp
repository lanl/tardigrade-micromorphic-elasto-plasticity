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
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 6) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dFdcJ2[j][i] ) ){
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 7) & False\n";
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
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 8) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dFdRCG[j][i] ) ){
                results << "test_computeHigherOrderDruckerPragerYieldEquation (test 9) & False\n";
                return 1;
            }
        }
    }



    results << "test_computeHigherOrderDruckerPragerYieldEquation & True\n";
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

    //Close the results file
    results.close();

    return 0;
}
