"""
Collection of routines that test the validity of a given parameter set.

Note, this assumes that the linear elasticity models is used and is located in
the same home directory as micromorphic_elasti_plasticity
"""

import os, sys
import numpy as np
import inspect

sys.path.insert( 0, os.path.abspath( '../../../micromorphic_linear_elasticity/src/python/' ) )

import linear_elastic_parameter_constraint_equations as le_pce

def check_elasticity_parameters( parameterDict ):
    """
    Check the elasticity parameters to make sure they are 
    allowable.

    :param dict parameterDict: A dictionary of parameters
    """

    #Create a list of constraint equations
    constraintEquations = [ le_pce.evaluate_g1,  le_pce.evaluate_g2,  le_pce.evaluate_g3,
                            le_pce.evaluate_g4,  le_pce.evaluate_g5,  le_pce.evaluate_g6,
                            le_pce.evaluate_g7,  le_pce.evaluate_g8,  le_pce.evaluate_g9,
                            le_pce.evaluate_g10, le_pce.evaluate_g11, le_pce.evaluate_g12,
                            le_pce.evaluate_g13 ]

    #Create a list of input parameters for the equations
    argSpecs = [ inspect.getfullargspec( ce ) for ce in constraintEquations ]
    
    #Ignore terms that have default arguments. They can't be material parameters.
    parameterNames = [ aS[ 0 ] if ( aS[ 3 ] is None ) else aS[ 0 ][ :( len( aS[ 0 ] ) - len( aS[ 3 ] ) ) ] for aS in argSpecs ]

    #Check that all of the called out parameters are in the parameter dictionary
    if not all ( [ all ( [ ( p in parameterDict.keys() ) for p in pN ] ) for pN in parameterNames ] ):
        outstr = ""
        for i, pN in enumerate( parameterNames ):
            missingValues = [ p for p in pN if not ( p in parameterDict.keys() )]

            if ( len( missingValues) > 0 ):
                outstr += "In equation " + constraintEquations[ i ].__name__ + " missing variables:\n    ";
                for mV in missingValues:
                    outstr += mV + ", ";
                outstr = outstr[: -2 ] + "\n";
        raise ValueError( "Missing variables in parameterDict\n" + outstr )

    #Evaluate the constraint equations
    return dict( [ ( cE.__name__, cE( **dict( [ ( p, parameterDict[ p ] ) for p in pN ] ) )[ 0 ] > 0 ) for pN, cE in zip( parameterNames, constraintEquations) ] )
