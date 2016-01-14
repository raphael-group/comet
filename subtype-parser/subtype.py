#!/usr/bin/python

# Load required modules
import os, sys

def get_parser():
    # Parse arguments
    import argparse
    description = 'Converts subtype and core event input data into'\
                  'appropriate input for subtype-comet'
    parser = argparse.ArgumentParser(description=description)

    # Mutation data
    parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-sub', '--subtype', required=True,
                        help='File with a list of subtype for performing subtype-comet.')
    parser.add_argument('-ce', '--core_events', default=None,
                        help='File with a list of core events for performing subtype-comet.')
    parser.add_argument('-td', '--temp_directory', default='./',
                        help='Directory for storing the temporary output of this conversion.')
    return parser


def run( args ):
    # Parsing arguments
    mutationMatrix = args.mutation_matrix
    tempDirectory = args.temp_directory
    subtypeList = args.subtype

    # Parsing core events list, if available
    if args.core_events:
        with open(args.core_events) as f:
            coreEvents = set( l.rstrip() for l in f )
    else:
        coreEvents = set()

    # Parsing subtype information
    from collections import defaultdict
    subtypeDict = defaultdict( lambda: None, dict() )
    subtypes = set()
    with open( subtypeList ) as f:
        sts = [ l.rstrip().split( '\t' ) for l in f if not l.startswith( '#' ) ]
        for p, s in sts:
            subtypeDict[p] = s
            subtypes.add( s )

    # Parsing mutation matrix
    with open( mutationMatrix ) as f:
        MM = [ l.rstrip().split( '\t' ) for l in f if not l.startswith( '#' ) ]

    # Generating output for new mutation matrix with added subtype information
    if tempDirectory[-1] != '/':
        tempDirectory += '/'
    with open( tempDirectory + 'temp.snv', "w" ) as out:
        output = []
        for p in MM:
            p.extend( list( subtypes.difference( set([subtypeDict[p[0]]]) ) ) )
            output.append( '\t'.join( p ) )
        out.write( '\n'.join( output ) )

    # Generating output for combined subtype and core event list
    with open( tempDirectory + 'temp-subtypes.txt', "w" ) as out:
        out.write( '\n'.join( subtypes | coreEvents ) )

if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]) )
