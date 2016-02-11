# Mappings from program arguments to chars passed into C program as flags
EXACT         = 'e'
BINOM         = 'b'
DENDRIX       = 'd'
PERMUTATIONAL = 'p'
weightFunctionChars = dict(exact=EXACT, binom=BINOM, dendrix=DENDRIX, permutational=PERMUTATIONAL)
weightFunctionNames = dict( (v, k) for k, v in weightFunctionChars.iteritems() )
