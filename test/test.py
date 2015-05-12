#!/usr/bin/python

"""Test the current version of CoMEt against previous output
   we've shown to be correct."""

if __name__ == "__main__":
	# Load required modules
	import sys, os, json
	sys.path.append('../')
	import comet as C, run_comet as RC, run_exhaustive as RE

	# Hard-coded results of earlier runs
	trueMCMC = [["CDK4(A),CDKN2A(D),TP53 EGFR,PDGFRA(A),PTEN(D)", 36, 30.44954813369114], ["CDK4(A),CDKN2A(D),TP53 EGFR,PTEN,PTEN(D)", 47, 30.31560877273872], ["CDK4(A),CDKN2A(D),TP53 EGFR,PDGFRA(A),RB1(D)", 6, 28.487778643951266], ["CDK4(A),CDKN2A(D),TP53 PDGFRA(A),PTEN,PTEN(D)", 2, 27.986374758844306], ["CDK4(A),CDKN2A(D),TP53 EGFR,PTEN(D),RB1(D)", 1, 27.207736477342962], ["CDK4(A),CDKN2A(D),TP53 EGFR,MDM4(A),RB1(D)", 1, 26.9059191024149], ["CDK4(A),CDKN2A(D),RB1(D) EGFR,PDGFRA(A),PTEN(D)", 3, 26.35048099140887], ["CDK4(A),CDKN2A(D),RB1(D) EGFR,PTEN,PTEN(D)", 1, 26.216541630456447], ["CDK4(A),CDKN2A(D),TP53 MDM4(A),PTEN(D),RB1(D)", 1, 24.99799710738864], ["CDK4(A),CDKN2A(D),TP53 EGFR,MDM4(A),PTEN", 1, 24.910765737533954], ["CDK4(A),CDKN2A(D),TP53 PDGFRA(A),PTEN,RB1(D)", 1, 24.633365961904023]]
	trueExhaustive = [[["CDK4(A)", "CDKN2A(D)", "TP53"], 0.0, 24.61051], [["CDK4(A)", "CDKN2A(D)", "RB1(D)"], 0.0, 20.51144], [["CDKN2A(D)", "RB1(D)", "TP53"], 0.0, 15.30173], [["CDKN2A(D)", "PDGFRA(A)", "TP53"], 0.0, 14.35459], [["CDK4(A)", "CDKN2A(D)", "PTEN(D)"], 1e-05, 11.92406], [["CDKN2A(D)", "PTEN(D)", "TP53"], 1e-05, 11.90466], [["CDK4(A)", "CDKN2A(D)", "EGFR"], 1e-05, 11.49421], [["CDKN2A(D)", "MDM4(A)", "TP53"], 2e-05, 10.68068], [["CDK4(A)", "CDKN2A(D)", "MDM4(A)"], 4e-05, 10.07425], [["CDKN2A(D)", "PTEN", "TP53"], 7e-05, 9.55559], [["CDKN2A(D)", "EGFR", "TP53"], 0.00015, 8.82747], [["CDK4(A)", "CDKN2A(D)", "PTEN"], 0.00035, 7.96747], [["CDK4(A)", "CDKN2A(D)", "PDGFRA(A)"], 0.00048, 7.6389], [["EGFR", "PDGFRA(A)", "PTEN(D)"], 0.00291, 5.83904], [["EGFR", "PTEN", "PTEN(D)"], 0.00333, 5.7051], [["EGFR", "PDGFRA(A)", "RB1(D)"], 0.02071, 3.87727], [["CDK4(A)", "PTEN", "PTEN(D)"], 0.02471, 3.70061], [["CDK4(A)", "EGFR", "RB1(D)"], 0.03291, 3.41404], [["PDGFRA(A)", "PTEN", "PTEN(D)"], 0.03419, 3.37587], [["CDK4(A)", "EGFR", "PTEN(D)"], 0.04911, 3.01377], [["PTEN", "PTEN(D)", "TP53"], 0.0636, 2.75514], [["EGFR", "PTEN(D)", "RB1(D)"], 0.07448, 2.59723], [["EGFR", "PTEN(D)", "TP53"], 0.07453, 2.59651], [["CDKN2A(D)", "EGFR", "RB1(D)"], 0.08983, 2.40983], [["MDM4(A)", "PTEN(D)", "TP53"], 0.09648, 2.33843], [["EGFR", "MDM4(A)", "RB1(D)"], 0.10072, 2.29541], [["CDKN2A(D)", "EGFR", "PTEN"], 0.12616, 2.0702], [["PTEN", "PTEN(D)", "RB1(D)"], 0.13519, 2.00109], [["EGFR", "MDM4(A)", "PDGFRA(A)"], 0.15149, 1.88723], [["CDK4(A)", "EGFR", "PDGFRA(A)"], 0.17147, 1.76336], [["EGFR", "MDM4(A)", "TP53"], 0.18851, 1.66862], [["MDM4(A)", "PTEN", "PTEN(D)"], 0.19338, 1.64308], [["CDKN2A(D)", "MDM4(A)", "RB1(D)"], 0.20638, 1.57806], [["CDK4(A)", "EGFR", "PTEN"], 0.20825, 1.56901], [["PDGFRA(A)", "PTEN(D)", "RB1(D)"], 0.21409, 1.54137], [["CDK4(A)", "PTEN(D)", "RB1(D)"], 0.23869, 1.43258], [["PDGFRA(A)", "PTEN(D)", "TP53"], 0.24857, 1.39204], [["CDKN2A(D)", "PDGFRA(A)", "PTEN"], 0.27821, 1.27939], [["MDM4(A)", "PDGFRA(A)", "TP53"], 0.29567, 1.21851], [["PTEN(D)", "RB1(D)", "TP53"], 0.3062, 1.18352], [["CDKN2A(D)", "PDGFRA(A)", "RB1(D)"], 0.31956, 1.14081], [["CDK4(A)", "EGFR", "MDM4(A)"], 0.3235, 1.12854], [["CDKN2A(D)", "PTEN", "RB1(D)"], 0.41747, 0.87355], [["CDK4(A)", "EGFR", "TP53"], 0.44446, 0.81089], [["MDM4(A)", "PDGFRA(A)", "RB1(D)"], 0.4468, 0.80566], [["CDK4(A)", "MDM4(A)", "RB1(D)"], 0.47088, 0.75315], [["MDM4(A)", "RB1(D)", "TP53"], 0.47495, 0.74454], [["EGFR", "PDGFRA(A)", "TP53"], 0.53346, 0.62837], [["EGFR", "PDGFRA(A)", "PTEN"], 0.53346, 0.62837], [["MDM4(A)", "PDGFRA(A)", "PTEN(D)"], 0.54257, 0.61144], [["EGFR", "MDM4(A)", "PTEN(D)"], 0.58204, 0.54122], [["CDKN2A(D)", "PTEN(D)", "RB1(D)"], 0.59467, 0.51974], [["CDK4(A)", "PTEN(D)", "TP53"], 0.60326, 0.50541], [["MDM4(A)", "PDGFRA(A)", "PTEN"], 0.62282, 0.47349], [["MDM4(A)", "PTEN", "TP53"], 0.62456, 0.4707], [["CDK4(A)", "MDM4(A)", "TP53"], 0.65547, 0.42241], [["MDM4(A)", "PTEN(D)", "RB1(D)"], 0.67876, 0.38749], [["CDKN2A(D)", "MDM4(A)", "PTEN"], 0.68465, 0.37884], [["CDK4(A)", "PDGFRA(A)", "PTEN(D)"], 0.70878, 0.34421], [["CDKN2A(D)", "EGFR", "PDGFRA(A)"], 0.7178, 0.33157], [["CDK4(A)", "MDM4(A)", "PTEN(D)"], 0.71914, 0.3297], [["EGFR", "RB1(D)", "TP53"], 0.73762, 0.30432], [["EGFR", "MDM4(A)", "PTEN"], 0.74063, 0.30026], [["CDK4(A)", "PDGFRA(A)", "RB1(D)"], 0.75931, 0.27534], [["CDKN2A(D)", "PTEN", "PTEN(D)"], 0.77636, 0.25314], [["PDGFRA(A)", "RB1(D)", "TP53"], 0.80285, 0.21958], [["CDK4(A)", "MDM4(A)", "PTEN"], 0.82971, 0.18668], [["CDKN2A(D)", "EGFR", "PTEN(D)"], 0.83641, 0.17864], [["CDKN2A(D)", "EGFR", "MDM4(A)"], 0.83641, 0.17864], [["CDKN2A(D)", "MDM4(A)", "PTEN(D)"], 0.8707, 0.13846], [["CDKN2A(D)", "PDGFRA(A)", "PTEN(D)"], 0.8759, 0.1325], [["CDKN2A(D)", "MDM4(A)", "PDGFRA(A)"], 0.8759, 0.1325], [["EGFR", "PTEN", "RB1(D)"], 0.88503, 0.12214], [["EGFR", "PTEN", "TP53"], 0.88751, 0.11933], [["CDK4(A)", "PTEN", "TP53"], 0.88937, 0.11725], [["CDK4(A)", "PDGFRA(A)", "PTEN"], 0.91005, 0.09426], [["CDK4(A)", "PTEN", "RB1(D)"], 0.91538, 0.08841], [["PTEN", "RB1(D)", "TP53"], 0.95692, 0.04404], [["CDK4(A)", "MDM4(A)", "PDGFRA(A)"], 0.95942, 0.04142], [["CDK4(A)", "RB1(D)", "TP53"], 0.96379, 0.03688], [["PDGFRA(A)", "PTEN", "RB1(D)"], 0.9774, 0.02286], [["PDGFRA(A)", "PTEN", "TP53"], 0.98084, 0.01935], [["MDM4(A)", "PTEN", "RB1(D)"], 0.98887, 0.01119], [["CDK4(A)", "PDGFRA(A)", "TP53"], 0.9892, 0.01086]]

	# Run CoMEt MCMC using a fixed seed
	seed = 23
	mcmcArgs = [ "-m", "../example_datasets/gbm/GBM.m2",  "-o", "tmp-mcmc",
				 "-g", "../example_datasets/gbm/GBM.glst", "-N", "10000",
				 "--seed", str(seed), "-mf", "30", "-ks", "3", "3"]
	mcmcResults = RC.run(RC.get_parser().parse_args(mcmcArgs))
	os.unlink("tmp-mcmc.para.k33.10K.1.json")
	os.unlink("tmp-mcmc.sum.k33.10K.1.tsv")
		
	# Run exhaustive
	exhaustArgs = ["-m", "../example_datasets/gbm/GBM.m2", "-o", "tmp",
				   "-g", "../example_datasets/gbm/GBM.glst", "-k", "3",
				   "-mf", "30"]
	exhaustResults = RE.run(RE.get_parser().parse_args(exhaustArgs))
	exhaustResults = [ (genes, round(phi, 5), round(score, 5)) for genes, phi, score in exhaustResults]

	os.unlink("tmp-k3-exact-exhaustive.tsv")

	# Check the results
	assert(json.dumps(exhaustResults) == json.dumps(trueExhaustive))

	# Don't test the MCMC due to different PRNGs used by different Python versions
	# assert(json.dumps(mcmcResults) == json.dumps(trueMCMC))

	print "PASS"