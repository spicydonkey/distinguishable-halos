# DISTINGUISHABLE HALO ANALYSIS
2017.06.21 - David Shin

## TODO
- [] comment functions
- [] unify function io behaviour:
    - work by single-shot: input is an array of Nx3 (for zxy), except for g2 analysis, etc. which require shot-to-shot identification

## GETTING STARTED
1. Configure the analysis by creating a Matlab script using config_090617_1" as an example template
	* turn the correlation analysis flag (+ savedata) off when getting started with a new dataset or experimenting with halo captures
2. Open Matlab and run:
	[halo_k,corr,efit,halo,txy,fout]=run_dist_halo(PATH_TO_CONFIG);

## METHODS
### g2
#### x-state (2 component); BB; Cartesian
1. Correlated pairs
	1.1. get all same-shot, unordered, unique, and x-state pairs 
	1.2. evaluate the delta (dk) for BB correlation; BB condition requires simple sum of a pair of vectors, giving dk:=k_{1}+k_{2}
	1.3. histogram the deltas (3D vector; 0 is perfect BB) --> G2_corr(dk)
	1.4. normalise histogram wrt number of pairs
2. Uncorrelated pairs
	2.1. get all diff-shot, unordered, unique, and x-state pairs
	2.2. evaluate the delta (dk) for BB correlation (same as 1.2)
	2.3. histogram the deltas --> G2_uncorr(dk)
	2.4. normalise histogram wrt number of pairs
3. Evaluate normalised g2
	3.1. g2(dk):=G2_corr(dk)/G2_uncorr(dk)

## BUGS
* 2017.06.21 - FIXED - formula for g2 normalisation incorrect by a factor of 2 - g2 plateauxing at 1/2 was a danger sign 
