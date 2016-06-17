# jInv Examples used in Manuscript "jInv - a flexible Julia package for PDE parameter estimation" (2016)

Julia files used to produce results reported in the above manuscript. 

## Full Waveform Inversion


## Joint Reconstruction 

## Weak Scaling Experiments


### DC Resistivity on Cloud Engine
For testing the scalability of the DC resistivity we used the file `runWeakScalingDivSigGrad.jl` which can be run from the commandline using 
```
usage: runWeakScalingDivSigGrad.jl [--n1 N1] [--n2 N2] [--n3 N3]
                        [--srcSpacing SRCSPACING]
                        [--nthreads NTHREADS] [--out OUT]
                        [--solver SOLVER] [-h]
```
We used Amazon EC2 cloud using a generated AMI and initialized 50 instances of which 49 got the tag "workers". Then, we ran a shell script similar to:
```
#!/bin/bash
for np in {1..50}
do
	# get machine file
	ec2-describe-instances --filter "tag:Name=workers" | grep running | awk '{ print $5; }' > mfile
	for s in 3 2 1 # use all solvers
	do
		for k in {1..5} # repeat five times
		do
			julia --machinefile mfile runWeakScalingDivSigGrad.jl --nthreads=2 --n1 48 --n2 48 --n3 24 --srcSpacing 5 --out weakScalingMUMPSAmazon.csv --solver $s
		done
	done
	# terminate one instances and remove from machine file
	ec2-describe-instances  --filter "tag:Name=workers"| grep running | awk '{ print $2; }' > idfile
	tail -1 idfile | xargs -I % sh -c  'ec2-terminate-instances  %'
done

```


## Strong Scaling Experiments