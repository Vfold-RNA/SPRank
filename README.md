# SPRank - A knowledge-based scoring function for RNA-ligand pose prediction and virtual screening

## Platform Requirements (Tested)
The following are tested system settings.

For compiling and running */bin/sprank*
* GNU/Linux x86_64 (Ubuntu 22.04 LTS)
* GNU Make 4.3
* gcc/g++ (version 11.4.0)

For running */bin/check_atom_order*
* Python 3.9.21
* NumPy 2.0.2

For running */bin/ambertools_prepare_rec* and */bin/ambertools_prepare_cpd*
* AmberTools22

For running random forest model attached in the **Releases**
* scikit-learn 1.5.0
* Polars 1.17.1
* NumPy 2.0.2
* Matplotlib 3.9.4

## Installation

#### Compile and install SPRank
Clone this repository on your local machine and run setup script
```
git clone https://github.com/Vfold-RNA/SPRank.git ${HOME}/SPRank
```
compile the code
```
cd ${HOME}/SPRank && make
```
put the following environment variable to your *.bashrc*
```
echo "export SPRANK_HOME=${HOME}/SPRank" >> ${HOME}/.bashrc
```
and source it
```
source ${HOME}/.bashrc
```

## Using SPRank and SPRank-RF

#### Check SPRank options
```
${SPRANK_HOME}/bin/sprank -h
```

#### or check SPRank-RF options
```
${SPRANK_HOME}/bin/sprank-rf -h
```

#### Run SPRank and SPRank-RF for example cases
```
cd ${SPRANK_HOME}/example && chmod +x ./run_example && ./run_example
```
The predicted scores will be saved in the corresponding folders as 
`score.dat` for *sprank* and `score_random_forest.dat` for *sprank-rf*.
By default, the script does not run AmberTools22 to prepare the input receptor and compound.
You can remove the comments in `./run_example` to prepare the input files using AmberTools22.

#### Run conversion scripts for example cases
```
cd ${SPRANK_HOME}/example && chmod +x ./run_conversion_test && ./run_conversion_test
```
This script will run `convert_rdock_pose` and `convert_vina_pose` to convert 
rDock and AutoDock Vina generated poses to mol2 format compatible with `sprank` and `sprank-rf`. 
After running the script, the predicted scores will be saved in the corresponding folders as 
`*_pose_score.dat` for `sprank` and `*_pose_random_forest_score.dat` for `sprank-rf`, 
where * will be rdock and vina. 
By default, the script does not run AmberTools22 to prepare the input receptor and compound.
You can remove the comments in `./run_conversion_test` to prepare the input files using AmberTools22.


## SPRank command line arguments
```
-r <receptor>         # path to target RNA (in amber mol2 format, must contain hydrogens)
-c <target compound>  # path to target compound (in amber mol2 format, must contain hydrogens 
                        and bond table (i.e., "@<TRIPOS>BOND" record))
-p <compound poses>   # path to poses sampled by docking software,
                        to be scored by SPRank (in mol2 format,
                        the order of the heavy atoms should be same as the target compound)
```

## SPRank-RF command line arguments
```
-r <receptor>         # path to target RNA (in amber mol2 format, must contain hydrogens)
-c <target compound>  # path to target compound (in amber mol2 format, must contain hydrogens 
                        and bond table (i.e., "@<TRIPOS>BOND" record))
-p <compound poses>   # path to poses sampled by docking software,
                        to be scored by SPRank (in mol2 format,
                        the order of the heavy atoms should be same as the target compound)
-o <output>           # path to save the RandomForest predicted scores
```

## Download data
The training set, pose sets, affinity sets, random forest model, amber atom types, potentials, 
HIV-1 TAR ensemble and compound library can be downloaded from the **Releases** or through the following commands:
```
mkdir -p ${SPRANK_HOME}/data/
```
```
for name in "checksum.txt" "training-set.tar.gz" "pose-sets.tar.gz" "affinity-sets.tar.gz" "random-forest.tar.gz" "amber-types.tar.gz" "potentials.tar.gz" "HIV-1-TAR.tar.gz"
do
    wget https://github.com/Vfold-RNA/SPRank/releases/download/data/${name} -O ${SPRANK_HOME}/data/${name}
done
```
Check the integrity of the files:
```
cd ${SPRANK_HOME}/data/
```
```
sha256sum --check checksum.txt
```

## References

[1] Zhou Y, Jiang YW, Chen SJ.
**SPRank - A Knowledge-Based Scoring Function for RNA-Ligand Pose Prediction and Virtual Screening.**
*Journal of Chemical Theory and Computation.*
2024 Aug. doi: 10.1021/acs.jctc.4c00681.
