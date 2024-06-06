# SPRank - A knowledge-based scoring function for RNA-ligand pose prediction and virtual screening

## Platform Requirements (Tested)
The following are tested system settings.

For compiling and running */bin/sprank*
* GNU/Linux x86_64 (Ubuntu 22.04.3 LTS kernel 5.15.0-91-generic)
* GNU Make 4.3
* gcc/g++ (version 11.4.0)

For running */util/check_atom_order*
* Python 3.11.5
* NumPy 1.26.0

For running */util/ambertools_prepare_rec* and */util/ambertools_prepare_cpd*
* AmberTools22

For running random forest model attached in the **Releases**
* scikit-learn 1.5.0
* Polars 0.20.4
* NumPy 1.26.0
* Matplotlib 3.7.2

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

## Using SPRank

#### Check SPRank options
```
${SPRANK_HOME}/bin/sprank -h
```

#### Run SPRank for example cases
```
cd ${HOME}/SPRank/example && chmod +x ./run_example && ./run_example
```
The predicted scores will be saved in the corresponding folders as `score.dat`.
By default, the script does not run AmberTools22 to prepare the input receptor and compound.
You can remove the comments in `./run_example` to prepare the input files using AmberTools22.


## SPRank command line arguments
```
-r <receptor>         # path to target RNA (in amber mol2 format, must contain hydrogens)
-c <target compound>  # path to target compound (in amber mol2 format, must contain hydrogens 
                        and bond table (i.e., "@<TRIPOS>BOND" record))
-p <compound poses>   # path to poses sampled by docking software,
                        to be scored by SPRank (in mol2 format,
                        the order of the heavy atoms should be same as the target compound)
```

## Download data
The training set, pose sets, affinity sets, random forest model, amber atom types, potentials, 
HIV-1 TAR ensemble and compound library can be downloaded from the **Releases** or through the following commands:
```
mkdir -p ${HOME}/SPRank/data/
```
```
for name in "checksum.txt" "training-set.tar.gz" "pose-sets.tar.gz" "affinity-sets.tar.gz" "random-forest.tar.gz" "amber-types.tar.gz" "potentials.tar.gz" "HIV-1-TAR.tar.gz"
do
    wget https://github.com/Vfold-RNA/SPRank/releases/download/data/${name} -O ${HOME}/SPRank/data/${name}
done
```
Check the integrity of the files:
```
cd ${HOME}/SPRank/data/
```
```
sha256sum --check checksum.txt
```

## References

[1] Yuanzhe Zhou, Yangwei Jiang, and Shi-Jie Chen. 
SPRank - A knowledge-based scoring function for RNA-ligand pose prediction and virtual screening. 
2024. *submitted*
