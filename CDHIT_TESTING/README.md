# Clustering Repository
The goal for this folder is to test and cluster proteins based on the amount of shared
identity for all the proteins in a fasta. The overall goal is to save commands for 
using CS-Hit and create an adapter for CD-HIT. 

# Obtaining test proteins
* The following curl command downloads proteins from uniprot
```
curl "https://www.uniprot.org/uniprot/?query=Phage&format=fasta&limit=10000&sort=score" > protein_test.fa;
```

# Using CD-HIT
* General Use (These are the most useful commands)
```
cd-hit -i <input fasta file> -o <Output files> -c <cutoff> -T <Threads>
```

* Specific use
```
cd-hit -i protein_test.fa -o output.txt -c 0.65 -T 3
```

# using JPRED (API - yet finicky)

```
python3 -m pip install jpredapi
```

1. The command to run is:
```
jpredapi.submit(mode="single", user_format="raw", seq="MQVWPIEGIKKFETLSYLPP")
```

2. Retrieve the JOBID
Currently its printed to STDOUT, so it's not clear how to do this?
 One idea is to use a giant fasta or to retrieve STDOUT. 

3. To get the results
```
jpredapi.get_results(jobid="jp_Is5LeY6", results_dir_path="jpred_sspred/results")
```

3. Then the file is located in:
```
cat jpred_sspred/results/jp_Is5LeY6/jp_Is5LeY6.jnet
```

OUtput: 
```
jnetpred:-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,H,H,H,H,H,H,H,H,H,H,H,H,-,-,-,-,-,-,-,-,-,-,H,H,H,H,H,H,-,-,-,-,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,-,-,-,H,H,H,H,H,H,H,-,-,-,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,-,-,-,-,-,-,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,-,-,-,-,-,-,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,-,-,-,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,-,-,-,-,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,-,-,-,-,-,-,-,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,E,E,E,E,-,-,-,-,-,-,-,-,-,-,-,-,
JNETCONF:8,9,8,8,7,7,7,7,7,7,7,7,7,7,7,7,6,4,2,8,9,9,9,9,9,9,9,8,6,0,3,6,7,8,8,7,5,2,4,0,1,2,3,1,0,3,0,6,6,4,6,8,9,9,9,9,9,9,9,9,9,9,8,7,4,1,2,5,5,7,7,5,1,2,2,3,7,8,8,8,6,1,4,6,1,8,9,9,9,9,9,9,9,9,9,9,9,9,8,8,4,8,9,9,9,9,9,9,9,9,8,6,1,5,8,7,4,2,7,8,9,9,9,9,9,9,9,9,4,9,9,9,9,9,9,9,9,8,7,3,0,3,3,6,7,6,4,7,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,8,7,4,4,7,1,7,8,9,9,9,9,9,9,9,9,9,8,8,6,4,7,8,8,8,9,9,9,9,9,9,9,9,9,8,5,1,6,7,8,1,1,2,3,2,7,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,8,7,0,1,6,8,7,6,2,7,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,8,5,1,5,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,5,3,3,3,5,6,6,7,7,7,7,7,7,6,5,1,1,0,1,5,6,3,1,6,6,0,2,3,3,3,6,6,7,7,7,8,8,9,
JNETSOL25:-,-
```

# using PORTER5

0. get the tool
```
git clone https://github.com/mircare/Porter5.git
```

1. get Uniprot20
```
wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/uniprot20_2016_02.tgz
```
* untar 
```
tar -zxf uniprot20_2016_02.tgz
```

2. HHBlitz
```
conda install -c bioconda hhsuite
```

3. Run and add path to Uniprot20 and command for hhsuite is `hhblits`
* to set paths
```
 python Porter5/Porter5.py -i test.fa --fast --setup
```
add:
```
Please insert the absolute path to psiblast (e.g., /home/username/psiblast): /Users/dreyceyalbin/miniconda/envs/DebruijnExtend/bin/psiblast
Please insert the absolute path to uniref90 (e.g., /home/username/UniProt/uniref90.fasta): /Users/dreyceyalbin/Dropbox/Fall2020-classes/Algorithms/project/DebruijnExtend/CDHIT_TESTING/data_for_tools/uniprot20_2016_02
Please insert the call to HHblits (e.g., hhblits): /Users/dreyceyalbin/miniconda/envs/DebruijnExtend/bin/hhblits
Please insert the absolute path to uniprot20 - DATABASE NAME INCLUDED (e.g., /home/username/uniprot20_2016_02/uniprot20_2016_02): /Users/dreyceyalbin/Dropbox/Fall2020-classes/Algorithms/project/DebruijnExtend/CDHIT_TESTING/data_for_tools/uniprot20_2016_02
```

* To run the tool: 
```
python Porter5/Porter5.py -i test.fa --fast 
```
Output:
```
#	AA	SS	Helix	Sheet	Coil
1	M	C	0.0002	0.0001	0.9997
2	T	C	0.0341	0.0162	0.9496
3	A	C	0.1036	0.0125	0.8839
4	P	C	0.1309	0.016	0.8531
5	T	C	0.176	0.055	0.7689
6	L	C	0.1752	0.1433	0.6815
```

# SPIDER2

0. get SPIDER2
```
wget https://servers.sparks-lab.org/downloads/SPIDER2_local.tgz;
tar -zxvf SPIDER2_local.tgz;
```

1. Connect psi blast
* download from conda
```
conda install -c bioconda blast
```
* modify BASH
```
export blastpgp=psiblast
```

2. Download NR
* command for downloading NR data
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
```
* decompress file
```
gzip -d nr.gz
```

* update BASH
```
export NR='${PWD}/data_for_tools/nr.gz'
```

2.5. Change the SPIDER script to use `python2`
```
$xdir/pred_pssm.py $pro1.pssm 
```
to
```
python2 $xdir/pred_pssm.py $pro1.pssm
```

3. Run an example
```
bash ../SPIDER2/misc/run_local.sh example_protein.seq
```

4. Output saved as <name>.spd3
```
cat protein_test.spd3
```

Output:
```
#	AA	SS	ASA	Phi	Psi	Theta(i-1=>i+1)	Tau(i-2=>i+1)	P(C)	P(E)	P(H)
1	M	C	127.1	-100.9	 139.3	 121.6	-143.9	0.973	0.021	0.003
2	P	C	107.9	 -66.7	 140.0	 109.9	-116.3	0.960	0.030	0.010
3	V	C	 62.1	-105.2	 133.7	 118.3	-154.9	0.897	0.096	0.008
```

# S4Pred (WORKING!)

1. Obtain the git repository.
```
git clone https://github.com/psipred/s4pred
```

2. Next, download the pretrained weights.
```
cd s4pred/; mkdir weights; cd weights;
wget http://bioinfadmin.cs.ucl.ac.uk/downloads/s4pred/weights_1.pt;
wget http://bioinfadmin.cs.ucl.ac.uk/downloads/s4pred/weights_2.pt;
wget http://bioinfadmin.cs.ucl.ac.uk/downloads/s4pred/weights_3.pt;
wget http://bioinfadmin.cs.ucl.ac.uk/downloads/s4pred/weights_4.pt;
wget http://bioinfadmin.cs.ucl.ac.uk/downloads/s4pred/weights_5.pt;
```

3. Using the tool.
```
python ./s4pred/run_model.py --outfmt fas test.fa
```

* outout for (`--outfmt fas`)
```
>1QYS_1|Chain A|TOP7|Computationally Designed Sequence
MGDIQVQVNIDDNGKNFDYTYTVTTESELQKVLNELMDYIKKQGAKRVRISITARTKKEAEKFAAILIKVFAELGYNDINVTFDGDTVTVEGQL
CCCEEEEEEECCCCCEEEEEEEECCHHHHHHHHHHHHHHHHHCCCCEEEEEEEECCHHHHHHHHHHHHHHHHHCCCCEEEEEEECCEEEEEEEC
```

* outout for (`--outfmt ss2`)
```
# PSIPRED VFORMAT (S4PRED V1.0)

   1 M C   0.997  0.000  0.002
   2 G C   0.984  0.000  0.016
   3 D C   0.812  0.001  0.187
   4 I E   0.351  0.000  0.648
   5 Q E   0.084  0.000  0.915
   6 V E   0.042  0.000  0.958
   7 Q E   0.021  0.000  0.979
   8 V E   0.019  0.000  0.981
```

# DeepCNF (still working on)

1. Use the git downlaod
```
git clone https://github.com/realbigws/Predict_Property
```

2. Compile
```
cd source_code/
    make
cd ../
```

*NOTE*: If you run into problems with openmp, do the following:

* update compiler using HomeBrew
```
brew install libomp
```

* use updated compilter in the Make file for `DeepCNF_Pred`
```
CXX = g++
```
to
```
CXX = /usr/local/bin/g++-10
```