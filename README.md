<img src="/images/enfusionNAME.png" width="450">   

## Sections

- [Introduction](#introduction)
- [Build the Docker image](#build-the-docker-image)
- [Run the Docker image](#run-the-docker-image)
- [Output](#output)

## Introduction

<ins>En</ins>semble <ins>Fusion</ins> (EnFusion) merges fusion output data from [Arriba](https://github.com/suhrig/arriba), [CICERO](https://github.com/stjude/CICERO), [FusionMap](http://www.arrayserver.com/wiki/index.php?title=Oshell#OmicScript_for_FusionMap), [FusionCatcher](https://github.com/ndaniel/fusioncatcher), [JAFFA](https://github.com/Oshlack/JAFFA/wiki), [MapSplice](http://www.netlab.uky.edu/p/bioinfo/MapSplice2), and [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki). 
  
To learn more about this approach, visit our BioRxiv article: [Discovery of Clinically Relevant Fusions in Pediatric Cancer](https://www.biorxiv.org/content/10.1101/2021.03.11.435013v1)

<img src="/images/biorxiv_lahaye.png" width="650">   

## Build the Docker image

We tested the instructions on Mac. They may work on Linux or Windows with or
without variation.

### 1. Install requirements

- git
- Docker

### 2. Download EnFusion from GitHub

Use `git` to clone the repo or navigate to the [EnFusion](https://github.com/nch-igm/EnFusion) GitHub page, download, and unzip the code.

```bash
git clone https://github.com/nch-igm/EnFusion.git
```

### 3. Build the EnFusion Docker image

Before building the image, if the user would like to use a `known fusion list`, please upload this and save it within the `SCRIPTS` directory. This file should be a txt file where each fusion paid (separrated by a --) is on a single line. No header should be added. Follow example 

Navigate to the EnFusion directory and run Docker build.

```bash
cd EnFusion
docker build . -t enfusion
```

#### Test the Docker image

Get a help message from the entrypoint.

```bash
docker run enfusion -h
```

Useful flags

- `--rm` deletes the container after is stops running. You can use the command `docker container ls --all` to view stopped containers that have not been deleted.
- `-it` allows for an interactive session.
- `--entrypoint "/bin/bash"` overwrites the entrypoint with the bash binary.

#### Please note that for the overlap script to run properly, there is an expected file structure hierarchy, however, as long as the below files are included within the mounted volume to the Docker, the overlap script will recursively search all directories under the `input_location` for fusion detection result files. The default search looks for the following expected file names:

arriba        = `fusions.tsv`  
cicero        = `annotated.fusion.txt`   
fusionCatcher = `final-list_candidate-fusion-genes.txt`   
fusionMap     = `FusionDetection.FusionReport.Table.txt`   
jaffa         = `jaffa_results.csv`   
mapSplice     = `fusions_well_annotated.txt`   
starFusion    = `star-fusion.fusion_predictions.abridged.tsv`   

## Run the Docker image

#### We must first mount our host directory that contains the fusion detection results as a volume to our Docker container:

##### To practice with this test data, save test data directory to local machine and then mount it as a volume:

Test data is located here: `EnFusion/test_data/test`

The test data contains output from 5 callers, and upon downloading this data, the directory structure should look like this:

``` 
├── test_data   
    ├── samples   
    └── test   
        ├── fusioncatcher   
        │   └── final-list_candidate-fusion-genes.txt   
        ├── fusionmap   
        │   └── results   
        │       └── FusionDetection.FusionReport.Table.txt   
        ├── jaffa   
        │   └── jaffa_results.csv   
        ├── mapsplice   
        │   └── fusions_well_annotated.txt   
        └── starfusion   
            └── star-fusion.fusion_predictions.abridged.tsv  
```
In this folder we also have a ```samples``` file which lists all samples included in the test_data directory. The sample listed in ```samples``` matches to the name of the other subdirectory ```test```. More than 1 sample can be in the ```samples``` directory, as long as those same samples are included as directories (samples). In this case above, the Patient ID is ```test_data``` and the Sample ID is ```test```.

#### Use Docker run to invoke EnFusion
```
docker run -v /~localpath/EnFusion/test_data:/SCRIPTS/test_data enfusion -o SCRIPTS/test_data/test -s test_data -p test -f 0.2
```
The `-v` flag will mount a host directory as a data volume to the docker container.   
The first part (before the `:`) `/~localpath/EnFusion/test_data` needs to be the absolute path to where you have your data.  
The second part (after the `:`)`/SCRIPTS/test_data` is where the data will be written to within the Docker container.
  
The following arguments are passed to the overlap script:  
`-o` output_location  
`-s` sample_ID  
`-p` patient_ID  
`-f` frequency cutoff to use (default is 0.10 ~ 10%, enter a decimal value here or leave blank for 10% cutoff)

When running the test data you will see the following on your screen:
  
 ```
SCRIPTS/test_data/test
test_data
/working_dir
cp: cannot stat '/SCRIPTS/test_data/test/arriba/fusions.tsv': No such file or directory
cp: cannot stat '/SCRIPTS/test_data/test/cicero/annotated.fusion.txt': No such file or directory
[1] "/working_dir"
[1] "/SCRIPTS/R"
Assembing a list of files to merge...
[1] "starFusion"    "fusionMap"     "fusionCatcher" "jaffa"        
[5] "mapSplice"    
```
Because we only have 5 outputs in the test data, you will see an error messages about the missing data from ```arriba``` and ```cicero```. This does not affect the ability of the overlap analysis to run, and is instead a note to let the user know that output from all 7 tools was not provided. For the overlap algorithm to run, output is only required from 2 callers (at minimum). The ```list of files to merge``` lets the user know which outputs were identified (of the 7 possible).
   
Next you will see a list of all unique fusions identified:
```
 [1] "KIF5B+RET"      "ALK+EML4"       "ETV6+NTRK3"     "LMNA+NTRK1"    
 [5] "FGFR3+TACC3"    "NCOA4+RET"      "NTRK1+TPM3"     "PAX8+PPARG"    
 [9] "BRAF+SLC45A3"   "BAIAP2L1+FGFR3" "ROS1+SLC34A2"   "CD74+ROS1"     
[13] "ERG+TMPRSS2"    "EGFR+SEPTIN14"  "PLEC+TSPAN4"    "NEK1+SNX25"    
[17] "TASOR+UBE2K"    "SNX29+TXNDC11"  "SAMD5+SASH1"    "CBX3+CCDC32"   
[21] "NAIP+OCLN"      "NCOR2+UBC"      "STK3+VPS13B"  
```
Then, printed to screen are the overlap results:
```
[1] "# Sample: test_data\n# NumToolsAggregated: \t5\n# - starFusionCalls = \t37\n# - fusionMapCalls = \t32\n# - fusionCatcherCalls = \t355\n# - jaffaCalls = \t491\n# - mapSpliceCalls = \t26\n# filtered_overlap = \t14"
   UnorderedFusion   OrderedFusion KnownFusion NumTools   GenePairFrequency
1   BAIAP2L1+FGFR3 FGFR3>>BAIAP2L1         yes        5  0.0103675777568332
3        CD74+ROS1      CD74>>ROS1         yes        5 0.00942507068803016
6       ETV6+NTRK3     ETV6>>NTRK3         yes        5  0.0141376060320452
7        KIF5B+RET      KIF5B>>RET         yes        5  0.0103675777568332
8       LMNA+NTRK1     LMNA>>NTRK1         yes        5 0.00754005655042412
9        NCOA4+RET      NCOA4>>RET         yes        5  0.0113100848256362
10      NTRK1+TPM3     TPM3>>NTRK1         yes        5  0.0103675777568332
11      PAX8+PPARG     PAX8>>PPARG         yes        5  0.0131950989632422
12    ROS1+SLC34A2   SLC34A2>>ROS1         yes        5 0.00942507068803016
2     BRAF+SLC45A3   SLC45A3>>BRAF         yes        4 0.00848256361922714
5      ERG+TMPRSS2    TMPRSS2>>ERG         yes        3 0.00188501413760603
4    EGFR+SEPTIN14  EGFR>>SEPTIN14         yes        2                   0
```
The above print to screen can be quite long.
  
You will also see any errors or warning printed to screen.

## Output  
  
### Upon completion of the script, 4 files will be output and stored to your output_location
  
`overlap_$sample_name.tsv`: full list of all intersecting fusion across callers, no filtering of any kind
  
`filtered_overlap_knownfusionlist_3callers_$sample_name.tsv`: full list of all fusions passing filter, each line for each fusion is associated with the contributing caller  
  
`collapse_filtered_overlap_knownfusionlist_3callers_$sample_name.tsv`: *collapsed* list, where each row is associated with a single fusion (not a single caller) this output does not include all individual caller information  
  
`Singleton_KnownFusions_$sample_name.tsv`: this output includes any fusions on the `known fusion list` that were only identified by a single caller
  
    
    
#### For more information about the individual scripts and to run the scripts without Docker, please see this GitHub repo: https://github.com/nch-igm/nch-igm-ensemble-fusion-detection
