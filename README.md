<img src="/images/enfusionNAME.png" width="450">   

## Sections

- [Introduction](#introduction)
- [Installation and running instructions for individual callers](#Installation-and-running-instructions-for-individual-callers)
- [Build the Docker image](#build-the-docker-image)
- [Run the Docker image](#run-the-docker-image)
- [Output](#output)

## Introduction

<ins>En</ins>semble <ins>Fusion</ins> (EnFusion) merges fusion output data from [Arriba](https://github.com/suhrig/arriba), [CICERO](https://github.com/stjude/CICERO), [FusionMap](http://www.arrayserver.com/wiki/index.php?title=Oshell#OmicScript_for_FusionMap), [FusionCatcher](https://github.com/ndaniel/fusioncatcher), [JAFFA](https://github.com/Oshlack/JAFFA/wiki), [MapSplice](http://www.netlab.uky.edu/p/bioinfo/MapSplice2), and [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki). 
  
To learn more about this approach, visit our BioRxiv article: [Discovery of Clinically Relevant Fusions in Pediatric Cancer](https://www.biorxiv.org/content/10.1101/2021.03.11.435013v1)

<img src="/images/biorxiv_lahaye.png" width="650">   

#### Information about the "recurrent fusion list", described in the above manuscript
This list can be found here: `SCRIPTS/R/GenePairCounts_2021-08-05.tsv` <br>
This is an internally generated list from the [The Steve and Cindy Rasmussen Institute for Genomic Medicine](https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine) at Nationwide Children's Hospital. This list contains fusion partner frequencies collected from de-identified RNA-seq data from our [Comprehensive Cancer Protocol](https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/translational-genomic-protocols) which is comprised of pediatric tumors (including central nervous system tumors, solid tumors, and hematologic malignancies). The utilization of this list allows for filtering out of recurrent and likely artifactual fusions. We will update this list biannually and will timestamp it by its release date. The below running instructions describe its use as a filtering mechanism. <br><br> 


# Installation and running instructions for individual callers

### Arriba

**Installation instructions:** Installation instructions available at the Arriba GitHub: https://github.com/suhrig/arriba  <br>
**Publication:** Uhrig, S., et al. (2021) Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Research. [manuscript link](https://genome.cshlp.org/content/31/3/448) <br>
**Running instructions:** We utilize default settings for Arriba  <br>
**Version used in publication:** v1.2.0 <br>
**Latest version validated:** v1.2.0 <br><br>

### CICERO

**Installation instructions:** Installation instructions available at the CICERO GitHub: https://github.com/stjude/CICERO <br>
**Publication:** Tian, L., et al. (2020) CICERO: a versatile method for detecting complex and diverse driver fusions using cancer RNA sequencing data. Genome Biology. [manuscript link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02043-x) <br>
**Running instructions:** We utilize default settings for CICERO <br>
**Version used in publication:** v0.3.0 <br>
**Latest version validated:** v0.3.0 <br><br>

### FusionMap

**Installation instructions:** FusionMap is a part of the Oshell toolkit and installation instructions are available at: http://www.arrayserver.com/wiki/index.php?title=Oshell#OmicScript_for_FusionMap <br>
**Publication:** Ge, H., et al. (2011) FusionMap: detecting fusion genes from next-generation sequencing data at base-pair resolution. Bioinformatics. [manuscript link](https://academic.oup.com/bioinformatics/article/27/14/1922/194689) <br>
**Running instructions:** We utilize default settings for FusionMap <br>
**Version used in publication:** v mono-2.10.9 <br>
**Latest version validated:** v mono-2.10.9 <br><br>

### FusionCatcher

**Installation instructions:** Installation instructions available at the FusionCatcher GitHub: https://github.com/ndaniel/fusioncatcher <br>
**bioRxiv Preprint:** Nicorici, D., et al. (2014) FusionCatcher – a tool for finding somatic fusion genes in paired-end RNA-sequencing data. bioRxiv. [preprint link](https://www.biorxiv.org/content/10.1101/011650v1.full.pdf+html) <br>
**Running instructions:** We utilize default settings for FusionCatcher <br>
**Version used in publication:** v0.99.7c <br>
**Latest version validated:** v0.99.7c <br><br>

### JAFFA

**Installation instructions:** Installation instructions available at the JAFFA GitHub: https://github.com/Oshlack/JAFFA/wiki/Download <br>
**Publication:** Wang, K., et al. (2010) JAFFA: High sensitivity transcriptome-focused fusion gene detection. Genome Medicine. [manuscript link](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0167-x)<br>
**Running instructions:** We utilize default settings for JAFFA direct <br>
**Version used in publication:** direct v1.09 <br>
**Latest version validated:** direct v1.09 <br><br>

### MapSplice

**Installation instructions:** Installation instructions available at the MapSplice GitHub: https://github.com/LiuBioinfo/MapSplice <br>
**Publication:** Davidson, N.M., et al. (2015) MapSplice: Accurate mapping of RNA-seq reads for splice junction discovery. Nucleic Acids Research. [manuscript link](https://academic.oup.com/nar/article/38/18/e178/1068935)<br>
**Running instructions:** We utilize default settings for MapSplice <br>
**Version used in publication:** v2.2.1 <br>
**Latest version validated:** v2.2.1 <br><br>

### STAR-Fusion

**Installation instructions:** Installation instructions available at the STAR-Fusion GitHub: https://github.com/STAR-Fusion/STAR-Fusion/wiki <br>
**Publication:** Haas, B.J., et al. (2019) Accuracy assessment of fusion transcript detection via read-mapping and de novo fusion transcript assembly-based methods. Genome Biology. [manuscript link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1842-9) <br>
**Running instructions:** STAR-Fusion parameters were altered to reduce the stringency setting for the fusion fragments per million reads (FFPM) to 0.02. <br> `--min_FFPM 0.02` <br>
**Version used in publication:** v1.6.0 <br>
**Latest version validated:** v1.6.0 <br><br>
<br> <br>  



# Instructions to build and run EnFusion overlap Docker Image 

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

#### Before building the image, if the user would like to use a `known fusion list`, please upload this and save it as `known_fusion_list.txt` in the `SCRIPTS` directory. Fusion partners listed in this file are not filtered out, regardless of level of support for the fusion. The `known fusion list` file should be a txt file where each fusion pair (gene names ordered alphabetically and genes separrated by a +) is on a single line. No header should be added. Follow example `known_fusion_list.txt` file (located in `SCRIPTS`)

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

#### Please note that for the overlap script to run properly, there is an expected file structure hierarchy (refer to how test data is stored in the `test_data` directory for an example, however, as long as the below files are included within the mounted volume to the Docker, the overlap script will recursively search all directories under the `input_location` for fusion detection result files. The default search looks for the following expected file names (which are automatically generated by the individual fusion detection algorithms):

arriba        = `fusions.tsv`  
cicero        = `annotated.fusion.txt`   
fusionCatcher = `final-list_candidate-fusion-genes.txt`   
fusionMap     = `FusionDetection.FusionReport.Table.txt`   
jaffa         = `jaffa_results.csv`   
mapSplice     = `fusions_well_annotated.txt`   
starFusion    = `star-fusion.fusion_predictions.abridged.tsv`   

## Run the Docker image

#### We must first mount our host directory that contains the fusion detection results as a volume to our Docker container:

##### To do this this with the test data, save test the data directory to your local machine and then mount it as a volume:

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

#### Note that the `-f` argument utilizes the "recurrent fusion list" which has been generated from an internal RNA-seq cohort and is used to identify recurrent and likely artifactual fusions. This list can be found here `SCRIPTS/R/GenePairCounts_2021-08-05.tsv` and will be updated biannually and timestamped by date of collection. Please note that if you would like to prevent any filtering, you can set this argument to `0`.

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
  
