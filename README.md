<img src="/images/EnFusion_name.png" width="450">   

## Sections

- [Introduction](#introduction)
- [Build the Docker image](#build-the-docker-image)
- [Run the Docker image](#run-the-docker-image)

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

Navigate to the EnFusion directory and run Docker build.

```bash
cd EnFusion
docker build . -t enfusion
```

#### Test the Docker image

Get a help message from the entrypoint.

```bash
docker run enfusion --help
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

### 4. Run the Docker image

#### We must first mount our host directory that contains the fusion detection results as a volume to our Docker container:

```bash
docker run -v /Users/sdl002/EnFusion/test_data:/SCRIPTS/test_data enfusion -o SCRIPTS/test_data/test -s test_data -p test 
```

The `-v` flag will mount a host directory as a data volume to the docker container
The following arguments are passed to the overlap script:  
`-o` output location  
`-s` sample ID  
`-p` patient ID  

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
[1] "# Sample: test_data\n# NumToolsAggregated: \t5\n# - starFusionCalls = \t37\n# - fusionMapCalls = \t32\n# - fusionCatcherCalls = \t355\n# - jaffaCalls = \t491\n# - mapSpliceCalls = \t26\n# filtered_overlap = \t13"
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
4    EGFR+SEPTIN14  EGFR>>SEPTIN14          no        2                   0
```
The above print to screen can be quite long.
  
You will also see any errors or warning printed to screen.

### 5. Output  
  
    

# TODO: nch-igm-ensemble-fusion-detection

#### In order for the overlap script to run properly, there is an expected file structure hierarchy, please see
#### the Example below for assemble_results or view example data in "test_data" directory

### kickoff_overlap.sh
Run 'kickoff_overlap.sh' to run the necessary R scripts (described in detail below)

```bash
overlap/kickoff_overlap.sh -h
```
The result:
```txt
[USAGE]: This script kicks off a script that merges fusion detection results from many tools.

-h * [None] Print this help message.
-p * [path1,path2,...] paths  where fusion results are where the script will be kicked off.
-s * [samples_file1,samples_file2,...] Files paths that contain the lists of samples to kickoff scripts for. The number of samples files provided must be equal to the number of paths, and they must be listed in the same respective order. If not specifically added, the script will search through sample names that are listed in the samples file within the path (-p) given[Default: samples,samples,...].
```
Note that -s is not required, and if not added the script will automatically run on all samples listed in your samples file


#### Try kickoff_overlap.sh with test data:

```bash
overlap/kickoff_overlap.sh -p /test_data
```

Results (will show you which callers it found output for):
```txt
Test
WARNING: ignoring environment value of R_HOME
Assembling a list of files to merge...
[1] "starFusion"    "fusionMap"     "fusionCatcher" "jaffa"        
[5] "mapSplice"     "dragen"       
```

#### Explanation of and option to run R code directly:

Run `assemble_results.R` to merge results and generate the overlap report. See the help message for `assemble_results.R` by running

```bash
R/assemble_results.R -h
```

or

```bash
Rscript R/assemble_results.R -h
```

The result

```txt
Usage: assemble_results.R --sample s1

assemble_results.R reads a set of fusion detection result files and produces an aggregated report with the overlapping fusion events predicted in the input files. The program searches recursively from the current directory to find known files. The input files can also be specified manually.

Options:
	--sample=SAMPLE
		Name of the sample. (required)

	--baseDir=BASEDIR
		Base directory to search for input files. [default = .]

	--outReport=OUTREPORT
		Location of the output report. [default = overlap_$sample.tsv]

	--collapseoutReport=COLLAPSEOUTREPORT
		Location of the output report. [default = collapsed_3callers_$sample.tsv]

	--foutReport=FOUTREPORT
		Location of the output report. [default = filtered_overlap_2callers_$sample.tsv]

	--foutReport3=FOUTREPORT3
		Location of the output report. [default = filtered_overlap_3callers_$sample.tsv]

	--outSingleton=OUTSINGLETON
		Location of the Singleton output report. [default = Singleton_KnownFusions_$sample.tsv]

	--outBreakpoints=OUTBREAKPOINTS
		Location of the breakpoint report. [default = breakpoints_$sample.tsv]

	--dragen=DRAGEN
		Path to the dragen results file. [default = search for file named like 'DRAGEN.fusion_candidates.final']

	--fusionCatcher=FUSIONCATCHER
		Path to the fusionCatcher results file. [default = search for file named like 'final-list_candidate-fusion-genes.txt']

	--fusionMap=FUSIONMAP
		Path to the fusionMap results file. [default = search for file named like 'FusionDetection.FusionReport.Table.txt']

	--jaffa=JAFFA
		Path to the jaffa results file. [default = search for file named like 'jaffa_results.csv']

	--mapSplice=MAPSPLICE
		Path to the mapSplice results file. [default = search for file named like 'fusions_well_annotated.txt']

	--soapFuse=SOAPFUSE
		Path to the soapFuse results file. [default = search for file named like '.final.Fusion.specific.for.genes']

	--starFusion=STARFUSION
		Path to the starFusion results file. [default = search for file named like 'star-fusion.fusion_predictions.abridged.tsv']

	--tophatFusion=TOPHATFUSION
		Path to the tophatFusion results file. [default = search for file named like 'result.txt']. potential_fusion.txt must be in the same directory as result.txt for the import to work.

	--arriba=ARRIBA
		Specific location of the arriba results file (fusions.tsv).

	--cicero=CICERO
		Specific location of the cicero results file (annotated.fusion.txt).

	-h, --help
		Show this help message and exit```

### `--sample`

`--sample` is the only required argument for `assemble_results.R`.

```bash
R/assemble_results.R --samples s1
```

`--sample` is used to insert a comment into the top of the output files

```txt
# Sample: s1
# NumToolsAggregated: 5
# - starFusionCalls = 20
# - fusionCatcherCalls = 552
# - jaffaCalls = 948
# - mapSpliceCalls = 23
# - dragenCalls = 9
```

and name the output files

```txt
--outReport      	= overlap_$sample.tsv
--collapsedoutReport	= collapsed_3callers_$sample.tsv
--foutReport     	= filtered_overlap_2callers_$sample.tsv
--foutReport3    	= filtered_overlap_3callers_$sample.tsv
--outSingleton		= Singleton_KnownFusions_$sample.tsv
--outBreakpoints 	= breakpoints_$sample.tsv
```

You can override the default output file names with the `--*out*` parameters.

### `--baseDir`

The `--baseDir` argument sets the location of the input data. By default, `--baseDir` is set to the current directory `.`, but you can override the default location

```bash
R/assemble_results.R --sample s1 --baseDir fusion_results
```

`kickoff_overlap.R` recursively searches all directories under the `baseDir` for fusion detection results files. The default search looks for

```txt
--arriba        = fusions.tsv
--cicero        = annotated.fusion.txt
--dragen        = DRAGEN.fusion_candidates.final
--fusionCatcher = final-list_candidate-fusion-genes.txt
--fusionMap     = FusionDetection.FusionReport.Table.txt
--jaffa         = jaffa_results.csv
--mapSplice     = fusions_well_annotated.txt
--soapFuse      = .final.Fusion.specific.for.genes
--starFusion    = star-fusion.fusion_predictions.abridged.tsv
--tophatFusion  = result.txt (potential_fusion.txt must be in the same directory as result.txt)
```

For each tool, overlap uses the first occurence of a file that matches the pattern associated with that tool. If there are no files that match for the tool, then overlap assumes that the tool was not run.

_Further info: `kickoff_overlap.R` uses the R function `list.files` with the `pattern` parameter set to the pattern. The pattern is used as a regular expression used to match file names. The `soapFuse` pattern `.final.Fusion.specific.for.genes` starts with a period because the actual output file is named `[sample].final.Fusion.specific.for.genes`. Our pattern will match `soapFuse` output files with different samples names, like `s1.final.Fusion.specific.for.genes` or `random_sample_name.final.Fusion.specific.for.genes`._

### Example

Given the following directory structure

```txt
/data/fusion_output_1/
    final-list_candidate-fusion-genes.txt
    FusionDetection.FusionReport.Table.txt
    star-fusion.fusion_predictions.abridged.tsv
    DRAGEN.fusion_candidates.final
```

run the command

```bash
cd /data/fusion_output_1
~/igm-ensemble-fusion-detection/R/assemble_results.R --sample s1
```

or

```bash
cd ~/igm-ensemble-fusion-detection
R/assemble_results.R --sample s1 --baseDir /data/fusion_output_1
```

and overlap will be able to find fusion output files for `fusionCatcher`, `fusionMap`, `starFusion`, and `dragen`. Overlap will be able to find the correct files even if they are under subdirectories since the search is recursive.

```txt
/data/fusion_output_1/
    fusioncatcher/
        final-list_candidate-fusion-genes.txt
    fusionmap/
        results/
            FusionDetection.FusionReport.Table.txt
    starfusion/
        star-fusion.fusion_predictions.abridged.tsv
    dragen/
        DRAGEN.fusion_candidates.final
```

### `--starFusion` override example

You can override the default fusion output file search with direct fusion output file names.

```bash
R/assemble_results.R --samples s1 --starFusion path/to/custom-star-fusion-output.tsv
```

This command will recursively search under the current directory `.` for the `dragen`, `fusionCatcher`, `fusionMap`, etc., default files, but it will not perform the same search for the `starFusion` file. It will use the exact location to the `starFusion` file as provided by the user `path/to/custom-star-fusion-output.tsv`.

The overlap code only _knows_ how to parse the data from the output files in the default list. Make sure your custom locations point to files that hold the same data. There may be multiple output files from each fusion caller. Our overlap code is only parsing certain files for certain data to generate the overlap report.


### Analysis

The database uses the term 'analysis' to refer to a fusion caller execution on a single sample. A `STAR-Fusion` run on `sample1` is considered an analysis, say `analysis1`. A `JAFFA` run on `sample1` is considered a different analysis, say `analysis2`.

#### Analysis upload process

When uploading data to the database, the upload script uses the concept of an analysis to decide whether the data already exists in the database. If the analysis already exists in the database, then the upload script will skip the upload for that analysis. If there is already `STAR-Fusion` data for `sample1` in the database (`analysis1`), then any additional attempt to upload or replace `STAR-Fusion` data for `sample1` will fail.

- Upside - Since an analysis will not be uploaded twice, you can re-run the upload script on the same results directory and not have to worry about overwriting data. This can be useful if, for example, you run `STAR-Fusion` and `JAFFA` on `sample1` and upload the results. At a later time, you decide to run `FusionMap` on the same sample. When you want to add the `FusionMap` results to the database, you can run the upload script on the `sample1` results directory and it will only add the new `FusionMap` analysis to the database.
- Downside - You cannot replace analysis data using the upload script. If, for example, you run a new version of `STAR-Fusion` on `sample1`, you might want to replace the old `STAR-Fusion` results in the database. The upload script does not handle this scenario and you will have to perform a manual replacement.

### Upload with R instead of bash

Use `upload_fusion_results.R` for a greater level of control over the database upload (or if you don't have `bash` on your system to run `kickoff_upload.sh`).

```bash
R/upload_fusion_results.R -h
```

The result

```txt
Usage: R/upload_fusion_results.R [options]

This script acts as an importer for fusion detection results to the fusion detection database.

Options:
        --subject=SUBJECT
                Subject ID, e.g. subject1. (required)
        --sample=SAMPLE
                Sample ID, e.g. sample1. (required)
        --phenotype=PHENOTYPE
                Phenotype of the subject. (optional)
        --storage=STORAGE
                Storage method for the sample, such as 'Frozen' or 'FFPE'. (optional)
        --reads=READS
                Number of read pairs for the sample. (optional)
        --tool=TOOL
                Custom tool for which the results are being entered, possible values are 'dragen', 'FusionCatcher', 'FusionMap', 'Jaffa', 'MapSplice', 'SoapFuse', 'StarFusion', and 'TophatFusion'. (optional)
        --results=RESULTS
                Custom location of the results file for the associated tool, e.g. 'star-fusion.fusion_predictions.abridged.tsv'. (optional)
        -h, --help
                Show this help message and exit
```

`upload_fusion_results.R` must be run from the results directory. It searches for results files under `.`. Alternatively, you can upload results one-by-one by providing the `--tool` and `--results` parameters.
