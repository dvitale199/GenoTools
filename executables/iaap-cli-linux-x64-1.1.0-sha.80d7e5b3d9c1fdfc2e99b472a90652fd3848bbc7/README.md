# IAAP-CLI

## Overview
IAAP-CLI is a command line interface application intended for cross-platform data analysis of Illumina microarray data.
Currently it only supports the genotyping (gencall) command which converts raw intensity files (IDATs) to multiple formats with genotypes (e.g. GTC).

## Prerequisites
There may be additional shared libraries that need to be installed on your [system]("https://docs.microsoft.com/en-us/dotnet/core/linux-prerequisites?tabs=netcore2x") related to the .NET core runtime for Linux.
For any installation problems see the **Troubleshooting** section or contact Illumina support.

## Installation
1. The un-tarred package should look something like this:
```
|-- iaap-cli
|   |-- lots of .dll files
|   |-- iaap-cli <-- executable
|-- licenses.txt
|-- README.md
|-- README.pdf
```

2. **(Optional)** Add the installation directory to your $PATH environmental variable. Example for bash:
```
echo "export PATH=$PATH:/path/to/iaap-cli-linux-x64-${VERSION}" >> ~/.bash_profile
source ~/.bash_profile
```

## Examples
- **View help for the `gencall` command**  
```
iaap-cli gencall -h
```  
- **Convert some GSA v1 raw intensity samples (IDAT files) to GTC files:**  
```
iaap-cli gencall /path/to/GSA-24v1-0_A1.bpm /path/to/GSA-24v1-0_A1_ClusterFile.egt /path/to/my_output -f /path/to/GSA_idats -g
```  
- **Convert some GSA v1 raw intensity samples (IDAT files) to GTC files using a Sample Sheet:**  
```
iaap-cli gencall /path/to/GSA-24v1-0_A1.bpm /path/to/GSA-24v1-0_A1_ClusterFile.egt /path/to/my_output -s /path/to/GSA_SampleSheet.csv -g
```  
- **Convert some GSA v1 raw intensity samples (IDAT files) to GTC files using 4 parallel threads:**  
```
iaap-cli gencall /path/to/GSA-24v1-0_A1.bpm /path/to/GSA-24v1-0_A1_ClusterFile.egt /path/to/my_output -f /path/to/GSA_idats -g -t 4
```
- **Convert some GSA v1 raw intensity samples (IDAT files) to PED files:**  
```
iaap-cli gencall /path/to/GSA-24v1-0_A1.bpm /path/to/GSA-24v1-0_A1_ClusterFile.egt /path/to/my_output -f /path/to/GSA_idats -p
```  
- **Convert some GSA v1 raw intensity samples (IDAT files) to GTC files using an include list of loci:**  
```
iaap-cli gencall /path/to/GSA-24v1-0_A1.bpm /path/to/GSA-24v1-0_A1_ClusterFile.egt /path/to/my_output -f /path/to/GSA_idats -g -inc /path/to/include_list.csv
```

### Sample Sheet

IAAP v1.1 adds the “sampleSheetPath” parameter as an alternative input to the “idatFolder” parameter. This parameter should be the full path location to a sample sheet file. The sample sheet should be a comma-separated (csv) file and is intended to be backwards-compatible with existing GenomeStudio-based sample sheets. However, only the [Data] section is used in IAAP and the [Header] and [Manifests] sections will be ignored.  An example sample sheet is shown below:
```
[Data]					
Sample_ID,SentrixBarcode_A,SentrixPosition_A,Path
SampleA,9348171001,R01C02,/full/path/to/idats/A/9348171001
SampleB,9348171002,R01C01,/full/path/to/idats/B/9348171002
```
Currently the `Sample_ID,SentrixBarcode_A,SentrixPosition_A,Path` columns are required.
Extra columns such as `Aux` or `Gender` will not result in a parsing failure.
If the sample sheet was generated in GenomeStudio, the `Path` must be mapped from the Windows paths
to their Unix paths.

### Include File

IAAP v1.1 adds the “includeFilePath” optional parameter to enable analysis on only the “Markers of Interest” from a product. It will allow filtering the SNPs/loci reported in output files (e.g GTC and PED) to only be those that appear in the “include file”. the file path should be the full path location to a csv file where each row lists a unique LociName that appears in the product’s manifest as shown below:
e.g.
```
rs1234
rs1235
rs1236
```
There are a couple of additional outcomes possible when passing an include file parameter compared to when doing a regular analysis.
- If the file includes loci that do not appear in the manifest, the analysis will continue but the application will output a file called badLoci.txt. The output file contains any loci not found. 
- The application will output an extra manifest (bpm) that is a subset/filtered version of the original manifest based on the loci in the include list (e.g. GSA-24v1-0_A1_filtered_read_only.bpm). This filtered manifest is what should be used during secondary analysis of the filtered GTC file (e.g with the [BeadArrayFiles] (https://github.com/Illumina/BeadArrayFiles) or GTCtoVCF libraries), instead of the full manifest.

**A MAJOR CAVEAT FOR THE “FILTERED” MANIFEST: This manifest is NOT intended to be re-used for primary analysis (ie, generation of GTC files, loading projects in GenomeStudio, etc.). The filtered manifest is to be used only for secondary analysis! Because of some details associated with the GenTrain algorithm, GTC files produced using the full manifest vs. the filtered manifest could result in discordance between base calls for the same samples.**

## Troubleshooting
    
### Missing icu library
```
$ ./iaap-cli
FailFast: Couldn't find a valid ICU package installed on the system. Set the configuration flag System.Globalization.Invariant to true if you want to run with no globalization support.

   at System.Environment.FailFast(System.String)
   at System.Globalization.GlobalizationMode.GetGlobalizationInvariantMode()
   at System.Globalization.GlobalizationMode..cctor()
   at System.Globalization.CultureData.CreateCultureWithInvariantData()
   at System.Globalization.CultureData.get_Invariant()
   at System.Globalization.CultureData.GetCultureData(System.String, Boolean)
   at System.Globalization.CultureInfo.InitializeFromName(System.String, Boolean)
   at System.Globalization.CultureInfo.Init()
   at System.Globalization.CultureInfo..cctor()
   at System.StringComparer..cctor()
   at System.AppDomainSetup.SetCompatibilitySwitches(System.Collections.Generic.IEnumerable`1<System.String>)
   at System.AppDomain.PrepareDataForSetup(System.String, System.AppDomainSetup, System.String[], System.String[])
Aborted (core dumped)
```
**Solution** 
Install `icu` examples:
- For Ubuntu run, `apt-get install icu-devtools`
- For CentOS/RHEL run, `yum install libicu`

Or, to avoid installation of `icu`, you can disable [globalization support]("https://github.com/dotnet/corefx/blob/master/Documentation/architecture/globalization-invariant-mode.md#enabling-the-invariant-mode") for .NET core via setting the evironmental variable `DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=true`.

### Out of memory exception
```
Killed.
```

**Solution**
This error is the result of the "OOM killer" of the linux kernel killing iaap-cli because it requested too much memory from the system.
This can occur if `--num-threads` is set to high for your system or the system does not have enough memory to analyze the product.
Reduce the `--num-threads` or remove the option to default to 1 thread.
