## Log for BB2490 Project

### February 17th 2016
We settled on GitHub both for hosting code and logbooks. Since the rest of the group were not very familiar with Git, we had a little workshop this afternoon where we set up the GitHub repo. We also decided on Python as our go-to scripting language and set up the IDE PyCharm on Amy's and Tobias's computers. We also made sure that UppMax has Git installed. 

### February 19th 2016
Meeting with Olof after seminar. 
To reduce memory footprint and computational time we will restrict to 3 individuals, rather than 8.
First step will be to remove adapters, check quality and other preprocessing. 
We discussied how we would do subsampling. Current idea is that sampling will be done after processing but before mapping.
Each sample was sequenced on 8 different lanes, four from each flow cell.

### February 23d 2016
I worked on preprocessing with Sailendra this afternoon. We started by becoming familiar [with the data](https://github.com/aerugo/BB2490-RNASeq-Project/wiki/Overview-of-data). It is structured as follows:

```
|-- Flowcell 1: 130104_SN866_0197_AC1DLVACXX
|   |-- SN_10_LPS: Sample from person 1, treated with LPS
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 1, adapter with tag CGATGT.
|   |-- SN_10_UNST: Sample from person 1, untreated
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 1, adapter with tag ATCACG.
|   |-- SN_11_LPS: Sample from person 2, treated with LPS
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 1, adapter with tag TGACCA.
|   |-- SN_11_UNST: Sample from person 2, untreated
|   |   `-- Paired end reads, fastq format from four lanes of flowcell 1, adapter with tag TTAGGC.
|   |-- SN_12_LPS: Sample from person 3, treated with LPS
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 1, adapter with tag GCCAAT.
|   |-- SN_12_UNST: Sample from person 3, untreated
|   |   `-- Paired end reads, fastq format from four lanes of flowcell 1, adapter with tag ACAGTG.
`-- Flowcell 2: 130104_SN866_0198_BC1DAYACXX
|   |-- SN_10_LPS: Sample from person 1, treated with LPS
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 2, adapter with tag CGATGT.
|   |-- SN_10_UNST: Sample from person 1, untreated
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 2, adapter with tag ATCACG.
|   |-- SN_11_LPS: Sample from person 2, treated with LPS
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 2, adapter with tag TGACCA.
|   |-- SN_11_UNST: Sample from person 2, untreated
|   |   `-- Paired end reads, fastq format from four lanes of flowcell 2, adapter with tag TTAGGC.
|   |-- SN_12_LPS: Sample from person 3, treated with LPS
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 2, adapter with tag GCCAAT.
|   |-- SN_12_UNST: Sample from person 3, untreated
|   |   `-- Paired end reads, fastq format from four lanes of flowcell 2, adapter with tag ACAGTG.
```
