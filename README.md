# TRGN-510 Final Project - "Identifying key structure targets in a cohort of RNA-binding proteins"

_Section 1_ Description of final project:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As a PhD candidate under supervision of Dr. Zhipeng Lu, I focus on improving various crosslinking tools to study RNA structures and RNA-protein interactions which are understudied. There are classes of proteins referred to as “RNA-binding proteins (RBPs)” that bind single-stranded RNA (ssRNA) and double-stranded RNA (dsRNA) important in basic biology and genetic/autoimmune diseases. My final project involves identification of structure targets for a cohort of these RBPs that bind either ssRNA, dsRNA, or both. Currently, crosslinking immunoprecipitation (CLIP) techniques allow identification of RBPs but lacks crosslinking efficiency. By application of RBP database “enhanced crosslinking immunoprecipitation (eCLIP)” (Van Nostrand et al. 2016) and new processing pipeline parameters our lab has generated, we can improve understanding of structure targets for a subset of RBPs. 

_Section 2_ Datasets involved in final project:

__eCLIP__ 

is an RNA dataset that the general public of researchers may peruse via ENCODE database at <http://encodeproject.org>. There are filters that one can explore based on the research question being addressed such as assay type, status of research published on RNA/RBP, project (based on assay type), genome assembly, etc. For this project we will work with the following filters to analyze specific RBPS -

- Assay type - RNA binding
- Assay title - eCLIP
- Project - ENCODE
- Genome Assembly - GRCh38
- Target category - control, RNA binding protein
- Organism - _Homo sapiens_

An example of exploring an RBP in the ENCODE database page: 

![image](https://github.com/maiacorpuz/Final-Project/blob/master/00_eCLIP_RBP_ENCODE_ex.png)

![image](https://github.com/maiacorpuz/Final-Project/blob/master/00_eCLIP_RBP_ENCODE_ex.2.png)

We will focus on nine RBPs including DGCR8, ILF3, TARDBP, HNRNPU, PCPB2, ZNF622, PTBP1, EWSR1, and LARP7. For each RBP, analysis will be based on fastq files, a text file that contains the sequence data from clusters that pass filter on a flow cell (<https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html>, Illumina®). Fastq files for each RBP are studied in either HEPG2 (liver cancer cells) or K562 (leukemia cancer cells). I will select a total of four fastq sequences for each RBP, two being a pair for a control/mock input and another pair for the experiment condition. To obtain control sequences, navigate the ENCODE page under "__Summary__" lists "__Controls:__" which is a hyperlink to files for control fastq sequences.

<p align=center>
<img src="https://github.com/maiacorpuz/Final-Project/blob/master/00_mock_input_sequence.png" width="400" height="300">
</p>

From this new page, we should see options of mock input for the RBP of interest. 

![image](https://github.com/maiacorpuz/Final-Project/blob/master/00_mock_input_sequence_2.png)

We navigate toward the bottom of the page to see that there are a pair of fastq files hyperlinked in the __association graph__.

![image](https://github.com/maiacorpuz/Final-Project/blob/master/00_mock_input_sequence_3.png)

Download both fastq files to begin working on control sample. Similarly, within this page you may find 2-3 replicates that each have two downloadable fastq files for experiment condition. 

_Section 3_ Proposed analysis for project:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Identifying structure targets based on the characteristics describing RNA binding to RBPs include a processing of raw fastq sequences. The pre-processing of these sequences will include removal of adapter sequences and spliced reads. To do this step I will extract fastq sequences from eCLIP database and work in simple Linux environment. Following this, we can post-process the formatted sequences to generate a summary of gapped reads via Commandline and STAR 2.7.3a. STAR, a Commandline-run package, is a short read mapper that allows alignment to reference human genomes synthesized by researchers. STAR has built-in capacity to handle many customized parameters for alignment of these sequences with gapped and chimeric reads (Dobin et al. 2013). Also, examination of which RNA are top targets based on most "hits" via written custom scripts. Plot of gap length distribution (N) on intron and exon locations will begin as custom scripts in Linux and extracted to build RNA models in Integrated Genome Viewer (IGV). This serves to make qualitative and quantitative interpretations of the summaries gathered on gapped reads and locations.

_Section 4_ Proposed Timeline & major milestones:

The following timeline will be based on three major milestones that require outlining in this markdown. Overall, the project will be carried out over the course of one month. The completion date is set as a Part A: Beta Test of program and Part B: Submission of program to Dr. Craig.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;* ~~Milestone 1 (11/14/19)~~ - Pre-processing fastq sequences. Pre-processing fastq sequence data as mentioned above will likely take some time to understand and conduct. This is mainly due to understanding the function of demultiplexing sequences and how to remove adapter sequences using the CutAdapt tool used in eCLIP method.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _Known issues for Milestone 1_ - In completing Milestone 1, the major setbacks were understanding which programs to install for demultiplexing and CutAdapt steps of the eCLIP processing pipeline. In addition, it was difficult at first to determine whether fastq sequences downloaded from ENCODE database were already demultiplexed. This is mainly due to inexperience in understanding contents within a given fastq file. A benefit of learning the demultiplexing step was that I got a feel for Java and common workflow (CWL). These programming languages are generally used in computational analysis across many fields beside biomedical research.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;* ~~Milestone 2 (11/21/19)~~ - Summarization of gapped reads. Following removal of adapter sequences and spliced reads, it will be important to map these sequences to the genome assemblies "hg38" and "RfamhumanrnamRNA" in STAR genome mapper. Our lab has written custom scripts to perform these tasks. Summarization of gapped reads will be performed in Unix language that can be displayed in a table format.

__Summary__: I completed a first round of mapping alignments for LARP7 against RfamhumanrnaMrna and hg38 reference genomes. Below is a table of data summarized from aligning a mock input of LARP7 with hg38. 

<p align="center">
<img src="https://github.com/maiacorpuz/Final-Project/blob/master/2_STAR_out_LARP7mock_hg38.png" width="300" height="400">
</p>  
This summary output from STAR alignment holds many characteristics of mapped data but we can focus on the number of uniquely mapped reads and calculate the total mapped reads based on unique + multi-mapped reads.

![image](https://github.com/maiacorpuz/Final-Project/blob/master/02_Mapping_stats_LARP7.csv)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _Known issues for Milestone 2_ - Though mapping of aligned reads was executed, there is a lot of data to interpret which may be improved by adjusting STAR parameters to re-align. Ensuring quality of the performed mapping is essential. I will need to perform further processing of these mapped reads that are now in a file format known as 'sam'. This includes removing PCR duplicates and spliced reads via custom scripts. 
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;* ~~Milestone 3 (11/27/19)~~ - Examine top ranked RNA. ~~Once the pre-processed RNA has been aligned in STAR according to custom parameters, we will need to use SAMTools to write a script to process GTF and SAM files. In doing so, we will transfer top ranked examples to display in IGV.~~ I have never worked with STAR in Unix or IGV so I expect a lot of time spent understanding and applying knowledge into performing these steps.

__Summary__: I learned how to use slurm commands on USC's HPC (High Performance Computing) to conduct Cutadapt (removal of adapter sequences), STAR mapping to reference genomes, filter out spliced reads (python script), and use Samtools to sort and index the output sam files after STAR aligner. In addition, I produced bam files from converting the sam files (important as a file format readable in IGV) to view in Integrated Genome Viewer and begin analyzing the sequenced data. While troubleshooting the execution of splicefilter.py (python script to remove spliced reads written by Dr. Lu), I learned how to manage a virtual python environment in HPC server to allow completion.

<p align="center">
<img src="https://github.com/maiacorpuz/Final-Project/blob/master/3_IGV_bam_file_Rfam_chr22.png" width="400" height="400">
</p>

Having produced bam files from the mapped data against reference genomes hg38 and RfamhumanrnaMrna (collection of Rfam and RNA transcripts put together by Dr. Lu), I am able to use Samtools to get quick statistics of the mapped data. Running, 
` samtools idxstats 20191119_LARP7_mock_Aligned_hg38_nosplice_Sortedindex.bam `
I am able to identify the reference sequence name, sequence length, # of reads mapped, and # of reads unmapped from this Samtools option. I can then pipe the results into a tsv (comma-delimited) file by running
` samtools idxstats 20191119_LARP7_mock_Aligned_hg38_nosplice_Sortedindex.bam > idxstats_20191119_LARP7_mock_Aligned_hg38.tsv `
These tsv files, one for each reference genome and condition (mock/experiment), can be further sorted in an Excel sheet by listing highest reads at the top of the csv file. All four tsv files are available in this directory that serve as tables to show the __top-ranked RNA__ when mapped against RfamhumanrnaMrna and hg38. 

![image](https://github.com/maiacorpuz/Final-Project/blob/master/4.1.1_LARP7_mock_hg38.md)

I have uploaded slurm files to execute steps of the pipeline we are proposing such as Cutadapt and STAR mapping to reference genomes. They are written in bash scripts and can be run on USC's HPC server. Files with numbered prefix in this directory 'Final Project' include: [Cutadapt](https://github.com/maiacorpuz/Final-Project/blob/master/01_Cutadapt%20(removal%20of%20spliced%20reads)), [STAR aligner](https://github.com/maiacorpuz/Final-Project/blob/master/02_STAR%20aligner), and [Samtools](https://github.com/maiacorpuz/Final-Project/blob/master/03_Samtools). 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _Known issues for Milestone 3_ - Using slurm commands have a learning curve since it is based on trial and error. If running a job on an interactive node, there is no log output of the job executed unless you specify a new log file be created. Therefore, if a job is taking an extensive amount of time to run and you would like a log of the final output and/or errors, it's better to use a remote node that can be run in the background without you actively logged on to HPC's server. In addition, it takes some time to understand how to prepare your sam file output from STAR aligner using Samtools and finally inputting into IGV.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;* Milestone 4 (__after scope of course__) - Continuation and more summarization of sequence data. ~~The process of examining top ranked RNA based on scripts written in Unix may take more than one week so it may continue to be worked on in Milestone 3~~. Though I proposed to completely process eCLIP data for the nine RBP listed in the project description: DGCR8, ILF3, TARDBP, HNRNPU, PCPB2, ZNF622, PTBP1, EWSR1, and LARP7, only LARP7 has been completely processed through the pipeline. Running through each step of this processing pipeline to a complete extent took much more time and priority rather than analyzing all RBP data in each processing step one at a time. The next and last major task will be to plot gap length distribution (N) on intron and exon locations. This is also Unix-based so writing custom scripts that will allow exporting tables and figures will be the mode of output.

_Section 5_ User interface:

A wire diagram that illustrates the pipeline process to include eCLIP dataset, pre-processing, post-processing by our lab, and summary data of sequences to outline structure targets for the subset of RBPs listed above (Lu et al. 2016, Van Nostrand et al. 2016).

<p align="center">
<img src="https://github.com/maiacorpuz/Final-Project/blob/master/eCLIP_PARIS_RBP_processing_pipeline_v1.3.png" width="470" height="600">
</p>

` Citations: `

1. Davis, C. A., Hitz, B. C., Sloan, C. A., Chan, E. T., Davidson, J. M., Gabdank, I., ... & Onate, K. C. (2017). The Encyclopedia of DNA elements (ENCODE): data portal update. Nucleic acids research, 46(D1), D794-D801.

2. Van Nostrand, E. L., Pratt, G. A., Shishkin, A. A., Gelboin-Burkhart, C., Fang, M. Y., Sundararaman, B., ... & Stanton, R. (2016). Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP). Nature methods, 13(6), 508.

3. Lu, Z., Zhang, Q. C., Lee, B., Flynn, R. A., Smith, M. A., Robinson, J. T., ... & Mesirov, J. P. (2016). RNA duplex map in living cells reveals higher-order transcriptome structure. Cell, 165(5), 1267-1279.
