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

![image](https://github.com/maiacorpuz/Final-Project/blob/master/eCLIP%20RBP%20ENCODE%20ex.png)

![image](https://github.com/maiacorpuz/Final-Project/blob/master/eCLIP%20RBP%20ENCODE%20ex.2.png)

We will focus on nine RBPs including DGCR8, ILF3, TARDBP, HNRNPU, PCPB2, ZNF622, PTBP1, EWSR1, and LARP7. For each RBP, analysis will be based on fastq files, a text file that contains the sequence data from clusters that pass filter on a flow cell (<https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html>, Illumina®). Fastq files for each RBP are studied in either HEPG2 (liver cancer cells) or K562 (leukemia cancer cells). I will select a total of four fastq sequences for each RBP, two being a pair for a control/mock input and another pair for the experiment condition. To obtain control sequences, navigate the ENCODE page under "__Summary__" lists "__Controls:__" which is a hyperlink to files for control fastq sequences.

![image](https://github.com/maiacorpuz/Final-Project/blob/master/mock_input_sequence.png)

From this new page, we should see options of mock input for the RBP of interest. 

![image](https://github.com/maiacorpuz/Final-Project/blob/master/mock_input_sequence_2.png)

We navigate toward the bottom of the page to see that there are a pair of fastq files hyperlinked in the __association graph__.

![image](https://github.com/maiacorpuz/Final-Project/blob/master/mock_input_sequence_3.png)

Download both fastq files to begin working on control sample. Similarly, within this page you may find 2-3 replicates that each have two downloadable fastq files for experiment condition. 

_Section 3_ Proposed analysis for project:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Identifying structure targets based on the characteristics describing RNA binding to RBPs include a processing of raw fastq sequences. The pre-processing of these sequences will include removal of adapter sequences and spliced reads. To do this step I will extract fastq sequences from eCLIP database and work in simple Linux environment. Following this, we can post-process the formatted sequences to generate a summary of gapped reads via Commandline and STAR 2.7.3a. STAR, a Commandline-run package, is a short read mapper that allows alignment to reference human genomes synthesized by researchers. STAR has built-in capacity to handle many customized parameters for alignment of these sequences to handle gapped and chimeric reads (Dobin et al. 2013). Also, examination of which RNA are top targets based on most "hits" via written custom scripts. Plot of gap length distribution (N) on intron and exon locations will begin as custom scripts in Linux and extracted to build RNA models in Integrated Genome Viewer (IGV). This serves to make qualitative and quantitative interpretations of the summaries gathered on gapped reads and locations.

_Section 4_ Proposed Timeline & major milestones:

The following timeline will be based on three major milestones that require outlining in this markdown. Overall, the project will be carried out over the course of one month. The completion date is set as a Part A: Beta Test of program and Part B: Submission of program to Dr. Craig.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;* Milestone 1 - Pre-processing and summarization of gapped reads. Pre-processing sequence data as mentioned above will likely take some time to understand and conduct. This is mainly due to replicating work produced by Dr. Gene Yeo and Dr. Eric Van Nostrand labs. Following removal of adapter sequences and spliced reads, it will be important to map these sequences to the genome assemblies "hg38" and "RfamhumanrnamRNA" in STAR genome mapper. Our lab has written custom scripts to perform these tasks. Summarization of gapped reads will be performed in Unix language that can be displayed in a table format.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;* Milestone 2 - Examine top ranked RNA. Once the pre-processed RNA has been aligned in STAR according to custom parameters, we will need to use SAMTools to write a script to process GTF and SAM files. In doing so, we will transfer top ranked examples to display in IGV. I have never worked with STAR in Unix or IGV so I expect a lot of time spent understanding and applying knowledge into performing these steps.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;* Milestone 3 - Continuation and more summarization of sequence data. The process of examining top ranked RNA based on scripts written in Unix may take more than one week so it may continue to be worked on in Milestone 3. The last major task will be to plot gap length distribution (N) on intron and exon locations. This is also Unix-based so writing custom scripts that will allow exporting tables and figures will be the mode of output.


_Section 5_ User interface:

A wire diagram that illustrates the pipeline process to include eCLIP dataset, pre-processing, post-processing by our lab, and summary data of sequences to outline structure targets for the subset of RBPs listed above.

![image](https://github.com/maiacorpuz/Final-Project/blob/master/eCLIP_PARIS_RBP_process_pipeline.png =200x300)
