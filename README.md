# Final-Project
##TRGN-510 Final Project - "Identifying key structure targets in a cohort of RNA-binding proteins"

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

We will focus on nine RBPs including DGCR8, ILF3, TARDBP, HNRNPU, PCPB2, ZNF622, PTBP1, EWSR1, and LARP7. For each RBP, analysis will be based on fastq files, a text file that contains the sequence data from clusters that pass filter on a flow cell (<https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html>, Illumina®). Fastq files for each RBP are studied in either HEPG2 (liver cancer cells) or K562 (leukemia cancer cells) 

_Section 3_ Proposed analysis for project:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Identifying structure targets based on the characteristics that describe the nature of RNA binding to RBPs.

_Section 4_ Proposed Timeline & major milestones:

_Section 5_ User interface:
