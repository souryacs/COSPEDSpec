*********************************
COSPEDSpec
*********************************

COSPEDSpec is a python based tool for computing species tree from a set of incongruent gene trees:
1. First a couplet based supertree method COSPEDTree (please check the reference at the end) is applied as its first stage, to generate a possibly 
non-binary species tree S' covering all the input taxa.
2. Subsequently we apply a refinement on S', with reducing the number of extra lineages (deep coalescence) with 
respect to the input gene tree set, for producing a strict binary species tree.

Input source trees can be either in NEWICK format or in NEXUS format. 
However, all the source trees should have identical input formats. They should be placed in a standard tree list file, 
according to the syntax of NEXUS or NEWICK formats. Such a tree list text file is to be provided as an input of this executable.

Source trees may or may not be weighted; however, species tree generation does not use any branch length information.

Output supertree is generated in the NEWICK format. It is unweighted.

*********************************
Dependencies
*********************************

COSPEDSpec is developed in Linux Systems (Ubuntu 14.04), using Python 2.7.

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default) 

Note: We have not tested the code on Python 3. Any user having Python 3 environment need to check the correct 
execution of our code, and optionally needs to upgrade it accordingly.

We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ ) 

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. We did not upgrade the 
code for Dendropy 4.0 support, so any user having this new version of Dendropy might need to check the functionalities of 
COSPEDSpec and possibly upgrade / replace / edit few dendropy related functions. So, we recommend users to 
use the earlier version of Dendropy, to avoid any conflict.

Support for Dendropy 4 and corresponding update of code will be done in a future release.

3) Numpy ( available on the link: http://www.numpy.org/ )

User can install Numpy using pip (python software downloader tool) module, which contains the latest Numpy module in it. 
We found that Numpy module in the traditional Apt-Get repository is of a previous version. So we recommend users to 
use pip downloader, instead.

***************
Command line options
****************

"-I", "--INPFILE"

          specifies name of the input file containing gene trees
			  
"-O", "--OUTFILE"

          specifies name of the output file to contain target species tree
  
"-p", "--inpform"
			  
        specifies the file format (input gene tree list containing file)
        its value can be 1 (corresponds to the newick format - default) or 2 (corresponding to the nexus format)

"-q", "--queues", \

        specifies the number of priority queues which will contain the support scores of individual couplets, 
        during execution of the supertree algorithm.
        
        Its value can be one of the following:
        a) Value of 1 indicates that only a single max priority queue is used for storing the support score values. 
        b) Value of 2 (default option) indicates that two separate queues are used to store the support scores. First queue is used to 
        store for non-conflicting couplets (all input trees having same relation for the couplet), whereas the second queue is 
        used to store the support scores of conflicting couplets (here at least two different relations of that couplet are supported by 
        the input trees).
        
        
An example of a command line option: 

./COSPEDSpec -I source_treelist.txt -p1 -O out_species_tree.txt

where, source_treelist.txt contains gene tree list,
and out_species_tree.txt will contain the derived species tree.

    
*********************
Citation
********************
If you use COSPEDSpec, please cite:

Sourya Bhattacharyya, Jayanta Mukhopadhyay, Couplet Supertree based Species Tree Estimation, ISBRA 2015, LNBI 9096, pp. 48–59, 2015.


Additionally, for details of our proposed couplet based supertree approach, you may check the following paper:

Sourya Bhattacharyya, Jayanta Mukherjee, COSPEDTree: COuplet Supertree by Equivalence Partitioning of taxa set and DAG formation, 
IEEE/ACM Transactions on Computational Biology and Bioinformatics. volume 12, number 3, pages 590-603, 2015.


*********************************
For any queries, please contact
*********************************

Sourya Bhattacharyya 
Department of Computer Science and Engineering
Indian Institute of Technology Kharagpur
<sourya.bhatta@gmail.com>



