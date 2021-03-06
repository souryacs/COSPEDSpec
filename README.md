*********************************
COSPEDSpec
*********************************

COSPEDSpec is a python based tool for computing species tree from a set of incongruent gene trees:

Here, species tree is generated by analyzing couplet based relations obtained from input gene trees. 
Such relations are syntesized to first produce  a possibly non-binary species tree S' covering 
all the input taxa. Subsequently we apply a refinement on S', which produces a strict binary species tree. The 
refinement uses a couplet based measure termed as 'excess gene leaf count'.

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

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. 
We did not upgrade the code for Dendropy 4.0 support. So, any user having 
this new version of Dendropy might need to check the functionalities of 
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

-I INP_FILENAME, --INPFILE=INP_FILENAME

                name of the input file containing gene trees

-O OUT_FILENAME, --OUTFILE=OUT_FILENAME

                name of the output file to contain target species tree

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

                1 - input file format is NEWICK (default)
                2 - input file format is NEXUS

-d DIST_MAT_TYPE --distmat=DIST_MAT_TYPE

                1 - Distance matrix entries for refinement of species tree are computed by taking average 
                of the excess gene leaf count values for individual couplets.
                2 - Distance matrix entries are computed by using the mean of (average and mode based 
                filtered average) of the excess gene leaf count values for individual couplets. (Default)

-r TAXON_NAME, --ROOT=TAXON_NAME

		User can specify a taxon name to root the output species tree with that specified taxon.
        
An example of a command line option: 

./COSPEDSpec -I source_treelist.txt -p1 -O out_species_tree.txt

where, source_treelist.txt contains gene tree list,
and out_species_tree.txt will contain the derived species tree.

    
*********************
Citation
********************
If you use COSPEDSpec, please cite:

Sourya Bhattacharyya, Jayanta Mukhopadhyay, Couplet Supertree based Species Tree Estimation, 
ISBRA 2015, LNBI 9096, pp. 48–59, 2015.

*********************************
For any queries, please contact
*********************************

Sourya Bhattacharyya 
Department of Computer Science and Engineering
Indian Institute of Technology Kharagpur
<sourya.bhatta@gmail.com>



