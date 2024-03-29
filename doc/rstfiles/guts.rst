The guts of sfs_coder
*********************

sfs_coder has an object oriented structure.  The main classes are Mutation,
Simulation, and SFSCommand.  The Mutation class stores detauls relevant to
each variable site in the simulation (including sites that have fixed in the
output).  The Simulation class stores details relevant to every mutation 
within the simulation or the simulation as a whole (such as the number 
of simulated loci and the distances between them) and provides a number of 
methods for analyzing data.  The SFSCommand class stores details that pertain 
to command lines and allows for building commands and analyzing data.

In the next few subsections, I discuss a few relevant details of the 
organization and use of each class, but please refer to the section
``sfs_coder Classes and methods'' section of this documentation for further 
details.

Mutation class
==============

Each Mutation instance has several attributes which correspond to the output
from SFS_CODE.  Briefly, SFS_CODE outputs a comma delimited string similar to
the following for each mutation within a simulation:

0,A,483,5574,5676,GGG,T,1,G,C,-4.613899e-03,5,2.1080,2.1850,2.2951,2.3299,2.4194;

Each field is explained in detail in the SFS_CODE documentation 
(http://sfscode.sourceforge.net). Suppose we have a Mutation object named mut. 
The fields correspond to:

1) locus number (0-based) -> mut.locus
2) Autosomal/X/Y -> mut.AXY
3) position in locus (0-based) -> mut.pos
4) Generation mutation occurred -> mut.t_int
5) Generation mutation fixed (time of sampling if not fixed) -> mut.t_fix
6) Ancestral tri-nucleotide -> mut.tri_nuc
7) Derived nucleotide -> mut.deriv_N 
8) Synonymouys (0) or non-synonymous (1) (0 is used for non-coding) -> mut.non_or_syn
9) Ancestral amino acid -> mut.ancest
10) Derived amino acid -> mut.deriv_aa
11) Fitness effect (s) -> mut.fit
12) Number of chromosomes that carry mutation -> mut.pops_numchr 
13) comma delimited list of chromosomes carrying mutation (-1 indicates all chromosomes)
    -> mut.chrs

Most of the attributes are pretty straightforward, but a few are stored as dictionaries.
The reason we've done this is because some data points can vary between populations 
(e.g., the time that a mutation fixes) and others are population invariant (e.g., position).

So if you want to access a fixation time, list of chromosomes carrying the mutation,
or the number of chromosomes carrying the mutation within a given population, all of these
attributes are accessed as mut.<attribute_name>[pop] where pop is the number of the 
population of interest.  For example, mut.chrs[0] is a dictionary of chromosomes that carry
the mutation in the 0th population.

Simulation class
================

The Simulation class keeps track of a few essential things, such as the command line
that was used to generate the data and the Mutation objects in the simulation.

Notably, Mutation objects are stored both as a list (self.muts) and a dictionary (self.loci).
The reason for this is that sometimes we don't want to go through every Mutation when 
computing something (e.g., the diversity within a particular locus).

The Simulation class provides methods for computing pi, ZnS, S, Tajima's D, theta_H, 
and Watterson's theta. Other statistics will be forthcoming. It is possible to compute these 
statistics either for every site in the simulation, for a specified list of loci, or 
for a specified set of genomic coordinates in the case of simulations that use the 
"genomic" method.

The Simulation class also provides the ability to simulate phenotypes under some canonical
models.  Details on this functionality will be forthcoming shortly.

SFSCommand class
================

This class stores and parses command lines.  The command line is necessary in order to
analyze the output of an SFS_CODE simulation becuse we need to know the number of simulated
sites, the number of sampled chromosomes, the number of populations, etc.


Troubleshooting
===============

1) The long form swtiches for SFS_CODE are NOT currently supported.  Thus, if you run 
a simulation with the option '-t 0.001' you're good to go but '--theta 0.001' could
cause you problems.  If you use sfs_coder to run your simulations you (hopefully!) should
not have any problems, but if you ran SFS_CODE on your own and used these options it may not 
be feasible to use sfs_coder to analyze the output.  

2) sfs_coder was designed with the hope that it will be easy to use.  It is not designed to be
fast.  Sorry.  I may try to make performance improvements over time but no promises.

3) If you're using the 'genomic' method to build your simulations, sfs_coder creates a bunch
of files that detail the start and end points of the simulation, the recombination map for 
the region of interest, and all the simulated genomic elements.  If you move these files
posthoc, sfs_coder will not be able to find them because it gets their file paths from
the sfs_code command line, which is the first line in the output file.  So don't move those
files.  

I admit this could be a problem if you like to do your simulations on one machine and your
analysis on another machine. If enough people complain about this I will add functionality
to allow the user to specify the locations of the files on the local machine.

Department of Complaints
========================

Feel free to email Lawrence at lawrence.uricchio@ucsf.edu about problems. We hope this 
software will be useful, so if we can facilitate your research then we are happy to help.


