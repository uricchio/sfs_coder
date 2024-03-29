Running SFS_CODE simulations
****************************

The basic structure of an sfs_coder script
==========================================

Running SFS_CODE simulations with sfs_coder requires only a few lines
of code.   

First, we import the command module, initialize an SFSCommand object, and tell
the software where the sfs_code binary is located.

.. code-block:: python

   from sfscoder import command
  
   com = command.SFSCommand()

   com.sfs_code_loc = '/path/to/sfs_code'

Next, we need to build a command line.  Although this process is
flexible (in fact we can build any command line that is accepted by SFS_CODE),
we have prepackaged several models that may be of general interest.  For
example, to simulate the model of Gutenkunst (2009, *PLoS Genetics*), we
call the following:

.. code-block:: python

   com.three_pop(model='gutenkunst')

   com.execute()   

And that's it!  Of course, there are many more options that can be altered to
modify the parameters of the simulation, such as the ancestral population size. 
Please see the "scripts" directory in the top level of sfs_coder for more 
examples. Below, we include a slight modification of the above 
script that demonstrates some basic functionality that may be useful to users,
as well as a few other examples.

.. code-block:: python

   from sfscoder import command
   import os
   from random import randint

   # initialize an SFS_CODE command, set the prefix of the output subdirectory
   com = command.SFSCommand(prefix='guten.N500')

   # build the command line for the Gutenkunst model, specifying some parameters
   com.three_pop(N=500,nsam=50,nsim=10,model='gutenkunst')

   # set the location of the sfs_code binary
   com.sfs_code_loc = os.path.join(os.path.expanduser('~'),
                          'path/to/sfs_code')

   # execute the command, supplying a random number
   com.execute(rand=randint(1,100000))

Examples
========

Adding selection
^^^^^^^^^^^^^^^^
.. code-block:: python
   
   from sfscoder import command
   import os
   from random import randint

   # initialize a new SFS_CODE command
   com = command.SFSCommand(prefix='tennessen.N1000')

   # build the command line for the tennessen model
   # a selection model is added with sel = sel=['-W','1','5','0','1']
   # this adds a type 1 selection model, with gamma =5, 
   # and the probability of negative seleciton set to 1.
   # for more on selection models in SFS_CODE, see the SFS_CODE handbook

   com.three_pop(N=1000,nsam=[50,50,0],nsim=10,model='tennessen',
                 L=['-L','1','100'],sel=['-W','1','5','0','1'])

   # set the location of sfs_code
   com.sfs_code_loc = os.path.join(os.path.expanduser('~'),
                       'path/to/sfs_code')

   # execute the command 
   com.execute(rand=randint(1,100000))

Simulations with realistic genomic structure and demography
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from sfscoder import command
   import os
   from random import randint

   # initialize a new SFS_CODE command
   com = command.SFSCommand(prefix='guten.lactase')

   # build the command line for the gutenkunst model in the lactase region
   com.genomic(N=100,model='gutenkunst',sel=False)

   # set the location of sfs_code
   com.sfs_code_loc = os.path.join(os.path.expanduser('~'),'path/to/sfs_code')

   # execute the command 
   com.execute(rand=randint(1,100000))


Using SGE
^^^^^^^^^

sfs_coder uses the sge_task_id system variable to number output files.  
If you submit an sfs_coder script to a cluster as an array job, it will take
care of all the work of numbering the output files for you.

For example, any of the above scripts can be sent to a cluster with the 
following header:

.. code-block:: none

   #!/usr/bin/python
   #$ -e sim.div.log
   #$ -o sim.div.log
   #$ -S /usr/bin/python
   #$ -cwd
   #$ -r yes
   #$ -l h_rt=240:00:00
   #$ -t 1-100
   #$ -l arch=linux-x64
   #$ -l mem_free=1G
   #$ -l netapp=1G

Simulations of phenotypes
^^^^^^^^^^^^^^^^^^^^^^^^^

Simulation of phenotypes is handled post-hoc to the simulation
of genotypes using either the genotypes or selection coefficients 
to pick effect sizes.  We provide three models, those of Wu (2011, 
AJHG), Eyre-Walker (2010, PNAS), and Simons (2014, Nature Genetics).

Briefly, the model of Wu takes the effect size of a variant to be 
proportional to log10 of the allele frequency in the sample.  Only 
5% of the variants under 3% frequency are taken as causal. The model of
Eyre-Walker takes effect size to be equal to the selection coefficient 
multiplied by (1 + e), where e is a normally distributed random variable.
For the model of Simons, we take the effect size to be proportional to 
the selection coefficient with probability rho, but randomly sample the 
effect size from the distribution of selection coefficients with probability
(1-rho).  There are several parameters that can be altered under each model.
Please see the documentation for the sim_pheno method for more information
on all the parameters that can be chosen.

For each method, the user can set the proportion of the variance in the 
phenotype that is explained by the sequence in question with the h_sq 
parameter.

.. code-block:: python

   #!/usr/bin/python
   #
   #

   import sys
   from sfscoder import sfs

   # an sfs_code output file that we will analyze
   f =sys.argv[1]

   # initialize a data object and set the file path
   data = sfs.SFSData(file=f)

   # get all the data from the simulations in the file

   data.get_sims()

   # simulate phenotypes using one of a few different models
   # Here we are using the 'EW' method, the model of Eyre-Walker
   # (2010, PNAS)

   for sim in data.sims:
       p = sim.sim_pheno(method='EW',pops=[0,1])

