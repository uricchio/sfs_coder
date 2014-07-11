Installation
************

sfs_coder does not require any installation of its own in order 
to get its basic functionality.  Just download the source and alter the 
PYTHONPATH variable on your machine so that python knows where to find the 
modules (see the section on "Importing" below). 

However, some methods within sfs_coder require external software in order to 
run.

Required dependencies
=====================

* `SFS_CODE`_
  
  .. _SFS_CODE: http://sfscode.sourceforge.net

SFS_CODE can be installed anywhere on the user's machine. The path
to the binary is supplied to sfs_coder (see the section on 
running SFS_CODE through sfs_coder).

* `python`_ (*2.7 or greater*)

  .. _python: https://www.python.org

Optional dependencies
=====================

* `mpmath`_ (*required for rescaled recurrent hitchhiking simulations*)

  .. _mpmath: https://code.google.com/p/mpmath/

* `scipy`_ (*required for rescaled recurrent hitchhiking simulations*)

  .. _scipy: http://www.scipy.org/

* `matplotlib`_ (*required for plotting the output*)

  .. _matplotlib: http://matplotlib.org/


Importing sfs_coder's modules             
=============================

Python uses the PYTHONPATH system variable to search for modules that are 
imported.  Suppose we download sfs_coder and store it in the directory 
'~/sfs_coder', and then we try to execute the following script called 
'basic.py':

.. code-block:: python

   import command

   com = SFSCommand()

If the directory that contains command.py ('~/sfs_coder/src' by default) is not 
included in the PYTHONPATH variable, this will result in an error similar to 
the following:

.. code-block:: python

   Traceback (most recent call last):
     File "basic.py", line 13, in <module>
       import command
   ImportError: No module named command

Python does not know where the command module is!  To fix this, we can add the
'~/sfs_coder/src' directory to the PYTHONPATH variable in a couple different 
ways.

Adding the sfs_coder source directory to PYTHONPATH in .bashrc
--------------------------------------------------------------

If you execute your scripts at the command line with a bash shell, you can add
a line to your .bashrc file that will fix this problem and allow you to run the
above script.  

.. code-block:: bash

   export PYTHONPATH=$PYTHONPATH:~/sfs_coder/src

The .bashrc file exists in your home directory and is read by bash every time 
you open a new shell.  If a file called .bashrc doesn't exist in your home 
directory you can create it.

Of course, if the path to your sfs_coder 'src' directory is different than
above you will need to provide the path to your copy of this directory. Note
that you will need to relaunch bash (open a new shell) once you have altered
the .bashrc.  

Note, this solution does not seem to work well with SGE (Sun Grid Engine) 
clusters, I use the following solution for scripts that I run on a
cluster.

Adding the path to sfs_coder's source directory within a python script
----------------------------------------------------------------------

You can also add the path to sfs_coder's 'src' directory to any python script
if for any reason you don't want to modify your .bashrc as above.  Assuming
the same directory layout as the above example, we can use:

.. code-block:: python

   import sys
   sys.path.append('~/sfs_coder/src')
   import command

   com = SFSCommand()

This adds '~/sfs_coder/src' to the PYTHONPATH variable within the script.
Note that this solution requires us to add this line of code to every python
script that imports something from sfs_coder, whereas the first solution allows 
us to import the modules just like any other python modules.

