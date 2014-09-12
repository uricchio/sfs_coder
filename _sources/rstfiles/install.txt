Downloading and Installing
**************************

A link to the sfs_coder tarball will appear here at release.

Installation
============

Installation is handled by python's distutils.

Using a terminal, execute the following commands:

1) cd into the sfs_coder source directory.

   .. code-block:: none
     
      cd /path/to/sfs_coder

2) run setup.py install

   .. code-block:: none
   
      python setup.py install
      
If you don't get any errors with the above steps,
you should be able to open python and import sfscoder.

    .. code-block:: none
      
      python
  
    .. code-block:: python

      import sfscoder

If you still don't have any errors, you're good to go.

Installing if you don't have admin status on your machine
---------------------------------------------------------

It's possible you will get a "permission denied" error after
step 2 if you don't have write access to the default install 
location on your machine.  In that case, try

     .. code-block:: none

        python setup.py install --user

which should install in a location that you have access to.
Repeat the last two steps above to make sure it's working
(i.e., try to import sfscoder).  If it's still not working
feel free to contact us.  We've only tested this on a few 
Unix based platforms and it's possible (likely?) that 
some system specific problems could pop up. 

Required dependencies
=====================

* `SFS_CODE`_
  
  .. _SFS_CODE: http://sfscode.sourceforge.net

SFS_CODE can be installed anywhere on the user's machine. The path
to the binary is supplied to sfs_coder (see the section on 
running SFS_CODE through sfs_coder).

* `python`_ (*2.7 or greater*)

  .. _python: https://www.python.org

* `numpy`_ 

  .. _numpy: http://www.numpy.org/

Optional dependencies
=====================

* `mpmath`_ (*required for rescaled recurrent hitchhiking simulations*)

  .. _mpmath: https://code.google.com/p/mpmath/

* `scipy`_ (*required for rescaled recurrent hitchhiking simulations*)

  .. _scipy: http://www.scipy.org/


Importing sfs_coder's modules             
=============================

All of the tools to execute simulations and analyze command
lines are in the command module.

.. code-block:: python

   from sfscoder import command

   com = command.SFSCommand()


All of the tools to analyze the output of SFS_CODE simulations
are in the sfs module.

.. code-block:: python

   from sfscoder import sfs

   mut = sfs.Mutation()

