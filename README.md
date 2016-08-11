# PostProcessing
Post process steps to be run using GWforce and GWMakeDatafiles
This code requires the following codes to run:
  diffimg.py 
  postproc.py 
  postproc.ini 
  diffimg_setup.sh
  FAKES_OVERLAID_KNFakes57279noHost.DAT
  HTML.py

First you must source diffimg_setup.sh in order to set up the environment necessary to run the postproc.py script

You can edit the postproc.ini file in order to change the default names of certain files that are to be produced and read, as well as the default output directory.

Whenever a change is made to the GWMakeDataFIles code that involves the output data, make sure you check diffimg.py to see if that code must be updated to read new output data.

TO TEST USE SEASON 106 and EXPNUMS 475914 475915 475916 482859 482860 482861

In order to run this code you must give the input arguments: Season (just the number) as well as the list of exposure numbers you wish to run over.

Example Syntax: python postproc.py --season 501 --expnums 445544 333444

You may also give the output directory that you wish for all plots to be made in as well as where the preliminary webpages are made.

The code is currently set up to handle any number/type of filter however the lightcureve section of the coding requires a title change depending on wht the filters being used are, as the current title only expects i and z. The code will still run without such a change but you won't be able to tell what color corresponds to what filter without the title as there is currently no legend provided (too dificult to automate for a variable number of filters).

The lightcurves are all placed in a subdirectory inside whatever output directory is chosen.

Efficiency plots are also produced and put in the output directory
