# PostProcessing
Post process steps to be run using GWforce and GWMakeDatafiles
This code requires the following codes to run:
  diffimg.py 
  postproc.py 
  postproc.ini 
  diffimg_setup.sh
  FAKES_OVERLAID_KNFakes57279noHost.DAT
  HTML.py


TO TEST USE SEASON 106 and EXPNUMS 475914 475915 475916 482859 482860 482861

In order to run this code you must give the input arguments: Season (just the number) as well as the list of exposure numbers you wish to run over.

Syntax: python postproc.py --season 501 --expnums 445544 333444

You may also give the output directory that you wish for all plots to be made in as well as where the preliminary webpages are made.

The code is currently set up to handle any number/type of filter however the lightcureve section of the coding requires a title change depending on wht the filters being used are, as the current title only expects i and z. The code will still run without such a change but you won't be able to tell what color corresponds to what filter without the title as there is currently no legend provided (too dificult to automate for a variable number of filters).

The lightcurves are all placed in a subdirectory inside whatever output directory is chosen.

Efficiency plots are also produced and put in the output directory
