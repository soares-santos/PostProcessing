README for code that creates lists of regions for viewing in ds9

This code creates a separate list file for each CCD containing regions for loading into ds9 .fits image viewing software.

If you have questions, please feel free to email me (Bobby Butler) at bbutler3@uchicago.edu. I would be happy to help, and am usually reasonably quick at responding.

*****
INPUT AND OUTPUT

This code takes the same input config.txt as the 'sb3.py' code. Note that the outdir will be the same, and it will create quite a few files (as many as there are usable CCDs). 

The output files will be named similarly to the following example: 

ds9regions_dp107_475914_i_01.reg

In this example, the season is 107, the exposure number is 475914, the band is i, and the CCD number is 01.

*****
COLOR CODING AND LABELING

There will in general be three different region types, with a fourth (*) if a cand_list is specified in the config.txt file (see sbREADME.txt for a bit more information on this):

1) RED BOXES: These are rejects from the filterObj.out file for each CCD. They will always be labeled with -9.
2) MAGENTA CIRCLES: These are accepted reals from the filterObj file. They will be labeled with their 'SNOBJID' from the filterObj file.
3) CYAN CIRCLES: These are accepted fakes from the filterObj file. They will be labeled with their 'SNFAKE_ID' from the filterObj file.
*4) GREEN CIRCLES: These are accepted reals that have passed an additional cut that was imposed in the creation of the cand_list file, if you've specified one. They will be labeled with their 'SNID' from that file. If you didn't specify a cand_list file (and put 'none' instead, as is convention here), these will just end up magenta with SNOBJID like all the other accepted reals.
