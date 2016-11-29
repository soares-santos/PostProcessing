README for SB code

This code computes the surface brightness of the background behind gravitational wave source candidates, and performs some other miscellaneous tasks. It outputs one .fits table per exposure it processes. The column names of this file are listed below in the OUTPUT section.

*****
this is undoubtedly unclear at some points and probably not comprehensive; if you have questions, please feel free to email me (Bobby Butler) at bbutler3@uchicago.edu. I would be happy to help, and am usually reasonably quick at responding.
*****

*small disclaimer: We haven't quite worked out the best way to compute the error in SB, so the SB_err column in the output .fits table should essentially be ignored. This is mentioned again far below in column descriptions.

-----
config.txt

1) The file should be formatted in the same way as the file /data/des41.b/data/rbutler/sb/example_config.txt is. The only space in each line should be between the colon and the start of the entry to the right (the space is the delimiter that the code uses to parse each line).

2) 'outdir' should ideally be an empty directory, or at least one where you don’t have anything valuable (The code deletes some files as it goes; it should only delete ones that it has created itself, but if you happen to have files that match or resemble the names in the same directory, they could get confused. Small chance, but fair warning.) 

3) MASTER FAKE/CANDIDATE LISTS: If you don't have a master fake list and/or a master list of 'good' candidates, you can just put 'none' where you would normally put 'path/filename'; the code can work around that. Specifying a master FAKE list will allow true_mag and true_flux for fakes to be included in the output file. Specifying a master CANDIDATE list will allow candidates which are detected in both the i and z band AND that have a machine learning score above a certain threshold (it was 0.5 for the one I used). 

If you're unsure of the format these files usually take, you can take a look at these examples (they're in the format that will work in this code):
FAKES: /data/des41.b/data/rbutler/10_19_2016/fake_list_dp107_sorted_uniq.tab
CANDIDATES: /data/des41.a/data/marcelle/G184098/analyze/107/output/candidates_table_real.txt

4) ccds: This will generally be 'all', but in special cases can be changed. In these cases, you can specify an array of candidates that is formatted in the same way as the expnums array in the first line of the file.

----- 
RUNNING THE CODE

0) If the normal conda environment from DES is already in use, the joblib package is available, and you've already cloned the 'sb3.py' code itself, you ***don't need the next step (1)***

1) The code can run using just the normal conda environment setup for DES, along with an extra package called joblib (for parallelizing). You can set all that up in your current directory by running

source /data/des41.b/data/rbutler/sbsetup

This will also put your very own copy of the SB code in your current directory (called 'sb3.py' right now; if that changes, I'll adjust this and the sbsetup file).

2) If you have the config.txt file (it can actually be called whatever you want) set up, you're ready to run the code. To do this, it's just 

python sb3.py config.txt

The code will spit out a few warnings at the start about non-integer indices, but these are nothing to worry about. If you don't want to see them, you can run

python -W ignore sb3.py config.txt

instead.

3) As it runs, the code will be printing quite a few things to the screen. You don't really have to pay attention to any of it. It's mostly there to help me see what's going on while I'm working on the code, but I find that it's good just to see that the code is running. If you've specified a master fake list, the code will print out quite a bit more things.

4) SPEED: The speed of the code can vary quite a bit, depending on several things: the number of input exposures, how many detections those exposures have, whether you specify a master fake list to match with (if yes, it's quite a bit slower), and how many cores you tell the code to use. As you probably know, the DES machines have 32 cores; I'd recommend using around 16-20 (the default is 16), to find the happy medium between speed and hogging all the cores. If you open the main code, the variable 'cpus' is at the top just under the package import section. That's what you'd change.

You can expect the code to take anywhere from 2 to 12 minutes per exposure, depending on how 'full' of detections they are. This is if you're matching fakes. If you're not, it's quite a bit faster--I've experienced ~2 minutes maximum per exposure, but I haven't tried it many times without a fake list.

-----
OUTPUT

1) While many smaller files will be created and destroyed along the way, the final output will be as follows: one log file for the entire run, and one .fits table per exposure processed. 

LOG FILE: The log file will be stored in the 'outdir' you specified in the config.txt file. The naming convention is 'log_year-month-day_hrmin.txt', where the time information corresponds to when the code began to run. So, for example, if you started the code at 11:26 PM on November 16, 2016, the name would be 'log_2016-11-16_2326.txt'. This is just to ensure no log files overwrite each other.
Contents: For each exposure, the log file will have exposure number, band, # of CCD processed, total unmasked area for each exposure (in square degrees), # of non-rejected real detections, and # of non-rejected real detections per square degree. These last two numbers become matched "good" candidates if you have specified a 'cand_list'.

.FITS TABLE: The fits table for each exposure will also be stored in the 'outdir' specified in the config.txt file. The naming convention is 'final_expnum.fits'. So for example, for the output of exposure 475914, it would just be 'final_475914.fits'.
***NOTE: If there is already a file in the 'outdir' with exactly the same name as one that will be produced by the code, the code WILL STOP and give you an error. Make sure that's not the case before you run the code, or you might be unhappy 20 minutes in when you have to start over.***

2) .FITS TABLE cont'd: List of columns in .fits output file:
[ID, candidate_ID, RA, DEC, radius, band, mag, mag_err, flux, flux_err, true_mag, true_flux, NITE, EXPNUM, CCDNUM, SEASON, SB, SB_err, fakeID]

There are several conventions used within the data of the fits tables. Some are the same as used in other files up the pipeline, and some I created myself. Below is a list of columns that you should know about.

ID: For rejected detections (already rejected in the filterObj file), the detection will have ID=-9. Otherwise it will be some 7-digit (or similar) valid number.

candidate_ID: If you didn't specify a 'cand_list' file, this is always -99. If you did, it's some positive integer ID (which will be *different* than the normal ID) ONLY for matched "good" candidates.

radius: This refers to the radius in arcsec of the circle drawn around the target RA and DEC in the template image used to compute the surface brightness. Right now, it's just set at a default of 2. When this number is 0, the detection is a reject.

true_mag and true_flux: As mentioned far above, these are only populated by a number other than 0 if the detection is a non-rejected fake. If you've specified a 'fake_list', non-rejected fakes will have a positive value for each. If you haven't, each defaults to -999.

SB: This is 0 when a background surface brightness is not computed (for rejects). When one is computed, it should be some reasonable positive number. For backgrounds with average pixel value <= 0, SB is fixed at 30.

SB_err: As of now, this number should be disregarded. We have not determined a reliable way of computing this value.

fakeID: If the detection is a fake, this is some 6-digit (or similar) number. If it's real, fakeID is 0.



