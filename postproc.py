import shutil
import numpy as np
import matplotlib.pyplot as plt 
import pyfits as py
import argparse
import ConfigParser
import glob, sys, datetime, getopt
import os
import subprocess
import diffimg

#Read User input#

print "Read user input"

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--expnums', metavar='e',type=int, nargs='+', help='List of Exposures', default= [476960])

parser.add_argument('--outputdir', metavar='d', type=str, help='Directory location of output files', default= "testevent")
#parser.add_argument('--season', help='season is required', default=107, type=int)

args = parser.parse_args()
print args.expnums
print args.outputdir

outdir = str(args.outputdir)

if not os.path.exists(outdir):
    os.mkdir(outdir) 

if not os.path.exists(outdir + '/' + 'stamps'):
    os.mkdir(outdir + '/' + 'stamps')

if not os.path.exists(outdir + '/' + 'plots'):
    os.mkdir(outdir + '/' + 'plots')

outplots = outdir + '/' + 'plots'
outstamps = outdir + '/' + 'stamps'

print "Read environment variables"

forcedir = os.environ["TOPDIR_SNFORCEPHOTO_IMAGES"]

season= os.environ["SEASON"]
run = "dp"+str(season)

#print season

#Read Config File#

print "Read config file"

config = ConfigParser.ConfigParser()
inifile = config.read('./postproc.ini')[0]

expdir = config.get('data', 'exp')
ncore = config.get('GWFORCE', 'ncore')
numepochs_min = config.get('GWFORCE', 'numepochs_min')
writeDB = config.get('GWFORCE', 'writeDB')

format= config.get('GWmakeDataFiles', 'format')
numepochs_min = config.get('GWmakeDataFiles', 'numepochs_min')
trigger = config.get('GWmakeDataFiles', '2nite_trigger')
outFile_stdoutreal = config.get('GWmakeDataFiles-real', 'outFile_stdout')
outDir_datareal = config.get('GWmakeDataFiles-real', 'outDir_data')
outFile_stdoutfake = config.get('GWmakeDataFiles-fake', 'outFile_stdout')
outDir_datafake = config.get('GWmakeDataFiles-fake', 'outDir_data')
fakeversion = config.get('GWmakeDataFiles-fake', 'version')


print "Check RUNMON outputs"

#expnumlist = args.expnums

#Read in and locate files#
for expnum in args.expnums: 
    e=str(expnum)
    print "Check dir content for exposure " +e
    d= expdir+"/*/"+e+"/"+run
    runmonlog=d+"RUNMON*.LOG"
    nfiles= len(glob.glob(runmonlog))
    if nfiles != 1:
        print "runmonlog not found"
###Think about what to do in the case that the file is not found###
        psf= forcedir+"/*/*"+e+"*.psf"
        diff = forcedir+"/*/*"+e+"*_diff.fits"
        diffmh = forcedir+"/*/*"+e+"*_diff_mh.fits"

        for filetype in (psf, diff, diffmh):
            if len(glob.glob(filetype)) == 0 :
                print "files not found"

print "Run GWFORCE"

#run "gwforce" section#
#forcePhoto_master.pl     \
#   -season         107   \
#   -numepochs_min  0     \
#   -ncore          4     \
#   -writeDB 

a= 'forcePhoto_master.pl ' + ' -season ' +str(season) + ' -numepochs_min ' +numepochs_min + ' -ncore ' +ncore

if writeDB == "on":
    a = a+ ' -writeDB ' 

print a
#subprocess.call(a, shell=True)


####run "gwhostmatch" section (if time allows)#


print "Run GWmakeDataFiles - real"

#run "Gwmakedatafiles" section#

#makeDataFiles_fromSNforce \
#   -format snana \
#   -season 201     \
#   -numepochs_min 0 \
#   -2nite_trigger iz \
#   -outFile_stdout  makeDataFiles_real.stdout  \
#   -outDir_data   GWevent2_numepoch1_iz_real_text \

b= 'makeDataFiles_fromSNforce' + ' -format ' +format + ' -numepochs_min ' +numepochs_min + ' -2nite_trigger ' +trigger + ' -outFile_stdout ' +outFile_stdoutreal + ' -outDir_data ' +outDir_datareal

print b 
#subprocess.call(b, shell=True)

print "Run GWmakeDataFiles - fake"

b= 'makeDataFiles_fromSNforce' + ' -format ' +format + ' -numepochs_min ' +numepochs_min + ' -2nite_trigger ' +trigger + ' -outFile_stdout ' +outFile_stdoutfake + ' -outDir_data ' +outDir_datafake + ' -fakeVersion ' +fakeversion

print b
#subprocess.call(b, shell=True)



print "Plot efficiency"

#Make plots Section#
###Plot1 Efficiency Plot ###
###Plot5 Magerror Distribution ###
###Plots should include all bands###

reals = diffimg.DataSet(outDir_datareal, label = 'reals')
fakes = diffimg.DataSet(outDir_datafake, label = 'fakes')
###Need to generate fakes input on own###
fakes.get_fakes_input(config.get('GWmakeDataFiles-fake', 'fake_input'))
truth = fakes.fakes_input

rdatag = reals.set_mask(PHOTFLAG_bit=4096)
fdatag = fakes.set_mask(PHOTFLAG_bit=4096)

bins = bins = np.arange(17,25,0.5)

###Generalize code to handle all bands/any combo of bands###

fmaski = (fdatag.BAND=='i')
fmaskz = (fdatag.BAND=='z')
tmaski = (truth.BAND == 'i')
tmaskz = (truth.BAND == 'z')

f_ihist, bin_edges = np.histogram(fdatag.SIMMAG[fmaski],bins=bins)
f_zhist, bin_edges = np.histogram(fdatag.SIMMAG[fmaskz],bins=bins)

t_ihist, bin_edges = np.histogram(truth.TRUEMAG[tmaski], bins=bins)
t_zhist, bin_edges = np.histogram(truth.TRUEMAG[tmaskz], bins=bins)

plt.figure()
plt.plot(bins[1:], f_ihist*100.0/t_ihist, label= 'i-band', lw=4, color='orange')
plt.plot(bins[1:], f_zhist*100.0/t_zhist,label='z-band',lw=4,color='darkblue')
plt.scatter(bins[1:], f_ihist*100.0/t_ihist,lw=4,color='orange')
plt.scatter(bins[1:], f_zhist*100.0/t_zhist,lw=4,color='darkblue')
plt.title('Efficiency: Blue = z  Orange = i')
plt.xlabel('Magnitude')
plt.ylabel('Percent Found')
plt.savefig(outplots + '/'+'efficiency.pdf')
plt.clf()

print "Plot DeltaMag/MAGERR histogram"

#Histogram of DeltaMag/MAGERR#

deltai = fdatag.MAG[fmaski] - fdatag.SIMMAG[fmaski]
deltaz = fdatag.MAG[fmaskz] - fdatag.SIMMAG[fmaskz]
deltaiovererr = deltai/(fdatag.MAGERR[fmaski])
deltazovererr = deltaz/(fdatag.MAGERR[fmaskz])
bins2 = np.arange(-30,30,.1)
iweights = np.ones_like(deltaiovererr)/float(len(deltaiovererr))
zweights = np.ones_like(deltazovererr)/float(len(deltazovererr))

deltaiovererr_hist, bin_edges= np.histogram(deltaiovererr, weights= iweights, bins=bins2)
deltazovererr_hist, bin_edges= np.histogram(deltazovererr, weights= zweights, bins=bins2)

plt.figure()
plt.plot(bins2[1:], deltaiovererr_hist, label= 'i-band', lw=3, color='orange')
plt.plot(bins2[1:], deltazovererr_hist, label= 'z-band', lw=3, color='blue')
plt.title('Delta Mag over Mag Error')
plt.ylabel('Percent of Total')
plt.savefig(outplots +'/'+'DeltaoverERR.pdf')
plt.clf()

print "Number of candidates per ccd"
#Plot Candidates per CCD #

x = np.zeros(len(np.unique(reals.data.SNID)))
y = np.unique(reals.data.SNID)
for i in np.arange(len(x)):
    x[i] =reals.data.CCDNUM[reals.data.SNID==y[i]][1]

plt.hist(x, bins= np.arange(min(reals.data.CCDNUM), max(reals.data.CCDNUM) +2 ,1), color='orange')
plt.title('Hist of real candidates per CCD')
plt.ylabel('Number of Real Candidates')
plt.xlabel('CCD Number')
plt.savefig(outplots +'/'+'Hist_of_real_candidates_per_CCD.pdf')
plt.clf()

x = np.zeros(len(np.unique(fakes.data.SNID)))
y = np.unique(fakes.data.SNID)
for i in np.arange(len(x)):
    x[i] =fakes.data.CCDNUM[fakes.data.SNID==y[i]][1]


plt.hist(x, bins= np.arange(min(reals.data.CCDNUM), max(reals.data.CCDNUM) +2,1),color = 'orange')
plt.title('Hist of fake candidates per CCD')
plt.ylabel('Number of Fake Candidates')
plt.xlabel('CCD Number')
plt.savefig(outplots +'/'+'Hist_of_fake_candidates_per_CCD.pdf')
plt.clf()

print "Save candidates info"
###Write data files for each candidate including info discussed###



rID= reals.data.SNID
urID= np.unique(rID)
numofcan = len(urID)
realss = reals.data

f= open(str(outdir)+'/'+'allcandidates.txt', 'w')
f.write('SNID, ' + ' RA, ' + ' DEC, ' + ' CandType' +  ' NumEpochs, ' + ' NumEpochsml, ' + ' LatestNiteml' + "\n")
for i in range(0,numofcan):
    line = str(reals.data.SNID[i]) + ", " + str(reals.data.RA[i]) + ", " + str(reals.data.DEC[i]) + ", " + str(reals.data.CandType[i]) + ", " + str(reals.data.NumEpochs[i]) + ", " + str(reals.data.NumEpochsml[i]) + ", " + str(reals.data.LatestNiteml[i]) + "\n"
    f.write(line)

f.close()
#Fit light curves section (If time allows)#

for i in range(0,numofcan):
    Cand =(reals.data.SNID == urID[i]) 
    filename = 'Candidate_'+str(int(urID[i])) + '.html'
    header = 'BAND ' + 'x ' + 'y ' + 'Mag ' + 'Nite ' + 'MJD ' + 'Season' 
    nobs = len(reals.data.BAND[Cand])
    seasoncol = np.ones((nobs,), dtype = np.int)*int(season)
    print seasoncol
    table = np.column_stack((reals.data.BAND[Cand], reals.data.XPIX[Cand], reals.data.YPIX[Cand], reals.data.MAG[Cand], reals.data.NITE[Cand], reals.data.MJD[Cand], seasoncol))
    np.savetxt(str(outdir)+'/'+filename, table, fmt = '%s', header = header)
    ###Collect Stamps for observations of this candidate###
    thiscand_stampsdir = outstamps + '/' + str(urID[i])
    if not os.path.exists(thiscand_stampsdir):
        os.mkdir(thiscand_stampsdir)

    for j in range(0,nobs):
#        thisobs_nite = str(int(reals.data.NITE[Cand][j]))
#        thisobs_band = reals.data.BAND[Cand][j]
#These lines are temporary and must be replaced with general versions###
        thisobs_nite = "20150917"
        thisobs_band = "z"
        ccdnum = 25
        season = 44
        expdir = "/data/des41.a/data/marcelle/diffimg/local-runs"
        thisobs_ID = 74
###End of Temporary Lines###
        a = expdir + '/' + thisobs_nite + '/*/dp' + str(season) + '/' + thisobs_band + '_' + str(ccdnum) + '/stamps*'
        print a
        thisobs_stampsdir = glob.glob(a)[0]
        print thisobs_stampsdir
        filenamediff = thisobs_stampsdir + '/diff' + str(thisobs_ID) + '.gif'
        filenamesrch = thisobs_stampsdir + '/srch' + str(thisobs_ID) + '.gif'
#        filenametemp = thisobs_stampsdir + '/temp' + str(thisobs_ID) + '.gif' 
        shutil.copy(filenamediff, thiscand_stampsdir)
        shutil.copy(filenamesrch, thiscand_stampsdir)
#        shutil.copy(filenametemp, thiscand_stampsdir)
                     
                       
    
print "SUCCESS"    





#Save output files in specific and organized directories#

###Setup a display of the search, template and difference images for each candidate. Automatically save for top n candidates. Give option to user to pick candidate by ID to display###
