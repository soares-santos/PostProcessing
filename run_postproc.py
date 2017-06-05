import os
import argparse
import ConfigParser
import sys
import postproc
import numpy as np
import pandas as pd

## Read config file
config = ConfigParser.ConfigParser()
if os.path.isfile('./postproc.ini'):
    inifile = config.read('./postproc.ini')[0]

## Read command line options
parser = argparse.ArgumentParser(description=__doc__, 
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--expnums', metavar='e',type=int, nargs='+', help='List of Exposures')
parser.add_argument('--outputdir', metavar='d', type=str, help='Location of output files')
parser.add_argument('--season', help='Season number', type=int)
parser.add_argument('--triggerid', help='LIGO trigger ID', type=str)
parser.add_argument('--mjdtrigger', type=float, help='MJD of LIGO trigger')
parser.add_argument('--ups', type=bool, help='ups mode: True/False')
args = parser.parse_args()

## Set ups mode: True/False
if args.ups == None:
    ups = config.getboolean('general','ups')
else:
    ups = args.ups

## If running in ups environment, replace the .ini file
if ups:
    cpath = os.environ["GWPOST_DIR"]
    inifile = config.read(os.path.join(cpath,"postproc.ini"))[0]

## Set outdir
if args.outputdir == None:
    outdir = config.get('general','outdir')
else:
    outdir = args.outdir

## Set season
if args.season == None:
    season = config.get('general','season')
else:
    season = str(args.season)

## Set triggerid
if args.triggerid == None:
    triggerid = config.get('general','triggerid')
else:
    triggerid = args.triggerid

## Set triggermjd
if args.mjdtrigger == None:
    triggermjd = config.getfloat('general','triggermjd')
else:
    triggermjd = float(args.mjdtrigger)

## Get the list of exposures
if args.expnums == None:
    expnums = []
    indir = config.get('general','indir')
    expnums_listfile = config.get('general','exposures_listfile')
    expnums_listfile = os.path.join(indir,expnums_listfile)
    try:
        explist = open(expnums_listfile,'r')
        expnums1 = explist.readlines()
        for line in expnums1:
            expnums.append(line.split('\n')[0])
            expnums = map(int,expnums)
    except:
        print "ERROR: List of exposures file not found or empty."
        expnums = []
        #sys.exit(1)
else:
    expnums = args.expnums

## Get remaining configuration from the .ini file 

rootdir = config.get('general','rootdir')
indir = config.get('general','indir')
expdir = os.path.join(rootdir,'exp')
forcedir = os.path.join(rootdir,'forcephoto')+'/images/dp'+season+'/*'

setupfile = config.get('general','env_setup_file')

triggerid = config.get('general','triggerid')
propid = config.get('general','propid')

mlscore_cut = config.getfloat('plots','mlscore_cut')

blacklist_file = config.get('masterlist', 'blacklist')
masterfile_1 = config.get('masterlist', 'filename_1')
masterfile_2 = config.get('masterlist', 'filename_2')

logfile = config.get('checkoutputs', 'logfile')
ccdfile = config.get('checkoutputs', 'ccdfile')
goodchecked = config.get('checkoutputs', 'goodfile')
steplist = config.get('checkoutputs', 'steplist')

ncore = config.get('GWFORCE', 'ncore')

numepochs_min_1 = str(config.getint('GWFORCE', 'numepochs_min'))
numepochs_min_2 = str(config.getint('GWmakeDataFiles', 'numepochs_min'))
#numepochs_min = str(min(numepochs_min_1,numepochs_min_2))

writeDB = config.getboolean('GWFORCE', 'writeDB')

version_hostmatch = config.get('HOSTMATCH', 'version')

db = config.get('general', 'db')
schema = config.get('general', 'schema')

filename = config.get('truthtable', 'filename')
truthplusfile = config.get('truthtable', 'plusname')

format = config.get('GWmakeDataFiles', 'format')

two_nite_trigger = config.get('GWmakeDataFiles', '2nite_trigger')

outFile_stdoutreal = config.get('GWmakeDataFiles-real', 'outFile_stdout')
outDir_datareal = config.get('GWmakeDataFiles-real', 'outDir_data')

combined_fits = config.get('GWmakeDataFiles-real', 'combined_fits')

outFile_stdoutfake = config.get('GWmakeDataFiles-fake', 'outFile_stdout')
outDir_datafake = config.get('GWmakeDataFiles-fake', 'outDir_data')

fakeversion = config.get('GWmakeDataFiles-fake', 'version')

## Make directory structure

if not os.path.isdir(outdir):
    os.mkdir(outdir) 

if not os.path.isdir(outdir + '/' + 'stamps'):
    os.mkdir(outdir + '/' + 'stamps')

if not os.path.isdir(outdir + '/' + 'plots'):
    os.mkdir(outdir + '/' + 'plots')

if not os.path.isdir(outdir + '/plots/' + 'lightcurves'):
    os.mkdir(outdir + '/plots/' + 'lightcurves')

outplots = outdir + '/' + 'plots'
outstamps = outdir + '/' + 'stamps'

#########
# STEP -1: Set up the environment
#########

print "Run STEP -1: Set up the environment"
postproc.prep_environ(rootdir,indir,outdir,season,setupfile,version_hostmatch,db,schema)
print

#########
# STEP 0: Create initial master list, check processing outputs
#########

print "Run STEP 0: Create initial master list, check processing outputs"
if len(expnums)>0:
    expniteband_df,master = postproc.masterlist(masterfile_1,blacklist_file,triggerid,propid,expnums)
else:
    print "No exposures specified by user. All exposures taken under trigger id "+str(triggerid)+" and prop id "+str(propid)+" will be used for the initial master list and the checkoutputs step."
    expniteband_df,master = postproc.masterlist(masterfile_1,blacklist_file,triggerid,propid)
print

### this method assumes .FAIL files are cleared out when a CCD is reprocessed
if len(expniteband_df)>0:
    expnums,a_blacklist,ccddf = postproc.checkoutputs(expniteband_df,logfile,ccdfile,goodchecked,steplist)
else:
    print "ERROR: No exposures provided, and no exposures found matching the trigger id and prop id provided in the .ini file:"
    print
    print "TRIGGER ID: "+str(triggerid)
    print "PROP ID: "+str(propid)
    print
    print "EXITING."
    print
    sys.exit()
print

#########
# STEP 1: Create final master list
#########

print "Run STEP 1: Create final master list"
expniteband_df,master = postproc.masterlist(masterfile_2,blacklist_file,triggerid,propid,expnums,a_blacklist)
print

expnums = expniteband_df['expnum'].tolist()
#sys.exit()

#########
# STEP 2: Forcephoto
#########

print "Run STEP 2: Forcephoto"
postproc.forcephoto(ncore,numepochs_min_1,writeDB)
print

#########
# STEP 3: Hostmatch
#########

print "Run STEP 3: Hostmatch"
#import desHostMatch
#desHostMatch.main()
print

#########
# STEP 4: Make truth table
#########

if len(expnums)>0:
    print "Run STEP 4: Make truth table"
    truthplus = postproc.truthtable(expnums,filename,truthplusfile)
else:
    print "WARNING: List of exposures is empty. Skipping STEP 4."
print

#########
# STEP 5: Make datafiles
#########

print "Run STEP 5: Make datafiles"
#postproc.makedatafiles(format,numepochs_min_2,two_nite_trigger,outFile_stdoutreal,outDir_datareal)

if not fakeversion=='KBOMAG20ALLSKY':
#    postproc.makedatafiles(format,numepochs_min_2,two_nite_trigger,outFile_stdoutfake,outDir_datafake,fakeversion)
    one=1
else:
    print "No datafiles made for fakes because fakeversion=KBOMAG20ALLSKY."
print                                                                                    

print "Run STEP 5b: Combine real datafiles"
fitsname = postproc.combinedatafiles(master,combined_fits,outDir_datareal)
print

#########
# STEP 6: Make plots
#########

skip=False
print "Run STEP 6: Make plots"
print
realdf = postproc.makeplots(ccddf,master,truthplus,fitsname,expnums,triggermjd,mlscore_cut,skip)
print

#########
# STEP 7: Make webpage
#########

print "Run STEP 7: Make webpage"
print "This is not yet implemented. Coming soon..."
postproc.createhtml(fitsname,realdf,master)
print






# ####run "gwhostmatch" section (if time allows)#
# if ups:
#     gwpostdir = os.environ['GWPOST_DIR']
#     deshostmatch = os.path.join(gwpostdir,'desHostMatch_v2.py')
# else:
#     deshostmatch = 'desHostMatch_v2.py'

# print '-'*50
# print 'Running',deshostmatch
# ##not this
# #print os.popen('python '+deshostmatch+' --season='+season).read()
# ###this
# #print os.popen('python '+deshostmatch+' '+season+' --username marcelle --password mar70chips --dbname destest --verbose --testdb').read()
# print 'Finished',deshostmatch
# print '-'*50

# #sys.exit()

# print "Run GWmakeDataFiles - real"

# #run "Gwmakedatafiles" section#

# #makeDataFiles_fromSNforce \
# #   -format snana \
# #   -season 201     \
# #   -numepochs_min 0 \
# #   -2nite_trigger iz \
# #   -outFile_stdout  makeDataFiles_real.stdout  \
# #   -outDir_data   GWevent2_numepoch1_iz_real_text \

# if not trigger=='null':
#     b= 'makeDataFiles_fromSNforce' + ' -format ' +format + ' -season '+ season  + '  -numepochs_min ' +numepochs_min + ' -2nite_trigger ' +trigger + ' -outFile_stdout ' +outFile_stdoutreal + ' -outDir_data ' +outDir_datareal
# else:
#     b= 'makeDataFiles_fromSNforce' + ' -format ' +format + ' -season '+ season  + '  -numepochs_min ' +numepochs_min + ' -outFile_stdout ' +outFile_stdoutreal + ' -outDir_data ' +outDir_datareal

# print 'B', b 
# #if running from ups we need to go to outdir because idk
# if ups:
#     gwpostdir = os.environ['GWPOST_DIR']
#     os.chdir(outdir)
#     os.system("cp "+gwpostdir+"/FAKES_OVERLAID_" + fakeversion + ".DAT " + outdir)
# #subprocess.call(b, shell=True)

# print "real MakeDataFiles complete"
# #sys.exit()
# #Run bobby's code here

# print "Run GWmakeDataFiles - fake"

# if not trigger=='null':
#     b= 'makeDataFiles_fromSNforce' + ' -format ' +format + ' -season ' + season + ' -numepochs_min ' +numepochs_min + ' -2nite_trigger ' +trigger + ' -outFile_stdout ' +outFile_stdoutfake + ' -outDir_data ' +outDir_datafake + ' -fakeVersion ' +fakeversion
# else:
#     b= 'makeDataFiles_fromSNforce' + ' -format ' +format + ' -season ' + season + ' -numepochs_min ' +numepochs_min + ' -outFile_stdout ' +outFile_stdoutfake + ' -outDir_data ' +outDir_datafake + ' -fakeVersion ' +fakeversion

# print b

# subprocess.call(b, shell=True)
# sys.exit()

# #if running from ups go back to ups dir
# if ups:
#     os.chdir(gwpostdir)

# #Produce Truth Table for Fakes#
# explist=','.join(map(str,oexpnums))

# # the database where diffimg outputs are stored                                 
# db='destest'
# schema = 'marcelle'

# # the query you want to run to get the truth table data                         
# query='select distinct SNFAKE_ID, EXPNUM, CCDNUM, TRUEMAG, TRUEFLUXCNT, FLUXCNT, BAND, NITE, MJD from '+ schema +'.SNFAKEIMG where EXPNUM IN ('+explist+') order by SNFAKE_ID'
# print query

# # the file where you want to save the truth table                              

 
# filename= config.get('GWmakeDataFiles-fake', 'fake_truth')
# filename=os.path.join(outdir,filename)
# connection=easyaccess.connect(db)
# connection.query_and_save(query,filename)
# connection.close()

# #sys.exit()
# ### FOR THE FIRST RUN EXIT HERE TO LEARN NAMES OF VALUES WE NEED###
# print "Data Made"
 

# print "Read Data"

# #Make plots Section#
# ###Plot1 Efficiency Plot ###
# ###Plot5 Magerror Distribution ###
# ###Plots should include all bands###

# #print os.path.join(outdir,outDir_datareal)
# #reals = diffimg.DataSet(os.path.join(outdir,outDir_datareal), label = 'reals')
# #fakes = diffimg.DataSet(os.path.join(outdir,outDir_datafake), label = 'fakes')
# reals = diffimg.DataSet(outDir_datareal, label = 'reals')
# fakes = diffimg.DataSet(outDir_datafake, label = 'fakes')
# ###Need to generate fakes input on own###
# os.system('mv fakes_truth.tab '+outdir)
# #fakes.get_fakes_input(os.path.join(outdir,config.get('GWmakeDataFiles-fake', 'fake_input')))
# #truth = fakes.fakes_input

# print '-----'
# #print reals
# print '-----'

# rdatag = reals.set_mask(PHOTFLAG_bit=4096)
# fdatag = fakes.set_mask(PHOTFLAG_bit=4096)
# colors = ['r','g','b','c','m','k','y']
# rID= reals.data.SNID
# urID= np.unique(rID)
# numofcan = len(urID)
# realss = reals.data
# bands = realss.BAND
# ubands = np.unique(bands)



# #bins = bins = np.arange(17,25,0.5)

# ###Generalize code to handle all bands/any combo of bands###

# #fmaski = (fdatag.BAND=='i')
# #fmaskz = (fdatag.BAND=='z')
# #tmaski = (truth.BAND == 'i')
# #tmaskz = (truth.BAND == 'z')
# #
# #print "Plot Efficiency"
# #
# #for i in range(0,len(ubands)):
# #    fmask= (fdatag.BAND == ubands[i])
# #    tmask = (truth.BAND == ubands[i])
# #    fhist, bin_edges = np.histogram(fdatag.SIMMAG[fmask],bins = bins)
# #    thist, bin_edges = np.histogram(truth.TRUEMAG[tmask], bins=bins)
# #    plt.plot(bins[1:], fhist*100.0/thist, label = str(ubands[i]), lw=4)
# #    plt.scatter(bins[1:], fhist*100.0/thist, lw=4)
# #    plt.title('Efficiency')
# #    plt.xlabel('Mag')
# #    plt.ylabel('Percent Found')
# #    plt.savefig(outplots + 'Efficiency for ' + str(ubands[i]) + '.png')
# #    plt.clf()
# #
# #
# #f_ihist, bin_edges = np.histogram(fdatag.SIMMAG[fmaski],bins=bins)
# #f_zhist, bin_edges = np.histogram(fdatag.SIMMAG[fmaskz],bins=bins)
# #
# #t_ihist, bin_edges = np.histogram(truth.TRUEMAG[tmaski], bins=bins)
# #t_zhist, bin_edges = np.histogram(truth.TRUEMAG[tmaskz], bins=bins)
# #
# #plt.figure()
# #plt.plot(bins[1:], f_ihist*100.0/t_ihist, label= 'i-band', lw=4, color='orange')
# #plt.plot(bins[1:], f_zhist*100.0/t_zhist,label='z-band',lw=4,color='darkblue')
# #plt.scatter(bins[1:], f_ihist*100.0/t_ihist,lw=4,color='orange')
# #plt.scatter(bins[1:], f_zhist*100.0/t_zhist,lw=4,color='darkblue')
# #plt.title('Efficiency: Blue = z  Orange = i')
# #plt.xlabel('Magnitude')
# #plt.ylabel('Percent Found')
# #plt.savefig(outplots + '/'+'efficiency.pdf')
# #plt.clf()
# #
# #print "Plot DeltaMag/MAGERR histogram"

# #Histogram of DeltaMag/MAGERR#

# #deltai = fdatag.MAG[fmaski] - fdatag.SIMMAG[fmaski]
# #deltaz = fdatag.MAG[fmaskz] - fdatag.SIMMAG[fmaskz]
# #deltaiovererr = deltai/(fdatag.MAGERR[fmaski])
# #deltazovererr = deltaz/(fdatag.MAGERR[fmaskz])
# #bins2 = np.arange(-30,30,.1)
# #iweights = np.ones_like(deltaiovererr)/float(len(deltaiovererr))
# #zweights = np.ones_like(deltazovererr)/float(len(deltazovererr))
# #
# #deltaiovererr_hist, bin_edges= np.histogram(deltaiovererr, weights= iweights, bins=bins2)
# #deltazovererr_hist, bin_edges= np.histogram(deltazovererr, weights= zweights, bins=bins2)
# #
# #plt.figure()
# #plt.plot(bins2[1:], deltaiovererr_hist, label= 'i-band', lw=3, color='orange')
# #plt.plot(bins2[1:], deltazovererr_hist, label= 'z-band', lw=3, color='blue')
# #plt.title('Delta Mag over Mag Error')
# #plt.ylabel('Percent of Total')
# #plt.savefig(outplots +'/'+'DeltaoverERR.pdf')
# #plt.clf()

# #sys.exit()

# print "Number of candidates per ccd"
# #Plot Candidates per CCD #
# #print reals.data

# x = np.zeros(len(np.unique(reals.data.SNID)))
# y = np.unique(reals.data.SNID)
# for i in np.arange(len(x)):
#     x[i] =reals.data.CCDNUM[reals.data.SNID==y[i]][1]

# plt.hist(x, bins= np.arange(min(reals.data.CCDNUM), max(reals.data.CCDNUM) +2 ,1), color='orange')
# plt.title('Hist of real candidates per CCD')
# plt.ylabel('Number of Real Candidates')
# plt.xlabel('CCD Number')
# plt.savefig(outplots +'/'+'Hist_of_real_candidates_per_CCD.pdf')
# plt.clf()

# sys.exit()

# x = np.zeros(len(np.unique(fakes.data.SNID)))
# y = np.unique(fakes.data.SNID)
# for i in np.arange(len(x)):
#     x[i] =fakes.data.CCDNUM[fakes.data.SNID==y[i]][1]


# plt.hist(x, bins= np.arange(min(reals.data.CCDNUM), max(reals.data.CCDNUM) +2,1),color = 'orange')
# plt.title('Hist of fake candidates per CCD')
# plt.ylabel('Number of Fake Candidates')
# plt.xlabel('CCD Number')
# plt.savefig(outplots +'/'+'Hist_of_fake_candidates_per_CCD.pdf')
# plt.clf()

# print "Save candidates info"
# ###Write data files for each candidate including info discussed###


# f1= open(str(outdir)+'/'+'allcandidates.txt', 'w')
# header1 = 'SNID, ' + ' RA, ' + ' DEC, ' + ' CandType,' +  ' NumEpochs, ' + ' NumEpochsml, ' + ' LatestNiteml' 
# f1.write(header1)
# for i in range(0,numofcan):
#     Cand =(reals.data.SNID == urID[i])
#     savedata(reals,urID[i],outdir,triggerid)
# #    if not debug:
# #        continue
#     line = str(urID[i]) + ", " + str(reals.data.RA[Cand][1]) + ", " + str(reals.data.DEC[Cand][1]) + ", " + str(reals.data.CandType[Cand][1]) + ", " + str(reals.data.NumEpochs[Cand][1]) + ", " + str(reals.data.NumEpochsml[Cand][1]) + ", " + str(reals.data.LatestNiteml[Cand][1]) + "\n"
#     table1 = np.array([[int(urID[i]),reals.data.RA[Cand][1],reals.data.DEC[Cand][1],int(reals.data.CandType[Cand][1]),int(reals.data.NumEpochs[Cand][1]),int(reals.data.NumEpochsml[Cand][1]),int(reals.data.LatestNiteml[Cand][1])]])
#     print table1
#     f1.write(line)
#     filename = 'Candidate_'+str(int(urID[i])) + '.txt'
#     htmlfilename = 'Candidate_'+str(int(urID[i])) + '.html'
#     header ='BAND ' + 'x ' + 'y ' + 'Mag ' + 'Nite ' + 'MJD ' + 'Season ' + 'Object ' + 'Exposure ' + 'Field ' + 'CCDNUM'
#     nobs = len(reals.data.BAND[Cand])
#     seasoncol = np.ones((nobs,), dtype = np.int)*int(season)
#     print seasoncol
#     table = np.column_stack((reals.data.BAND[Cand], reals.data.XPIX[Cand], reals.data.YPIX[Cand], reals.data.MAG[Cand], reals.data.OBSNITE[Cand], reals.data.MJD[Cand], seasoncol))
#     np.savetxt(str(outdir)+'/'+filename, table, fmt = '%s', header = header)
#     htmlcode = HTML.table(table.tolist(),header_row = header.split(' '))
#     htmlcode1 = HTML.table(table1.tolist(), header_row = header1.split(', '))
#     f = open(str(outdir)+ '/'+  htmlfilename, 'w')
#     f.write(htmlcode1)
#     f.write(htmlcode)
#     #Collect Stamps for observations of this candidate#
#     thiscand_stampsdir = outstamps  + '/' + str(int(urID[i]))
#     if not os.path.exists(thiscand_stampsdir):
#         os.mkdir(thiscand_stampsdir)
# #Create Stamps_table with 8 columns and nobs rows all empty#
#     stampstable = ([[None]*8])
#     stampsheader = 'Filter, ' + 'Object ID, ' + 'Nite, ' + 'MJD, ' + 'Search, ' + 'Template, ' + 'Difference, ' + 'AutoScan Score,'
#     for j in range(0,nobs):
#         thisobs_nite = str(int(realss.OBSNITE[Cand][j]))
#         thisobs_band = realss.BAND[Cand][j]
#         ccdnum = int(realss.OBSCCDNUM[Cand][j])
#         objexpnum = int(realss.EXPNUM[Cand][j])
#         if int(realss.OBSCCDNUM[Cand][j]) <10:
#             ccdnum = '0' + str(int(realss.OBSCCDNUM[Cand][j]))
# #        expdir = "/data/des41.a/data/marcelle/diffimg/local-runs"
#         thisobs_ID = realss.OBJID[Cand][j]
#         a = expdir + '/' + thisobs_nite + '/'+ str(objexpnum)+'/dp' + str(season) + '/' + thisobs_band + '_' + str(ccdnum) + '/stamps*'
#         print a
#         thisobs_stampsdir = glob.glob(a)[0]
#         print thisobs_stampsdir
# ### MUST UNDO THESE COMMENTS ONCE OBSID IS DEFINED TERM ###
#         filenamediff = thisobs_stampsdir + '/diff' + str(thisobs_ID) + '.gif'
#         filenamesrch = thisobs_stampsdir + '/srch' + str(thisobs_ID) + '.gif'
#         filenametemp = thisobs_stampsdir + '/temp' + str(thisobs_ID) + '.gif' 
#         if thisobs_ID != 0:
#             shutil.copy(filenamediff, thiscand_stampsdir)
#             shutil.copy(filenamesrch, thiscand_stampsdir)
#             shutil.copy(filenametemp, thiscand_stampsdir)
#             path1 = thiscand_stampsdir + '/srch' + str(thisobs_ID) + '.gif'
#             print path1
#             search= image('', 'stamps/' + str(int(urID[i]))  + '/srch' + str(thisobs_ID) + '.gif')
#             temp  = image('', 'stamps/' + str(int(urID[i])) + '/temp' + str(thisobs_ID) + '.gif')
#             diff  = image('', 'stamps/' + str(int(urID[i])) + '/diff' + str(thisobs_ID) + '.gif')
            
#         if thisobs_ID == 0:
#             search = 'no search'
#             temp = 'no temp'
#             diff = 'no diff'
# #Replace in the empty spaces in table with values/pictures#
#         stampstable[j][0] = realss.BAND[Cand][j]
#         stampstable[j][1] = realss.OBJID[Cand][j]
#         stampstable[j][2] = realss.OBSNITE[Cand][j]
#         stampstable[j][3] = realss.MJD[Cand][j]
#         stampstable[j][4] = search
#         stampstable[j][5] = temp
#         stampstable[j][6] = diff
#         stampstable[j][7] = realss.PHOTPROB[Cand][j]
#         stampstable.append([None] * 8)
        
#     htmlcode2 = HTML.table(stampstable, header_row= stampsheader.split(', '))
#     f.write(htmlcode2)
#     #Making Plot of Flux vs MJD for each Candidate#
#     fvmplottable = ([[None]*len(ubands)])
#     for b in range(0, len(ubands)):
#         Bandcand = realss.BAND[Cand]
#         Band = Bandcand == ubands[b]
#         Flux = realss.FLUXCAL[Cand][Band]
#         MJD = realss.MJD[Cand][Band]
#         Fluxerr = realss.FLUXCALERR[Cand][Band]
#         plt.scatter(MJD,Flux, color = colors[b])
#         plt.errorbar(MJD,Flux, yerr=Fluxerr,color=colors[b], ls = 'none')
#         plt.xlabel('MJD')
#         plt.ylabel('Flux')
#     plt.title('Flux vs. MJD for candidate'  + str(int(urID[i]))+ ' iband=r, zband=g')
#     plt.savefig(outdir +'/plots/lightcurves/FluxvsMJD_for_cand_' + str(int(urID[i]))+'.png')
#     plt.clf()
#     plotpic = image('', 'plots/lightcurves/FluxvsMJD_for_cand_' + str(int(urID[i])) + '.png')
#     fvmplottable[0][b] = plotpic
#     htmlcode3 = HTML.table(fvmplottable)
#     f.write(htmlcode3)

                     
# f1.close()  
    
# print "SUCCESS"    


             


# #Save output files in specific and organized directories#

# ###Setup a display of the search, template and difference images for each candidate. Automatically save for top n candidates. Give option to user to pick candidate by ID to display###
