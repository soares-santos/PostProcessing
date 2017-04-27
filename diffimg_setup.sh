
##### edit these lines before running the code:
#export ROOTDIR=/pnfs/des/persistent/gw
#export ROOTDIR2=/data/des41.a/data/marcelle/postproc/test
#export SEASON=300
#####


#for IFDH
export EXPERIMENT=des
export PATH=${PATH}:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/cpn/v1_7/NULL/bin:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/ifdhc/v1_8_10/Linux64bit-2-6-2-12/bin
export PYTHONPATH=/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/ifdhc/v1_8_10/Linux64bit-2-6-2-12/lib/python:${PYTHONPATH}
export IFDH_NO_PROXY=1
export IFDHC_LIB=/cvmfs/fermilab.opensciencegrid.org/products/common/prd/ifdhc/v1_8_10/Linux64bit-2-6-2-12/lib
export IFDH_CP_MAXRETRIES=2
#/cvmfs/grid.cern.ch/util/cvmfs-uptodate /cvmfs/des.opensciencegrid.org
source /cvmfs/des.opensciencegrid.org/2015_Q2/eeups/SL6/eups/desdm_eups_setup.sh
source /cvmfs/des.opensciencegrid.org/eeups/startup.sh
export EUPS_PATH=/cvmfs/des.opensciencegrid.org/eeups/fnaleups:$EUPS_PATH

#other setups
setup perl 5.18.1+6 # || exit 134
setup Y2Nstack 1.0.6+18
setup diffimg #gw0 
setup ftools v6.17 
export HEADAS=$FTOOLS_DIR
setup autoscan
setup easyaccess
setup extralibs 1.0
#setup html
echo "EUPS setup complete"

export DES_SERVICES=${HOME}/.desservices.ini 
export DES_DB_SECTION=db-sn-test
export DIFFIMG_HOST=FNAL
export SCAMP_CATALOG_DIR=$PWD/SNscampCatalog
export AUTOSCAN_PYTHON=$PYTHON_DIR/bin/python
export DES_ROOT=/data/des20.b/data/SNDATA_ROOT/INTERNAL/DES 
export TOPDIR_SNFORCEPHOTO_IMAGES=${ROOTDIR}/forcephoto/images/dp${SEASON} 
export TOPDIR_SNFORCEPHOTO_OUTPUT=${ROOTDIR2}/forcephoto/output/dp${SEASON}   
export TOPDIR_DATAFILES_PUBLIC=${ROOTDIR}/DESSN_PIPELINE/SNFORCE/DATAFILES_TEST
export TOPDIR_WSTEMPLATES=${ROOTDIR}/WSTemplates
export TOPDIR_TEMPLATES=${ROOTDIR}/WSTemplates
export TOPDIR_SNTEMPLATES=${ROOTDIR}/SNTemplates
export TOPDIR_WSRUNS=${ROOTDIR}/data/WSruns
export TOPDIR_SNRUNS=${ROOTDIR}/data/SNruns

# these vars are for the make pair function that we pulled out of makeWSTemplates.sh
TOPDIR_WSDIFF=${TOPDIR_WSTEMPLATES}
DATADIR=${TOPDIR_WSDIFF}/data             # DECam_XXXXXX directories
CORNERDIR=${TOPDIR_WSDIFF}/pairs          # output XXXXXX.out and XXXXXX-YYYYYY.out
ETCDIR=${DIFFIMG_DIR}/etc                 # parameter files
CALDIR=${TOPDIR_WSDIFF}/relativeZP        # relative zeropoints
MAKETEMPLDIR=${TOPDIR_WSDIFF}/makeTempl   # templates are made in here

XY2SKY=${WCSTOOLS_DIR}/bin/xy2sky
AWK=/bin/awk
export PFILES=${PWD}/syspfiles

##mkdir -p $TOPDIR_SNFORCEPHOTO_IMAGES $DES_ROOT $TOPDIR_SNFORCEPHOTO_OUTPUT $TOPDIR_DATAFILES_PUBLIC

export SNANA_DIR=/data/des41.b/data/kessler/snana/snana
export SNANA_ROOT=/data/des41.b/data/SNDATA_ROOT

## use Ken's development version of the diffimg code:
#export DIFFIMG_DIR=/data/des40.b/data/kherner/Diffimg-devel/diffimg-trunk
#export PATH=$DIFFIMG_DIR/bin:$PATH
