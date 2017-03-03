#!/bin/bash

##source /sps/lsst/Library/stack_v11_0/loadLSST.bash
#export LSSTSW=/sps/lsst/Library/lsstsw
#export EUPS_PATH=$LSSTSW/stack
#source $LSSTSW/bin/setup.sh
##setup pipe_tasks
#setup sconsUtils 10.0
#C++11

export PATH=/opt/rh/devtoolset-3/root/usr/bin:${PATH}

#export LSSTSW=/sps/lsst/Library/lsstsw
export LSSTSW=/sps/lsst/Library/new/lsstsw
export EUPS_PATH=$LSSTSW/stack
source $LSSTSW/bin/setup.sh

export PYTHONPATH=../lib/python2.7/site-packages:$PYTHONPATH

setup -k -r ../pykg_config
#setup -k -r ../healpy
setup -k -r ../python_future
setup -k -r ../sims_data
setup -k -r ../sims_utils
#setup -k -r ../sims_catalogs_generation
setup -k -r ../sims_sed_library
setup -k -r ../sims_dustmaps
setup -k -r ../throughputs
setup -k -r ../sims_maps
setup -k -r ../sims_photUtils
setup -k -r ../sims_coordUtils
setup -k -r ../sqlalchemy
setup -k -r ../sims_catalogs
setup -k -r ../sims_maf
#setup -k -r ../sims_operations