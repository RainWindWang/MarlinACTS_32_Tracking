echo "evnironment setup..."
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-04-12
source /Users/yufengwang/Work/key4hep/ACTs/install/bin/this_acts.sh
export CMAKE_PREFIX_PATH=/Users/yufengwang/Work/key4hep/ACTs/install:$CMAKE_PREFIX_PATH
export MARLIN_DLL=$MARLIN_DLL:/Users/yufengwang/Work/key4hep/LUXE/MarlinACTS_32_Tracking/install/lib/libMarlinACTSTracking.so
