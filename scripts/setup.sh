echo "evnironment setup..."
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-04-12
source /nfs/dust/ilc/user/wangyufe/ACTs/install/bin/this_acts.sh
#source /Users/yufengwang/Work/key4hep/ACTs/install/bin/this_acts.sh
export CMAKE_PREFIX_PATH=/nfs/dust/ilc/user/wangyufe/ACTs/install:$CMAKE_PREFIX_PATH
#export CMAKE_PREFIX_PATH=/Users/yufengwang/Work/key4hep/ACTs/install:$CMAKE_PREFIX_PATH
export MARLIN_DLL=$MARLIN_DLL:/nfs/dust/ilc/user/wangyufe/MarlinACTS_32_Tracking/install/lib/libMarlinACTSTracking.so
#export MARLIN_DLL=$MARLIN_DLL:/Users/yufengwang/Work/key4hep/LUXE/MarlinACTS_32_Tracking/install/lib/libMarlinACTSTracking.so
