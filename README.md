# MarlinACTSTracking
Marlin-based ACTS (v32) Tracking adapted from MuonColliderSoft

For automatic conversion of DD4hep geometry to ACTS Tracking geometry, the tracker needs to be defined as barrel/endcap. Use LUXETrackerAsEndcap.xml from `luxegeo` for the DD4hep geometry description.

To setup, do
```bash
source /cvmfs/ilc.desy.de/key4hep/luxe_setup.sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=../install
make install
cd ..
export MARLIN_DLL=${MARLIN_DLL}${PWD}'/install/lib/libMarlinACTSTracking.so:'
```
To run combinatorial Kalman Filter tracking:
```bash
Marlin actsseedckf_steer.xml
```

Local setup With ACTS v32 installed:
```bash
git clone git@github.com:RainWindWang/MarlinACTS_32_Tracking.git
cd MarlinACTS_32_Tracking
git checkout -b porting_v32 origin/porting_v32
source /cvmfs/sw.hsf.org/key4hep/setup.sh
source <path/to/acts>/install/bin/this_acts.sh
mkdir build
cd build
export CMAKE_PREFIX_PATH=<path/to/acts>/install:$CMAKE_PREFIX_PATH
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=../install
make install
cd ..
export MARLIN_DLL=${MARLIN_DLL}${PWD}'/install/lib/libMarlinACTSTracking.so:'
```
