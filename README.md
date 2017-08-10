# DTDPGAnalysis
Code for CMS DT offline analysis
This package contains the code needed for CMS DT Prompt Offline Analysis and for DT root-ple production.

To install it and run DTNtuple production:

```
cmsrel CMSSW_8_0_29
cd CMSSW_8_0_29/src/

cmsenv

git clone https://github.com/battibass/DTDPGAnalysis.git UserCode/DTDPGAnalysis

cd UserCode/DTDPGAnalysis/

git checkout dt_rpc_trigger_bari

cd ../../

scramv1 b -j 5

cd UserCode/DTDPGAnalysis/test

cmsRun RunTree_collisions_cfg.py 

```

