source /swshare/psit3/etc/profile.d/cms_ui_env.sh
source $VO_CMS_SW_DIR/cmsset_default.sh
cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src
cmsenv
git clone https://github.com/clseitz/ttbbAnalysis.git
scramv1 b
cd ttbbAnalysis/KinFitter/test
python runStuff.py