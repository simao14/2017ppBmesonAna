DOANALYSISPbPb_FULL_BP=0
DOANALYSISPbPb_BINNED_PT_BP=1
DOANALYSISPbPb_BINNED_Y_BP=0
DOANALYSISPbPb_BINNED_MULTI_BP=0

DOANALYSISPbPb_MATCHING_Bs_BINS_BP=0

DOANALYSISPbPb_BINNED_PT_BP_TRK=0
DOANALYSISPbPb_BINNED_Y_BP_TRK=1
DOANALYSISPbPb_BINNED_MULTI_BP_TRK=0

#Data and MC Samples
MC_BP="/data3/tasheng/presel/BPMC_nom.root"
Data_BP="/data3/tasheng/presel/BPData_nom.root"
#MC_BP="/data3/smcosta/data/BPMC_nom.root"
#Data_BP="/data3/smcosta/data/BPData_nom.root"

INPUTJPSI="/data3/tasheng/presel/jpsinp_nom.root"
#INPUTJPSI="/data3/smcosta/data/jpsinp_nom.root"
#Data and MC Samples

#NEW NMB from https://twiki.cern.ch/twiki/pub/CMS/HINUpsilonRaa2016/Jason_MinBiasCounting_2017-02-02.pdf

CUTPbPb="Bpt>0"
cut_trk_tight="(track>1)"

mkdir -p ROOTfiles/
OutputFile_BP_FULL="ROOTfiles/yields_Bp_full.root"
OutputFile_BP_BINNED_Y="ROOTfiles/yields_Bp_binned_y.root"
OutputFile_BP_BINNED_Mult="ROOTfiles/yields_Bp_binned_Mult.root"
OutputFile_BP_BINNED_PT="ROOTfiles/yields_Bp_binned_pt.root"
OutputFile_BP_BINNED_PT_trk="ROOTfiles/yields_Bp_binned_pt_trk.root"
OutputFile_BP_BINNED_Y_trk="ROOTfiles/yields_Bp_binned_y_trk.root"
OutputFile_BP_BINNED_Mult_trk="ROOTfiles/yields_Bp_binned_Mult_trk.root"

NPROOFIT_PbPb="yes"

if [ $DOANALYSISPbPb_FULL_BP  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','1','\"$Data_BP\"','\"$MC_BP\"','\"Bpt\"','\"$CUTPbPb\"','\"$OutputFile_BP_FULL\"','\"results/BP\"','\"$NPROOFIT_PbPb\"','\"$INPUTJPSI\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_PT_BP  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"Bpt\"','\"$CUTPbPb\"','\"$OutputFile_BP_BINNED_PT\"','\"results/BP/Bpt\"','\"$NPROOFIT_PbPb\"','\"$INPUTJPSI\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_PT_BP_TRK -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"Bpt\"','\"$cut_trk_tight\"','\"$OutputFile_BP_BINNED_PT_trk\"','\"results/BP/trk_tight\"','\"$NPROOFIT_PbPb\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_Y_BP  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"By\"','\"$CUTPbPb\"','\"$OutputFile_BP_BINNED_Y\"','\"results/BP/By\"','\"$NPROOFIT_PbPb\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_Y_BP_TRK  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"By\"','\"$cut_trk_tight\"','\"$OutputFile_BP_BINNED_Y_trk\"','\"results/BP/By\"','\"$NPROOFIT_PbPb\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_MULTI_BP  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"nMult\"','\"$CUTPbPb\"','\"$OutputFile_BP_BINNED_Mult\"','\"results/BP/nMult\"','\"$NPROOFIT_PbPb\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_MULTI_BP_TRK  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"nMult\"','\"$cut_trk_tight\"','\"$OutputFile_BP_BINNED_Mult_trk\"','\"results/BP/nMult\"','\"$NPROOFIT_PbPb\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi