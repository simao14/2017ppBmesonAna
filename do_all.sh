#! /bin/bash

maketnp () {
    ## run TnP and attach them to the unskimmed MC files
    pushd MakeMCTnP/
    root -q -b -l TnPWeight.C'(0)' &
    root -q -b -l TnPWeight.C'(1)' &
    wait
    popd
}

# new bpt reweighting uses FONLL and doesn't depend on acceptance
# function >> MCEff.C
bptshape () {
    pushd NewBptStudies
    python ReweightY.py &
    wait
    root -q -b -l ReweightBpt.C'(0)' &
    root -q -b -l ReweightBpt.C'(1)' &
    popd
}


yield () {
    ## yield extraction
    pushd henri2022
    # roofitB.C contains 2 By cuts
    bash bpdoRoofit.sh &
    bash bsdoRoofit.sh &
    wait
    popd
}

bpEff () {
    pushd BP/EffAna
    echo "Takes BPw.root as input"
    ls -l BDTWeights/BPw.root

    #root -b -l -q MCEff.C'(1,0)' > eff.log # >> bpsyst2d.root
    wait
    popd
    root -b -l -q CrossSectionAna.C'(1,0)'                              #UNIFY

    #root -b -l -q CrossSectionAnaMult.C'(1,0)'
    # >> BP/EffAna/FinalFiles/BPPPCorrYieldPT.root
}
bsEff () {
    pushd Bs/EffAna
    echo "Takes Bsw.root as input"
    ls -l BDTWeights/Bsw.root

    #root -b -l -q MCEff.C'(1,0)' > eff.log
    wait
    popd
    root -b -l -q CrossSectionAna.C'(1,1)'                               #UNIFY

    #root -b -l -q CrossSectionAnaMult.C'(1,0)'
    # >> Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root
    
}

syst () {
    pushd 2DMapSyst
    root -b -l -q CalEffSystB.C'(0)'                      # >> outfiles/BPsyst2d.root
    root -b -l -q CalEffSystB.C'(1)'                      # >> outfiles/Bssyst2d.root             
    root -b -l -q PlotEffSyst2D.C'(0)'                                                          # << outfiles/BPsyst2d.root
    root -b -l -q PlotEffSyst2D.C'(1)'                                                          # << outfiles/BPsyst2d.root
    popd
}

################# NOT SURE WHAT THIS IS FOR
## MC Stat Systematics
# takes 2D map eff as input
bpStat () {
    pushd MCStatSyst/BP
    root -b -l -q Generate2DMaps.C
    # more than 2hr
    root -b -l -q MCStatCal.C > mcstat.log
    popd
}
bsStat() {
    cd MCStatSyst/Bs
    root -b -l -q Generate2DMaps.C
    # ~1hr
    root -b -l -q MCStatCal.C > mcstat.log
    cd ../..
}
################# NOT SURE WHAT THIS IS FOR

comp () {
    # get pdf variation errors
    python master.py
    # Get pre-selection error
    python comppre.py                     #<-----------------------NOT RUNNING (FILE FROM CODE MISSING)

    cd BsBPFinalResults/BsBPRatio/
    root -b -l -q PlotBsBPRatio.C'(1)'
    cd ..

    cd Comparisons/Fiducial/
    root -b -l -q Bmeson_Comparisons.C'(0)'
    root -b -l -q Bmeson_Comparisons.C'(1)'
    #root -b -l -q BPNewFidNoScale.C                       #Unify w priotiy
    #root -b -l -q BsNewFidNoScale.C                       #Unify w priotiy (Bs results are visualy ugly)
    python syst_table.py                  #<-----------------------NOT RUNNING
    cd ../..

    cd RAA/
    root -b -l -q BPRAA.C                                 #UNIFY
    root -b -l -q BsRAA.C                                 #UNIFY
    cd ../..                                            
}

paperPlots () {
    # input
    pushd MakeFinalPlots/NominalPlots/CrossSection
    root -b -l -q plotPt.C'(1,1,0,1,1)'                 ########## UNIFY THESE WITH A LOT OF CARE
    cd ../RAA
    root -b -l -q plotPt.C'(1,1,0,1,1)'                 ########## UNIFY THESE WITH A LOT OF CARE
    popd
}

#UNCOMMENT ACORDINGLY
#(Run by THIS ORDER!)

#makeTnP
#bptshape

#yield
#wait

bpEff &
bsEff &
wait

#syst
#wait

# bpStat&
# bsStat&
# wait

#comp
#paperPlots
