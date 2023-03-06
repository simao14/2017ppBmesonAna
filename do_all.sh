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
    wait
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

    root -b -l -q MCEff.C'(1,0)' > eff.log # >> bpsyst2d.root
    wait
    root -b -l -q CrossSectionAna.C'(1)'

    # root -b -l -q CrossSectionAnaMult.C'(1)'
    # >> BP/EffAna/FinalFiles/BPPPCorrYieldPT.root
    popd
}

bsEff () {
    pushd Bs/EffAna
    echo "Takes Bsw.root as input"
    ls -l BDTWeights/Bsw.root

    root -b -l -q MCEff.C'(1,0)' > eff.log
    wait
    root -b -l -q CrossSectionAna.C'(1)'

    # root -b -l -q CrossSectionAnaMult.C'(1)'
    # >> Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root
    popd
}

nominal () {
    # central value
    cd BsBPFinalResults/BsBPRatio
    root -b -l -q PlotBsBPRatio.C'(1)'
    cd ../RAA
    root -b -l -q BPRAA.C
    root -b -l -q BsRAA.C
    cd ..

    cd Comparisons/Fiducial
    root -b -l -q Bmeson_Comparisons.C'(0)'
    root -b -l -q Bmeson_Comparisons.C'(1)'
    root -b -l -q BPNewFidNoScale.C
    root -b -l -q BsNewFidNoScale.C
    cd ../../..
}

syst () {
    pushd 2DMapSyst
    root -b -l -q CalEffSystBP.C # >> outfiles/bpsyst2d.root
    root -b -l -q CalEffSystBs.C
    root -b -l -q PlotEffSyst2D.C'(0)'
    root -b -l -q PlotEffSyst2D.C'(1)'
    # need to pass the above to cross calculation
    popd
}

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

comp () {
    # get pdf variation errors
    python master.py
    # Get pre-selection error
    python comppre.py
    # comparison plot again
    pushd BsBPFinalResults/Comparisons/Fiducial/
    # << BP/EffAna/FinalFiles/BPPPCorrYieldPT.root
    root -b -l -q Bmeson_Comparisons.C'(0)'
    root -b -l -q Bmeson_Comparisons.C'(1)'
    python syst_table.py
    cd ../../RAA/
    root -b -l -q BPRAA.C
    root -b -l -q BsRAA.C
    popd
}

paperPlots () {
    # input
    pushd MakeFinalPlots/NominalPlots/CrossSection
    root -b -l -q plotPt.C'(1,1,0,1,1)'
    cd ../RAA
    root -b -l -q plotPt.C'(1,1,0,1,1)'
    popd
}

#UNCOMMENT ACORDINGLY
#(Run by THIS ORDER!)

bptshape

yield
wait

bpEff &
bsEff &
wait

#syst
#wait

#nominal

# bpStat&
# bsStat&
# wait

#comp
#paperPlots
