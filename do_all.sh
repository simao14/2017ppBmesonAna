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
    bash bpdoRoofit.sh &
    bash bsdoRoofit.sh &
    wait
    popd
}

bpEff () {
    pushd EffAna
    echo "Takes BPw.root as input"
    ls -l BDTWeights/BPw.root

    root -b -l -q MCEff.C'(1,0,0)' > effbp.log                                          
    wait
    root -b -l -q CrossSectionAnaMult.C'(1,1,0)'  #(TnPon =1 noTnP=0 , pT=0,y=1,mult=2 , data =0 mc = 1)
    root -b -l -q CrossSectionAnaMult.C'(1,0,0)'
    popd
}

bsEff () {
    pushd EffAna
    echo "Takes Bsw.root as input"
    ls -l BDTWeights/Bsw.root

    root -b -l -q MCEff.C'(1,0,1)' > effbs.log                                          
    wait
    root -b -l -q CrossSectionAnaMult.C'(1,1,1)'  #(TnPon =1 noTnP=0 , pT=0,y=1,mult=2 , data =0 mc = 1)
    root -b -l -q CrossSectionAnaMult.C'(1,0,1)'
    popd
}

syst2D () {
    pushd 2DMapSyst
    root -b -l -q CalEffSystB.C'(0,0)'                      # >> outfiles/BPsyst2d.root
    root -b -l -q CalEffSystB.C'(0,1)'
    root -b -l -q CalEffSystB.C'(1,0)'                      # >> outfiles/Bssyst2d.root
    root -b -l -q CalEffSystB.C'(1,1)'             
    root -b -l -q PlotEffSyst2D.C'(0,0)'                                                          # << outfiles/BPsyst2d.root
    root -b -l -q PlotEffSyst2D.C'(0,1)'
    root -b -l -q PlotEffSyst2D.C'(1,0)'                                                          # << outfiles/BPsyst2d.root
    root -b -l -q PlotEffSyst2D.C'(1,1)'
    popd
}

syst1D () {
    pushd 1DMapSyst
    root -b -l -q CalEffSystB.C'(0,0)'                      # >> outfiles/BPsyst2d.root
    root -b -l -q CalEffSystB.C'(0,1)'
    root -b -l -q CalEffSystB.C'(1,0)'                      # >> outfiles/Bssyst2d.root
    root -b -l -q CalEffSystB.C'(1,1)'             
    root -b -l -q PlotEffSyst1D.C'(0,0)'                                                          # << outfiles/BPsyst2d.root
    root -b -l -q PlotEffSyst1D.C'(0,1)'
    root -b -l -q PlotEffSyst1D.C'(1,0)'                                                          # << outfiles/BPsyst2d.root
    root -b -l -q PlotEffSyst1D.C'(1,1)'
    popd
}

## MC Stat Systematics
bpStat () {
    pushd MCStatSyst/BP
    root -b -l -q Generate2DMaps.C
    root -b -l -q MCStatCal.C > mcstat.log
    popd
}
bsStat() {
    cd MCStatSyst/Bs
    root -b -l -q Generate2DMaps.C
    root -b -l -q MCStatCal.C > mcstat.log
    cd ../..
}

comp () {
    # get pdf variation errors
    python master.py "Bpt"
    python master.py "By"
    
    # Get pre-selection error
    python comppre.py                     #<----------------------- NOT RUNNING (FILE FROM CODE MISSING)

    cd BsBPFinalResults/BsBPRatio/
    root -b -l -q PlotBsBPRatio.C'(1)'
    cd ..

    cd Comparisons/Fiducial/
    root -b -l -q Bmeson_Comparisons.C'(0,0)'
    root -b -l -q Bmeson_Comparisons.C'(1,0)'
    root -b -l -q Bmeson_Comparisons.C'(0,1)'
    root -b -l -q Bmeson_Comparisons.C'(1,1)'

    python syst_table.py                  #<----------------------- NOT RUNNING
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

#maketnp
#bptshape




#yield
#wait

#bpEff &
#bsEff &
#wait

syst2D
wait
syst1D
wait

bpStat&
bsStat&
wait

comp
paperPlots
