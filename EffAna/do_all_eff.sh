root -b -l -q CrossSectionAnaMult.C'(1,0,0)'
wait
root -b -l -q CrossSectionAnaMult.C'(1,1,0)'
wait
root -b -l -q CrossSectionAnaMult.C'(1,0,1)'
wait
root -b -l -q CrossSectionAnaMult.C'(1,1,1)'
wait
rm -r BP/EffFinal
rm -r Bs/EffFinal
wait
root -b -l -q CrossSectionAnaY.C'(1,0,0)'
wait
root -b -l -q CrossSectionAnaY.C'(1,1,0)'
wait
root -b -l -q CrossSectionAnaY.C'(1,0,1)'
wait
root -b -l -q CrossSectionAnaY.C'(1,1,1)'


