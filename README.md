# honeycomb_scinti
Detector simulation with honeycomb structure for reactor neutrino 

## Compile
mkdir build
cd build/
cmake ../B4a
make

## Execute
cd macro
../build/exampleB4 -m run1.mac -f <filename>

## Convert data
root -l -q -b <filename> draw.C

## Analize data
root -l ana_<filename> plot.C
