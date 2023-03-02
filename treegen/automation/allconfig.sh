#!/ bin / bash

cd allLong
sh runbatch.sh
cd ../allShort
sh runbatch.sh
cd ../internalLong
sh runbatch.sh
cd ../internalShort
sh runbatch.sh
cd ..
python3 aggregateTopo.py