# example processing shell script for bash on unix/linux systems
# Daniel Buscombe, July 2014

# compile cython modules
bash compile_pyhum.sh

# test data directory
k=$PWD'/test_data'

# convert raw data
python pyhum_read.py -i test.DAT -s $k

# radiometrically correct
python pyhum_correct.py -i test.DAT -s $k

# estimate texture lengtn scales
python pyhum_texture.py -i test.DAT -s $k

