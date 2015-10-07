# $1 is w
# $2 is initial guess of shooting parameter
# $3 is r1
# $4 is rmax

echo 1.0 0.9 $1 $2 $3 $4 > parameters

python run.py

mkdir -p data/w=$1
mv {pars,integral-data,functions}.dat data/w=$1/
mv {f,g}.pdf data/w=$1/
mv parameters data/w=$1/
