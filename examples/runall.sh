ls -d python* | xargs -I{} echo "pushd {}; [[ ! -e report.pdf ]] && plot.py; popd;" | parallel -j 10

OUT=~/workspace/multicmt/doc/2016-10-24

cp python/report.pdf $OUT/fmk_walker.pdf
cp python_cmt/report.pdf $OUT/cmt_walker.pdf
cp python_nonlinear/report.pdf $OUT/nl_walker.pdf
cp python_nonlinear_gms/report.pdf $OUT/nl08_walker.pdf
