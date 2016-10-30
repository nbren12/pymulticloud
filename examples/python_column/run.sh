pushd ../../ ; ./compile.sh ; popd

cat << EOF > input.nml
&MCPARM
 alpha_s = .125,
 xic = -0.5,
 xis = 0.0,
 mu = 0.10,
/


EOF

ipython column\ model.py
