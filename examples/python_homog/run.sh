#clean.sh
ARGS=""
lastrestart=$(ls -rt restart_*.pkl | tail -n 1)

if [[ "$lastrestart" != "" ]]
then
   ARGS+=" -r $lastrestart"
fi

run_mc.py $ARGS -d 400 -i 1.0 --init-args='asst=0.0,dx=40./1500,n=500' --solver-args='dissipation=.12' > output.log 2>&1
#plot.py
python -m python.read data/ diags.pkl out.nc
