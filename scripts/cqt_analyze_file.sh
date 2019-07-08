sox $1 -t f64 -c 1 -r 16k /tmp/guit.f64 trim "$2" "$3"
[ -z $P_MIN ] && P_MIN=38
export P_MIN
PYTHONPATH=".:$PYTHONPATH" python3 scripts/cqt_file.py&
sox -t f64 -c 1 -r 16k /tmp/guit.f64 -d repeat 10000000
