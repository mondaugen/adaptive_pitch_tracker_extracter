" stuff that makes working here easier

" fill the line with only the environment variable name and then this will put
" the get_env function around it
nmap ,ge ^"zyiw$a=common.get_env('"zpa')i

" include present directory for python completion to work
python3 import sys
python3 sys.path.append(".")
