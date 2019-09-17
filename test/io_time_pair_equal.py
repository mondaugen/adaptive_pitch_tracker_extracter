from time_map_tstretch import *

a=io_time_pair(1,2)
b=io_time_pair(1,2)
c=io_time_pair(2,3)

assert(a==b)
assert(a==a)
assert(a!=c)
