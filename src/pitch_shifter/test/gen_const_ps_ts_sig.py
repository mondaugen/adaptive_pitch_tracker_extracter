import common
import numpy as np

PS=common.get_env('PS',default=1,conv=float)
TS=common.get_env('TS',default=1,conv=float)
N=common.get_env('N',default=100000,conv=int)

ps=np.zeros(N,dtype='uint32')
ps[:]=int(round(PS*(1<<16)))

ts=np.zeros(N,dtype='int32')
ts[:]=int(round(TS*(1<<16)))

ps.tofile('/tmp/ps_rate.u16q16')
ts.tofile('/tmp/ts_rate.s16q16')
