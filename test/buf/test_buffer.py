import numpy as np
from buf import inf_buf

s=np.arange(0,-5,-1,dtype='float')

ib=inf_buf(s,start=2,fill=1.)

def test_lookup():
    assert ib[1] == 1.
    assert ib[2] == 0.
    assert ib[3] == -1.
    assert ib[7] == 1.
    assert ib[8] == 1.
    assert np.all(ib[1:4] == np.array([1,0,-1],dtype='float'))
    assert np.all(ib[3:9] == np.array([-1,-2,-3,-4,1,1],dtype='float'))
    assert np.all(ib[4:7] == np.array([-2,-3,-4],dtype='float'))
    assert np.all(ib[0:8] == np.array([1,1,0,-1,-2,-3,-4,1],dtype='float'))
    assert np.all(ib[-2:2] == np.array([1,1,1,1],dtype='float'))
    assert np.all(ib[7:9] == np.array([1,1],dtype='float'))

