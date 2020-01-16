# check number of blocks calculation

def N_blocks(N,H,W):
    return max(0,(N-W)//H + 1)

H1=100
W1=500

# 0 blocks should fit
N1=499
print("true 0")
print(N_blocks(N1,H1,W1))

# 1 block should fit
N2=500
print("true 1")
print(N_blocks(N2,H1,W1))

# 1 block should fit
N3=599
print("true 1")
print(N_blocks(N3,H1,W1))

# 2 blocks should fit
N4=600
print("true 2")
print(N_blocks(N4,H1,W1))

# 2 blocks should fit
N5=100
print("true 0")
print(N_blocks(N5,H1,W1))

