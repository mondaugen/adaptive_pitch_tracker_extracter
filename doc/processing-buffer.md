# A buffer for processing

Often there might be a mis-match between the number of samples fed into a
routine each processing block and the number of samples required to make a
decision based on some analysis. For this a buffer is required that can hold
some number of samples required for analysis.

This is very simple. Writes and reads must happen at the same rate if the buffer
is not to overflow. If we assume the pattern write -> analyze -> read is always
followed, then for analysis window size A and read/write size B, the ringbuffer
is initialized with max(A-B,0) samples and a maximum size of max(A,B). Then B
samples are written, A samples analyzed, and B samples read out from the
beginning of the ring-buffer.
