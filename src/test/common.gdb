# print an array of float32
# arguments are expression-giving-pointer length
define p_f32
    eval "print *(float(*)[%d])%s", $arg1, "$arg0"
end

# dump array of float32
# arguments are expression-giving-pointer length output-file
define d_f32
    eval "dump binary memory %s %s &((float*)%s)[%d]", "$arg2", "$arg0", "$arg0", $arg1
end

# dump array of complex64 (two packed f32 floats) representing real data
# length is the array of the real sequence, so the complex sequence has length/2+1 values
# arguments are expression-giving-pointer length output-file
define d_rz64
    # twice the number of complex values, because we use the float32 routine to dump
    set $length=2*($arg1/2+1)
    d_f32 $arg0 $length $arg2
end

# append array of float32
# arguments are expression-giving-pointer length output-file
define a_f32
    eval "append binary memory %s %s &((float*)%s)[%d]", "$arg2", "$arg0", "$arg0", $arg1
end

# append array of complex64 (two packed f32 floats) representing real data
# length is the array of the real sequence, so the complex sequence has length/2+1 values
# arguments are expression-giving-pointer length output-file
define a_rz64
    # twice the number of complex values, because we use the float32 routine to append
    set $length=2*($arg1/2+1)
    a_f32 $arg0 $length $arg2
end
