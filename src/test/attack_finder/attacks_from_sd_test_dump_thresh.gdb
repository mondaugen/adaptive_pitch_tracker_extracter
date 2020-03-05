set pagination off

break discount_local_max_f32
commands
    finish
    step
    dump memory /tmp/thresh.f32 thresh thresh+sd->length
    continue
end 
