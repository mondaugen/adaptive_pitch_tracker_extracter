set pagination off

break spec_diff_finder_find
commands 1
    finish
    step
    dump memory /tmp/sd.f32 sd->spec_diff sd->specdiff+sd->length
    continue
end

break iir_avg
commands 2
    finish
    step
    dump memory /tmp/sd_smooth.f32 sd->spec_diff sd->specdiff+sd->length
    continue
end
