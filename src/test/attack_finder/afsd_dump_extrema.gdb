define dumpit
 dump binary memory /tmp/sd_mins.u32 sd_mins sd_mins+n_sd_mins
 dump binary memory /tmp/sd_maxs.u32 sd_maxs sd_maxs+n_sd_maxs
 dump binary memory /tmp/spec_diff.f32 sd->spec_diff sd->spec_diff+sd->length
end
