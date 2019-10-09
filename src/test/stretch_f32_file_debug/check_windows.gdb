define signal_stretcher_f32_process_break
    set $window_size=window_size
    eval "print *(float(*)[%d])signal_stretcher->pvs->analysis_window", $window_size
    eval "print *(float(*)[%d])signal_stretcher->pvs->synthesis_window", $window_size
end

define signal_stretcher_f32_dump_windows
    set $window_size=window_size
    eval "dump binary memory /tmp/an.f32 signal_stretcher->pvs->analysis_window &((float*)signal_stretcher->pvs->analysis_window)[%d]", $window_size
    eval "dump binary memory /tmp/sy.f32 signal_stretcher->pvs->synthesis_window &((float*)signal_stretcher->pvs->synthesis_window)[%d]", $window_size
end

break signal_stretcher_f32_process
commands
    finish
    signal_stretcher_f32_dump_windows
end

run
