# Before the divide is done, dump the 2 arrays to see what's going on.

source src/test/common.gdb

break pvoc_synth.c:106
commands
    printf "window length: %d", pvs->config.window_length

    # dump values
    d_rz64 pvs->z_input0 pvs->config.window_length /tmp/z_input0.z64
    d_rz64 pvs->z_inputH pvs->config.window_length /tmp/z_inputH.z64

    # remove old commands
    commands 1
    end

    # now append
    commands 1
    a_rz64 pvs->z_input0 pvs->config.window_length /tmp/z_input0.z64
    a_rz64 pvs->z_inputH pvs->config.window_length /tmp/z_inputH.z64
    continue
    end

    continue
end

run
