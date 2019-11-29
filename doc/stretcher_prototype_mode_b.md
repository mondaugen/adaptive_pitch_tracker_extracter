# Mode b: Input has been recorded, playback triggered, time-stretch/pitch-shift in real-time

## Parts to develop

- Attack detector
    - Spectral difference computer
    - Discounted max finder
- A way to play doing both time-stretching, pitch-shifting, and skipping the attacks
    - start point is fixed at the playback beginning
        - what this means is that when accumulating the read point using the
          time-stretch amount, the start of the accumulation is given by
          sampling a position signal and successive positions are given by
          accumulating the time-stretch 
        - we also want start point to stick to an attack: this is the
          responsibility of the machinery that converts a normalized signal to a
          position in the sound file: in the case that we don't stick to
          attacks, this normalized signal is simply multiplied by the length of
          the file, in the case that we stick to attacks, this signal is rounded
          to the nearest atttack point or something

## Interface
- The way I'm thinking so that we don't have to worry about (needless) thread issues is something like this:

    init() {
        /* The object that does the time stretching / pitch shifting */
        tsps_obj = init_time_stretch_pitch_shift();
        /* The object that manages recording into buffers */
        recorder_obj = init_recorder();
        ...
    }

    dsp_tick() {
        /* Called either by the callback filling the DMA buffer (real-time) or
         * by a "systick" that looks when the input ringbuffers holding signal
         * values have values to compute with (pseudo-real-time, more latency, a
         * little safer to prevent underruns) */
        settings_adjustments = get_settings();
        /* Has some kind of routing scheme that routes the control changes to
         * the correct object */
        apply_settings(settings_adjustments);
        input_signal = get_input_signal();
        recorder_record(recorder_obj,input_signal);
        /* Get signals like pitch-shift amount, time-stretch rate, etc. */
        control_signals = get_control_signals();
        tsps_obj_process(tsps_obj,control_signals);
        /* Alternatively, callbacks could be registered with tsps_obj and called
         * when you call tsps_obj_process() */
        /* Because DMA will shift out data without intervention from the CPU,
         * tsp_obj_process could have the signature
         * tsp_obj_process(tsps_obj,control_signals,output_buffer) where
         * output_buffer is the next buffer that is registered with DMA */
    }


