#ifndef PS_PROC_FILE_H
#define PS_PROC_FILE_H 

#include "pitch_shifter.h"
int
ps_proc_file(
    struct pitch_shifter *ps,
    const char *ps_file_path,
    const char *ts_file_path,
    const char *output_file_path);
#endif /* PS_PROC_FILE_H */
