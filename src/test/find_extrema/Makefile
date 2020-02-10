local_max_f32_test_deps = src/find_extrema.c \
                          src/find_extrema.h \
                          src/test/test_common.h

define local_max_f32_test_template =
mkdir -p src/test/bin/find_extrema
cc -DLOCAL_MAX_TYPE=$(1) -g -Isrc -Isrc/test $(filter %.c, $^) -o $@
endef

local_max_f32_test_types = both none left right

local_max_f32_test_target_stem = \
src/test/bin/find_extrema/local_max_f32_test_type_

run_local_max_f32_test_target_stem = \
run_local_max_f32_test_type_

# TODO: Make tests for other local max types
$(foreach type,$(local_max_f32_test_types),\
$(run_local_max_f32_test_target_stem)$(type)) : \
    $(run_local_max_f32_test_target_stem)% : \
    src/test/bin/find_extrema/local_max_f32_test_type_% \
    src/test/find_extrema/local_max_f32_type_%_gen_input.py \
    src/test/find_extrema/local_max_f32_check_output.py
	PYTHONPATH=. python3 \
    src/test/find_extrema/local_max_f32_type_$*_gen_input.py
	$(local_max_f32_test_target_stem)$*
	PYTHONPATH=. python3 \
	src/test/find_extrema/local_max_f32_check_output.py

$(foreach type,$(local_max_f32_test_types),\
$(local_max_f32_test_target_stem)$(type)) : \
    $(local_max_f32_test_target_stem)% : \
    src/test/find_extrema/local_max_f32_test.c \
    $(local_max_f32_test_deps)
	$(call local_max_f32_test_template,local_max_type_$*)

