/*
Pitch tracker written in C for speed.
*/

#include "Python.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/halffloat.h"

#define wrap_phase(x) \
while ((x) >= M_PI) { (x) -= (2*M_PI); }\
while ((x) < -M_PI) { (x) += (2*M_PI); }

typedef complex double zdouble;

static PyMethodDef PtrackerMethods[] = {
        {NULL, NULL, 0, NULL}
};

struct zdouble_lpfilt {
    zdouble a;
    zdouble b0;
    zdouble yn_1;
};

static inline zdouble
zdouble_lpfilt_filt(
    struct zdouble_lpfilt *dlpf,
    zdouble x)
{
    zdouble yn = x*dlpf->b0 - dlpf->a * dlpf->yn_1;
    dlpf->yn_1 = yn;
    return yn;
}

/* The loop definition must precede the PyMODINIT_FUNC. */

static void zdouble_ptrackerprod(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];

    /* inputs */
    char *x = *args++,
         *w0 = *args++,
         *a_lp = *args++,
         *alph_p = *args++,
         *alph_w = *args++,
         *g_thresh = *args++,
    /* outputs */
         *a = *args++,
         *phi = *args++,
         *w = *args++;

    /* input steps */
    npy_intp x_step = *steps++,
             w0_step = *steps++,
             a_lp_step = *steps++,
             alph_p_step = *steps++,
             alph_w_step = *steps++,
             g_thresh_step = *steps++,
    /* output steps */
             a_step = *steps++,
             phi_step = *steps++,
             w_step = *steps++;

    double p = 0,    /* power */
           phi_g = 0,/* the phase of the signal we use to heterodyne */
           w_ = 0,   /* the estimated instantaneous frequency */
           w0_ = *((double*)w0); /* the initial guess frequency */

    zdouble h = 1,   /* current heterodyned sample */
            hf = 1,  /* current heterodyned, low-pass filtered sample */
            hf_ = 1, /* last heterodyned, low-pass filtered sample */
            g = 1;   /* the the signal we use to heterodyne */

    /* The low-pass filter for heterodyning */
    struct zdouble_lpfilt hlpfilt = {
        .a = -1 * *((double*)a_lp),
        .b0 = 1. - *((double*)a_lp),
        .yn_1 = 0
    };

    for (i = 0; i < n; i++) {
        h = *((double*)x) * g;
        hf = zdouble_lpfilt_filt(&hlpfilt,h);
        p = (1 - *((double*)alph_p)) * cabs(hf) + *((double*)alph_p)*p;
        if (p >= *((double*)g_thresh)) {
            w_ = carg(hf/hf_);
            wrap_phase(w_);
            w0_ += w_ * *((double*)alph_w);
            hf_ = hf;
        }
        phi_g += w0_;
        wrap_phase(phi_g);
        g = cexp(-1*I*phi_g);
        *((double*)w) = w0_;
        *((double*)phi) = phi_g;
        *((double*)a) = p;

        x += x_step; 
        w0 += w0_step; 
        a_lp += a_lp_step; 
        alph_p += alph_p_step; 
        alph_w += alph_w_step; 
        g_thresh += g_thresh_step; 
        a += a_step; 
        phi += phi_step; 
        w += w_step;
    }
}


/*This a pointer to the above function*/
PyUFuncGenericFunction funcs[1] = {&zdouble_ptrackerprod};

/* These are the input and return dtypes of ptracker.*/

static char types[9] = {
NPY_DOUBLE,
NPY_DOUBLE,
NPY_DOUBLE,
NPY_DOUBLE,
NPY_DOUBLE,
NPY_DOUBLE,
NPY_DOUBLE,
NPY_DOUBLE,
NPY_DOUBLE,
};


static void *data[1] = {NULL};

const char ptracker_docstring[] =
"ptracker\n"
"\n"
"input arguments:\n"
"x: \n"
"    the signal containing the sinusoid to track\n"
"w0: \n"
"    the initial guess, in radians/cycle\n"
"a_lp: \n"
"    the radius of the low-pass filter pole\n"
"alph_p: \n"
"    the coefficient smoothing the power estimation\n"
"alph_w: \n"
"    the coefficient smoothing the frequency tracking\n"
"g_thresh: \n"
"    the power below which no frequency updates are made (the estimate stays at\n"
"    the previous frequency)\n"
"\n"
"output arguments\n"
"a:\n"
"    the instantaneous amplitudes of the tracked sinusoid\n"
"phi:\n"
"    the instantaneous phases of the tracked sinusoid\n"
"w:\n"
"    the instantaneous frequencies of the tracked sinusoid\n";

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "ptrackers",
    NULL,
    -1,
    PtrackerMethods,
    NULL,
    NULL,
    NULL,
    NULL
};


PyMODINIT_FUNC PyInit_ptrackers(void)
{
    PyObject *m, *ptracker, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    ptracker = PyUFunc_FromFuncAndData(funcs, data, types, 1, 6, 3,
                                    PyUFunc_None, "ptracker",
                                    ptracker_docstring, 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "ptracker", ptracker);
    Py_DECREF(ptracker);

    return m;
}
#else
PyMODINIT_FUNC initptrackers(void)
{
    PyObject *m, *ptracker, *d;


    m = Py_InitModule("ptrackers", PtrackerMethods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    ptracker = PyUFunc_FromFuncAndData(funcs, data, types, 1, 6, 3,
                                    PyUFunc_None, "ptracker",
                                    ptracker_docstring, 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "ptracker", ptracker);
    Py_DECREF(ptracker);
}
#endif
