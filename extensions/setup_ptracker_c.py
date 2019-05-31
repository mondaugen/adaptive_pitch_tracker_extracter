def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('ptracking',
                           parent_package,
                           top_path)
    config.add_extension(
    'ptrackers',
    ['ptracker_c.c'],
    #extra_compile_args=['-O0']
    )

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)

