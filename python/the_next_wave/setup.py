import os
from glob import glob
from setuptools import find_packages, setup

package_name = 'the_next_wave'

setup(
    name=package_name,
    version='0.1.0',
    packages=['the_next_wave'],
    data_files=[
        ('share/ament_index/resource_index/packages',
            ['resource/' + package_name]),
        ('share/' + package_name, ['package.xml']),
        (os.path.join('share', package_name, 'launch'), glob('launch/*.launch.py')),
        (os.path.join('share', package_name, 'config'), glob('config/*.yaml'))
    ],
    zip_safe=True,
    description='Deterministic ocean wave prediction from sparse buoy data. This methodology and codes generate phase-resolved reconstructions of ocean waves over short time-space scales using measurements of wave motion collected by sparse arrays of Surface Wave Instrument Floats w/ Tracking (SWIFTs)',
    tests_require=['pytest'],
    # NOTE: the following commented items are usually in setup.py for ROS 2 ament_python, but they
    #       are now located in pyproject.toml for compatibility with uv
    #install_requires=['setuptools'],
    #maintainer='anderson',
    #maintainer_email='anderson@mbari.org',
    #license='CC0-1.0',
    #entry_points={
    #    'console_scripts': [
    #        'the_next_wave_example = the_next_wave.example:main'
    #    ],
    #},
)
