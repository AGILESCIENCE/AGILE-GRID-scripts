#!/usr/bin/python

from setuptools import setup, find_packages

entry_points = {
	'console_scripts': [
		'analyzeAP          = ap.analyzeAP:main',
		'analyzeLoadAP4     = ap.analyzeLoadAP4:main',
		'analyzeResults     = ap.analyzeResults:main',
		'BBMainFunction     = ap.BBMainFunction:main',
		'display_lc         = ap.display_lc:main',
		'DisplayLC          = ap.DisplayLC:main',
		'drawLS             = ap.drawLS:main',
		'drawVM             = ap.drawVM:main',
		'Edp                = ap.Edp:main',
        'normalizeAP        = ap.normalizeAP:main',
		'runanalysisLoadAP4 = ap.runanalysisLoadAP4:main',
		'runcheckresults    = ap.runcheckresults:main',
		'runnormalizeAP     = ap.runnormalizeAP:main',
        'SimAP              = ap.SimAP:main'
    ]
}

setup( 
    name='ap',
    author='Bulgarelli Andrea, Parmiggiani Nicolò, Valentina Fioretti',
    author_email='andrea.bulgarelli@inaf.it, nicolo.parmiggiani@inaf.it, valentina.fioretti@inaf.it',
    packages=find_packages(),
    package_dir={ 'ap': 'ap' },
    entry_points=entry_points,
    include_package_data=True
)