finder_charts
======

**A Python package to build finder charts using visible and near-infrared astronomy surveys**

You can create a clone of this package with the following command::

    git clone https://github.com/jgagneastro/finder_charts.git

One the package is downloaded, it is *imperative* to create a Python environment dedicated to its use, and install the exact versions of the packages that were used in its development. This is true because several Python packages will change over time in a non-retrocompatible way.

In order to create your own Python environment, open a terminal, navigate somewhere you will remember (e.g. your finder_charts directory), and enter the following command::

    python -m venv finder_charts_env

This will create a finder_charts_env directory where the exact versions of the Python packages required by finder_charts will be stored without interfering with your system installations. Once this is done, you need to activate this virtual environment with the follwing terminal command::

    source finder_charts_env/bin/activate

Once your virtual environment is activated, you should see that your command line now stars with the "(finder_charts_env) " flag before the usual prompts. Once you are located in this environment, you need to install all packages with the following command (here I am assuming you have navigated to the finder_charts directory)::

    pip install -r requirements.txt

Instead of navigating to your finder_charts installation directory, you can also alternatively just download the requirements file wherever your environment was created with something like wget (or manually download requirements.txt from this repo) and then install the contents::

    wget https://raw.githubusercontent.com/jgagneastro/finder_charts/master/requirements.txt
    pip install -r requirements.txt

Once the packages are installed, you should be able to launch Python and use finder_charts normally. Note every time you need to use finder_charts, you should launch the same mocapy Python environment, by navigating wherever you have created it, and launching the same command again::

    source finder_charts_env/bin/activate

Documentation
-------------

The full documentation for this project is not yet available.

The following Python command lines will allow you to create a finder chart at the position of SIMP J013656.5+093347.3 (RA=24.235917 an DEC=9.5631389, both in decimal degrees)::
    
    #Import the finder routine
    from jg_finderchart import finder

    #Create a finder chart and save it to disk
    finder('24.235917 9.5631389',plot=False,savepdf=True,buffer=True,keepfiles=True,skipdownloads=False)

Creating a finder chart is relatively slow because it requires downloading FITS image cutouts. Running this code should generate a PDF finder chart similar to the one included in this repository. Surveys such as PSO, VHS or UKIDSS that do not cover both full hemispheres will not always appear in a finder chart.

License
-------

Copyright 2022 Jonathan Gagne, adapted from previous work by Davy Kirkpatrick and Adam Schneider.

finder_charts is free software made available under the MIT License. For details see
the LICENSE file.
