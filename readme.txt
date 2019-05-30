Description
===========
multi2birculate is a Python script that takes a matrix of multistate data and returns the same data recoded as binary characters. It accepts matrices in the NEXUS format and can export the recoded data into the NEXUS, TNT, and Phylip formats. Characters can be ordered or unordered, and contain polymorphisms.

Installation
============
This script requires Python 3 and numpy.

You can install Python from https://www.python.org/. This script was developed with Python version 3.6. It will *not* work with Python 2.

After installing Python, you can proceed to install numpy. The easiest way of doing this is opening the Command Prompt (on Windows) or the Terminal (on Mac) and executing the following command:

```
    pip install numpy
```

On Windows, you may need to open the Command Prompt with administrator privileges. This can be done by right-clicking on the Command Prompt icon and choosing "Run as administrator" from the context menu. If you find other problems, you can also download an installer for numpy from https://pypi.python.org/pypi/numpy.

This script was developed and tested with numpy 1.13.3.

With Python and numpy installed, you only need to keep the script file in any place of your convenience. If you are not very familiar with command line interfaces, you may find it useful to keep the script file in the same folder as the matrices that you want to recode.

Usage
=====
The script is excecuted with commands on the Command Prompt. The following examples are going to assume that you are keeping the script and the matrices that you want to recode in the same folder, and that this folder is your current directory.

The simplest way of using the script is with the following command:

```
python multi2birculate.py example.nex
```

Where "example.nex" is the NEXUS file containing the multistate character matrix to recode. With that command, the script is going to print the reocded matrix in NEXUS format in the Command Prompt. In order to save the recoded matrix to a file, one can add the -o option followed by a name for the output file ("example_out.nex"), like this:

```
python multi2birculate.py example.nex -o example_out.nex
```

If the "example_out.nex" file did not exist, it will be created. Otherwise, it will be overwritten without warning.

Additional options of the script can be displayed with the following command:

```
python multi2birculate.py --help
```

Which prints out the following information:

```
optional arguments:
  -h, --help            show this help message and exit
  -o [OUTPUT], --output [OUTPUT]
                        Name for the output file
  -f {nexus,phylip,tnt}, --format {nexus,phylip,tnt}
                        Output format. If nexus is selected, detailed recoding
                        information will be written out as a comment
  -O, --ordered         Force all characters to be treated as ordered
  -U, --unordered       Force all characters to be treated as unordered
  -a, --remove_autapomorphic
                        Remove autapomorphic characters after recoding
  -i, --remove_invariant
                        Remove invariant characters after recoding
  -m, --remove_missing  Remove characters without coded data (all missing or
                        inapplicable) after recoding
```

All these arguments may be combined (except -O and -U together) to obtain the desired output. For instance, the following command would be appropriate to produce a binary matrix for analysis with RAxML:

```
python multi2birculate.py example.nex -f phylip -o example_out.phy -i -m
```

That command tells the script to save the matrix in Phylip format ("-f phylip"), and to delete characters from the recoded matrix that happen to be invariant ("-i") or contain only missing data ("-m").

Note that the script will follow the indications from the ASSUMPTIONS block of the input matrix to set characters as ordered or unordered. If a character is not included in the ASSUMPTIONS block (or there is no ASSUMPTIONS block), it will be treated as unordered. This behaviour can be overriden with the "-O" and "-U" options.

Input format
============
This script targets the NEXUS format as implemented by Mesquite, with some limitations (see below). Moreover the 'standard' data format that is commonly used for morphological data. It will not treat apropriately data that is in 'DNA', 'protein', or other formats.

There are a few things that you should know about how this script parses NEXUS files:

1) Species names ("taxa") can contain only the following symbols:
0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-.,/_'

2) In the matrix block, species names must be separated from the data by white space. Note that Mesquite some times omits spaces if the species name is between single quote marks.

3) The matrix must be in the 'sequential' format. It cannot be interleaved.

4) Only the following symbols can used to represent character states: 0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ

Inapplicables are always represented by -, and missing data by ?

5) The states of ordered characters must be represented by the symbols shown in (3), *in that exact order*.

6) If you have intermediate states for which there are no observations in your matrix (e.g. only observations for states 0 and 2), the unobserved intermediate state (1, in the previous example) will be included in the recoded matrix. Those intermediate states can be removed from the recoded matrix with the -i option.

7) You can specify which characters are ordered or unordered in a standard ASSUMPTIONS block. This script can understand hyphens to denote character ranges (e.g. 1-3 means characters 1, 2, and 3), but not the '/' notation used sometimes by Mesquite. If '/' is found in the ASSUMPTIONS block, it will be ignored.

8) There should be only one CHARACTERS or DATA block, and only one ASSUMPTIONS block in your file. TAXA blocks are not necessary and are not read by this script.