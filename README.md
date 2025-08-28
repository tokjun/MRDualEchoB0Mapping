# MR B0 Mapping Module for 3D Slicer
[3D Slicer](https://slicer.org/) utility scripts for MR B0 Mapping


## Installation

Download the script as a zip or using the `git` command as:

~~~~
$ git clone https://github.com/tokjun/MRDualEchoB0Mapping
~~~~

Once the folder is deployed on the local computer, use the following steps to install the module in 3D Slicer (in case of Slicer 5.8.1)

1. `Edit` -> `Application Settings`
2. On the `Settings` dialog box, choose `Modules` from the left menu.
3. Click `>>` button on the right side of `Additional module paths:` list. `Add` and `Remove` buttons should appear.
4. Click `Add` button. On the file dialog box, choose `<path to the source code directory>/MRDualEchoB0Mapping/MRDualEchoB0Mapping`. (Make sure that MRDualEchoB0Mapping.py is under the selected directory.)
5. Restart 3D Slicer

If successful, the `MRDualEchoB0Mapping` should appear the `Quantification` section of the module menu.


## Usage

1. Load phase images for the first (TE1) and second (TE2) images.
2. Select the TE1 phase image from the `Baseline phase volume` pull-down menu.
3. Select the TE2 phase image from the `Reference phase volume` pull-down menu.
4. Click the `Output B0 Map` pull down menu, and choose `Create a new Volume` or `Create a new Volume as..` to create an output volume.
5. Set the B0 field strength under the `MR Parameters` section.
6. Set the echo times for TE1 and TE2 *in seconds* (not in milliseconds) under the `MR Parameters` section.
7. Click `Generate a B0 map`
8. To view the image, you may need to choose the output volume on the volume selector on each 2D viewer.

