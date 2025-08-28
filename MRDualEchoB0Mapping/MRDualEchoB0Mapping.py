import os
import unittest
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import SimpleITK as sitk
import sitkUtils
import numpy
import scipy
import math
import copy

#
# MRDualEchoB0Mapping
#

class MRDualEchoB0Mapping(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "MRDualEchoB0Mapping" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Quantification"]
    self.parent.dependencies = []
    self.parent.contributors = ["Junichi Tokuda (Brigham and Women's Hospital)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    MR dual-echo B0 field mapping for susceptibility correction and B0 field mapping.
    This module provides B0 field mapping capabilities using dual-echo sequences.
    """
    self.parent.acknowledgementText = """
    This module was developed based on a template created by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# MRDualEchoB0MappingWidget
#

class MRDualEchoB0MappingWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):

    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    # --------------------------------------------------
    # I/O Area
    # --------------------------------------------------
    ioCollapsibleButton = ctk.ctkCollapsibleButton()
    ioCollapsibleButton.text = "I/O"
    self.layout.addWidget(ioCollapsibleButton)

    ioFormLayout = qt.QVBoxLayout(ioCollapsibleButton)
    
    ioCommonFormLayout = qt.QFormLayout()
    ioFormLayout.addLayout(ioCommonFormLayout)

    # --------------------
    # Single Frame
    #
    singleFrameGroupBox = ctk.ctkCollapsibleGroupBox()
    singleFrameGroupBox.title = "Single Frame"
    singleFrameGroupBox.collapsed = False

    ioFormLayout.addWidget(singleFrameGroupBox)
    singleFrameFormLayout = qt.QFormLayout(singleFrameGroupBox)
    
    #
    # input volume selector
    #
    self.baselinePhaseSelector = slicer.qMRMLNodeComboBox()
    self.baselinePhaseSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.baselinePhaseSelector.selectNodeUponCreation = True
    self.baselinePhaseSelector.addEnabled = False
    self.baselinePhaseSelector.removeEnabled = False
    self.baselinePhaseSelector.noneEnabled = True
    self.baselinePhaseSelector.showHidden = False
    self.baselinePhaseSelector.showChildNodeTypes = False
    self.baselinePhaseSelector.setMRMLScene( slicer.mrmlScene )
    self.baselinePhaseSelector.setToolTip( "Pick the baseline phase map" )
    singleFrameFormLayout.addRow("Baseline Phase Volume: ", self.baselinePhaseSelector)

    #
    # input volume selector
    #
    self.referencePhaseSelector = slicer.qMRMLNodeComboBox()
    self.referencePhaseSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.referencePhaseSelector.selectNodeUponCreation = True
    self.referencePhaseSelector.addEnabled = False
    self.referencePhaseSelector.removeEnabled = False
    self.referencePhaseSelector.noneEnabled = True
    self.referencePhaseSelector.showHidden = False
    self.referencePhaseSelector.showChildNodeTypes = False
    self.referencePhaseSelector.setMRMLScene( slicer.mrmlScene )
    self.referencePhaseSelector.setToolTip( "Select a reference phase volume." )
    singleFrameFormLayout.addRow("Reference Phase Volume: ", self.referencePhaseSelector)

    #
    # B0 map volume selector
    #
    self.b0MapSelector = slicer.qMRMLNodeComboBox()
    self.b0MapSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.b0MapSelector.selectNodeUponCreation = True
    self.b0MapSelector.addEnabled = True
    self.b0MapSelector.removeEnabled = True
    self.b0MapSelector.noneEnabled = True
    self.b0MapSelector.renameEnabled = True
    self.b0MapSelector.showHidden = False
    self.b0MapSelector.showChildNodeTypes = False
    self.b0MapSelector.setMRMLScene( slicer.mrmlScene )
    self.b0MapSelector.setToolTip( "Select an output Delta B0 map." )
    singleFrameFormLayout.addRow("Output B0 Map (ppm): ", self.b0MapSelector)

    #
    # Apply Button
    #
    self.applyButtonSingle = qt.QPushButton("Generate a B0 Map")
    self.applyButtonSingle.toolTip = "Run the algorithm for the single frame input."
    self.applyButtonSingle.enabled = False
    singleFrameFormLayout.addRow(self.applyButtonSingle)

    # --------------------------------------------------
    # MR Parameters Area
    # --------------------------------------------------
    #
    mrParametersCollapsibleButton = ctk.ctkCollapsibleButton()
    mrParametersCollapsibleButton.text = "MR Parameters"
    self.layout.addWidget(mrParametersCollapsibleButton)

    mrParametersFormLayout = qt.QFormLayout(mrParametersCollapsibleButton)
    
    #
    # Gyromagnetic ratio (gamma)
    #
    self.gammaSpinBox = qt.QDoubleSpinBox()
    self.gammaSpinBox.objectName = 'gammaSpinBox'
    self.gammaSpinBox.setMaximum(100.0)
    self.gammaSpinBox.setMinimum(-100.0)
    self.gammaSpinBox.setDecimals(8)
    self.gammaSpinBox.setValue(42.576)
    self.gammaSpinBox.setToolTip("Gyromagnetic ratio / PI (rad / s T) (Default is proton )")
    mrParametersFormLayout.addRow("gamma (rad / s T): ", self.gammaSpinBox)

    #
    # Field strength (B0)
    #
    self.B0SpinBox = qt.QDoubleSpinBox()
    self.B0SpinBox.objectName = 'B0SpinBox'
    self.B0SpinBox.setMaximum(20.0)
    self.B0SpinBox.setMinimum(0.0)
    self.B0SpinBox.setDecimals(8)
    self.B0SpinBox.setValue(3.0)
    self.B0SpinBox.setToolTip("Static field strength (Tesla) (Used to compute ppm)")
    mrParametersFormLayout.addRow("B0 (Tesla): ", self.B0SpinBox)

    #
    # Echo time 1 (TE1)
    #
    self.TE1SpinBox = qt.QDoubleSpinBox()
    self.TE1SpinBox.objectName = 'TE1SpinBox'
    self.TE1SpinBox.setMaximum(10.0)
    self.TE1SpinBox.setMinimum(0.0)
    self.TE1SpinBox.setDecimals(8)
    self.TE1SpinBox.setValue(0.005)
    self.TE1SpinBox.setToolTip("Echo time (s)")
    mrParametersFormLayout.addRow("TE1 (s): ", self.TE1SpinBox)

    #
    # Echo time 2 (TE2)
    #
    self.TE2SpinBox = qt.QDoubleSpinBox()
    self.TE2SpinBox.objectName = 'TE2SpinBox'
    self.TE2SpinBox.setMaximum(10.0)
    self.TE2SpinBox.setMinimum(0.0)
    self.TE2SpinBox.setDecimals(8)
    self.TE2SpinBox.setValue(0.015)
    self.TE2SpinBox.setToolTip("Echo time (s)")
    mrParametersFormLayout.addRow("TE2 (s): ", self.TE2SpinBox)


    # --------------------------------------------------
    # Parameters Area
    # --------------------------------------------------
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # Interpolation mode
    #
    self.dispInterpFlagCheckBox = qt.QCheckBox()
    self.dispInterpFlagCheckBox.checked = 0
    self.dispInterpFlagCheckBox.setToolTip("If checked, voxels will be interpolated for display. ")
    parametersFormLayout.addRow("Display interpolation: ", self.dispInterpFlagCheckBox)

    #
    # Delta B0 color range
    #
    self.scaleRangeMaxSpinBox = qt.QDoubleSpinBox()
    self.scaleRangeMaxSpinBox.objectName = 'scaleRangeMaxSpinBox'
    self.scaleRangeMaxSpinBox.setMaximum(5000.0)
    self.scaleRangeMaxSpinBox.setMinimum(-5000.0)
    self.scaleRangeMaxSpinBox.setDecimals(4)
    self.scaleRangeMaxSpinBox.setValue(500.0)
    self.scaleRangeMaxSpinBox.setToolTip("Maximum value for the delta B0 map.")
    parametersFormLayout.addRow("Max. scale (ppm): ", self.scaleRangeMaxSpinBox)

    #
    # Lower threshold - We set threshold value to limit the range of intensity 
    #
    self.scaleRangeMinSpinBox = qt.QDoubleSpinBox()
    self.scaleRangeMinSpinBox.objectName = 'scaleRangeMinSpinBox'
    self.scaleRangeMinSpinBox.setMaximum(5000.0)
    self.scaleRangeMinSpinBox.setMinimum(-5000.0)
    self.scaleRangeMinSpinBox.setDecimals(4)
    self.scaleRangeMinSpinBox.setValue(-500.0)
    self.scaleRangeMinSpinBox.setToolTip("Minimum value for the delta B0 map.")
    parametersFormLayout.addRow("Min. scale (ppm): ", self.scaleRangeMinSpinBox)

    #
    # Check for automatic update
    #
    self.autoUpdateCheckBox = qt.QCheckBox()
    self.autoUpdateCheckBox.checked = 1
    self.autoUpdateCheckBox.setToolTip("Automatic Update: ")
    parametersFormLayout.addRow("Automatic Update", self.autoUpdateCheckBox)

    # connections
    
    self.applyButtonSingle.connect('clicked(bool)', self.onApplyButton)
    self.baselinePhaseSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.referencePhaseSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.b0MapSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    self.autoUpdateCheckBox.connect('toggled(bool)', self.onAutoUpdate)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

    self.tag = None


  def cleanup(self):
    pass


  def onSelect(self):
    self.applyButtonSingle.enabled = self.baselinePhaseSelector.currentNode() and self.referencePhaseSelector.currentNode() and self.b0MapSelector.currentNode()


  def onAutoUpdate(self):
    if self.autoUpdateCheckBox.checked == True:
      if self.referencePhaseSelector.currentNode():
        self.referencePhaseSelector.enabled = False
        refNode = self.referencePhaseSelector.currentNode()
        self.tag = refNode.AddObserver(vtk.vtkCommand.ModifiedEvent, self.onModelRefImageModifiedEvent)
      else: # Cannot set autoupdate 
        self.autoUpdateCheckBox.checked = False
    else:
      if self.tag:
        if self.referencePhaseSelector.currentNode():
          refNode = self.referencePhaseSelector.currentNode()
          refNode.RemoveObserver(self.tag)
        self.referencePhaseSelector.enabled = True

        
  def onModelRefImageModifiedEvent(self, caller, event):
    self.onApplyButton()

    
  def onApplyButton(self):
    logic = MRDualEchoB0MappingLogic()

    param = {}
    param['displayInterpolation']     = self.dispInterpFlagCheckBox.checked
    param['baselinePhaseVolumeNode']  = self.baselinePhaseSelector.currentNode()
    param['referencePhaseVolumeNode'] = self.referencePhaseSelector.currentNode()
    param['b0MapVolumeNode']          = self.b0MapSelector.currentNode()
    param['gamma']                    = self.gammaSpinBox.value
    param['B0']                       = self.B0SpinBox.value
    param['TE1']                      = self.TE1SpinBox.value
    param['TE2']                      = self.TE2SpinBox.value
    param['colorScaleMax']            = self.scaleRangeMaxSpinBox.value
    param['colorScaleMin']            = self.scaleRangeMinSpinBox.value

    logic.runSingleFrame(param)


  def onReload(self, moduleName="MRDualEchoB0Mapping"):
    # Generic reload method for any scripted module.
    # ModuleWizard will substitute correct default moduleName.

    globals()[moduleName] = slicer.util.reloadScriptedModule(moduleName)


#
# MRDualEchoB0MappingLogic
#

class MRDualEchoB0MappingLogic(ScriptedLoadableModuleLogic):

  def isValidInputOutputData(self, baselinePhaseVolumeNode, referencePhaseVolumeNode):
    """Validates if the output is not the same as input
    """
    if not baselinePhaseVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node for ParamA image defined')
      return False
    if not referencePhaseVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node for ParamB image defined')
      return False
    return True


  def generateDiskMask(self, refImage, center=[0.5,0.5,0.5], radius=0.5):

    refImageArray = sitk.GetArrayFromImage(refImage)
    dims = refImageArray.shape
    maxDim = numpy.double(numpy.max(dims))
    r = maxDim * radius
    
    cx = dims[0] * center[0]
    cy = dims[1] * center[1]
    cz = dims[2] * center[2]
    y, x, z = numpy.ogrid[-cx:dims[0]-cx, -cy:dims[1]-cy, -cz:dims[2]-cz]
    mask = x*x + y*y + z*z <= r*r
    voxels = numpy.zeros((dims[0], dims[1], dims[2]))
    voxels[mask] = 1
    
    maskImage = sitk.GetImageFromArray(voxels)
    maskImage.CopyInformation(refImage) 
    
    return maskImage

  
  def runSingleFrame(self, param): 
    """
    Run the actual algorithm
    """

    displayInterpolation     = param['displayInterpolation']
    baselinePhaseVolumeNode  = param['baselinePhaseVolumeNode']
    referencePhaseVolumeNode = param['referencePhaseVolumeNode']
    b0MapVolumeNode          = param['b0MapVolumeNode']
    gamma                    = param['gamma']
    B0                       = param['B0']
    TE1                      = param['TE1']
    TE2                      = param['TE2']
    colorScaleMax            = param['colorScaleMax']
    colorScaleMin            = param['colorScaleMin']


    if not self.isValidInputOutputData(baselinePhaseVolumeNode, referencePhaseVolumeNode):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    print(baselinePhaseVolumeNode.GetName())
    print(referencePhaseVolumeNode.GetName())

    imageBaseline = None
    imageReference = None

    # Check the scalar type (Siemens SRC sends image data in 'short' instead of 'unsigned short')
    baselineImageData = baselinePhaseVolumeNode.GetImageData()
    scalarType = ''
    if baselineImageData != None:
      scalarType = baselineImageData.GetScalarTypeAsString()

    imageBaseline  = sitk.Cast(sitkUtils.PullVolumeFromSlicer(baselinePhaseVolumeNode), sitk.sitkFloat64)
    imageReference  = sitk.Cast(sitkUtils.PullVolumeFromSlicer(referencePhaseVolumeNode), sitk.sitkFloat64)

    mask = None

    if b0MapVolumeNode:

      self.phaseDiff = None
      self.phaseDrift = None
      
      if scalarType == 'unsigned short':
        print('imageBaseline*numpy.pi/2048.0 - numpy.pi')
        imageBaselinePhase = imageBaseline*numpy.pi/2048.0 - numpy.pi
        imageReferencePhase = imageReference*numpy.pi/2048.0 - numpy.pi
      else:
        print('imageBaseline*numpy.pi/4096.0')
        imageBaselinePhase = imageBaseline*numpy.pi/4096.0
        imageReferencePhase = imageReference*numpy.pi/4096.0

      # Convert the phase images to complex images (as numpy arrays)
      arrayBaseline = sitk.GetArrayFromImage(imageBaselinePhase)
      arrayBaselineComplex = numpy.cos(arrayBaseline) + numpy.sin(arrayBaseline) * 1.0j
      arrayReference = sitk.GetArrayFromImage(imageReferencePhase)
      arrayReferenceComplex = numpy.cos(arrayReference) + numpy.sin(arrayReference) * 1.0j

      arrayPhaseDiffComplex = arrayReferenceComplex / arrayBaselineComplex
      arrayPhaseDiff = numpy.angle(arrayPhaseDiffComplex)

      # Calculate Delta B0
      deltaB0ppm = 1000000.0 * (arrayPhaseDiff / (gamma * (TE2-TE1)))/B0

      self.deltaB0ppm = sitk.GetImageFromArray(deltaB0ppm)
      self.deltaB0ppm.SetOrigin(imageBaseline.GetOrigin())
      self.deltaB0ppm.SetSpacing(imageBaseline.GetSpacing())
      self.deltaB0ppm.SetDirection(imageBaseline.GetDirection())

      print("(gamma, B0, TE1, TE2) = (%f, %f, %f, %f)" % (gamma, B0, TE1, TE2))

      sitkUtils.PushVolumeToSlicer(self.deltaB0ppm, b0MapVolumeNode.GetName(), 0, True)

      dnode = b0MapVolumeNode.GetDisplayNode()
      if dnode == None:
        dnode = slicer.mrmlScene.CreateNodeByClass('vtkMRMLScalarVolumeDisplayNode')
        slicer.mrmlScene.AddNode(dnode)
        b0MapVolumeNode.SetAndObserveDisplayNodeID(dnode.GetID())

      dnode.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileColdToHotRainbow.txt')
      dnode.SetAutoWindowLevel(0)
      dnode.SetWindowLevelMinMax(colorScaleMin, colorScaleMax)

      colorLegendDisplayNode = slicer.modules.colors.logic().AddDefaultColorLegendDisplayNode(b0MapVolumeNode)
      colorLegendDisplayNode.VisibilityOn()

      if displayInterpolation == True:
        dnode.SetInterpolate(1)
      else:
        dnode.SetInterpolate(0)

    logging.info('Processing completed')

    return True


  def ft3d(self, array):
    return scipy.fft.fftshift(scipy.fft.fftn(array))


  def ift3d(self, array):
    return scipy.fft.ifftn(scipy.fft.fftshift(array))


class MRDualEchoB0MappingTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_MRDualEchoB0Mapping1()

  def test_MRDualEchoB0Mapping1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    pass

