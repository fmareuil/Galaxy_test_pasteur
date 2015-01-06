from ccpnmr.format.general.userInteraction import setupMultiDialog
from ccpnmr.format.general.Util import createSelection

def setupDialog(dialog, gui):


  if not dialog:
    if gui is None:
      try:
        import Tkinter
        gui = Tkinter.Tk()
      except:
        gui = None

    dialog = setupMultiDialog(gui)
  
  return dialog, gui


def setupDialogAndSelect(object, dialog, gui, **kw):

  (dialog, gui)   = setupDialog(dialog, gui)
  (selList, dict) = createSelection(object)

  kw['selectionDict'] = dict

  interaction = dialog.SelectionList(gui, selList, **kw)

  return interaction, dialog, gui


def selectConstraintList(ccpProject, gui, restraint_type=None, msg=None):

  from aria.importFromCcpn import getObjectKeyString

  if ccpProject is None:
    return None, 'No CCPN model selected.'
  
  nmrProject = ccpProject.currentNmrProject

  constraintStores = list(nmrProject.nmrConstraintStores)
  
  key = None

  ccpn_name = None
  
  if constraintStores:
  
    if len(constraintStores) > 1:
      selList, selDict = createSelection(constraintStores)
      constraintStore = gui_select_constraintStore(selList, selDict, gui=gui)
    else:
      constraintStore = constraintStores[0]
  
    if not constraintStore:
      #print 'No %s Constraint Set selected.' % msg
      m = 'No %s Constraint Set selected.' % msg	
      return key, m, ccpn_name
      
    constraintLists = list(constraintStore.findAllConstraintLists(className=restraint_type))
    
    if constraintLists:
      if len(constraintLists) > 1:
        selList, selDict = createSelection(constraintLists)
        constraintList = gui_select_constraintList(selList, selDict, gui=gui)
      else:
        constraintList = constraintLists[0]
    
      if not constraintList:
        #print 'No %s Constraint List selected.' % msg
	m = 'No %s Constraint List selected.' % msg
        return key, m, ccpn_name
     
      key = getObjectKeyString(constraintList)
     
    else:
      #print 'No %s Constraint Lists found in store %s.' % (msg, constraintStore.serial)
      m = 'No %s Constraint Lists found in store %s.' % (msg, constraintStore.serial)	
      return key, m, ccpn_name
      
  else:
    #print 'No %s Constraint Sets found.' % msg
    m = 'No %s Constraint Sets found.' % msg
    return key, m, ccpn_name
  
  m = 'Constraint List %s (%s) selected.' % (constraintList.name, msg)
  
  ccpn_name = "%d:%d:%s" % (constraintStore.serial, constraintList.serial, constraintList.name)
  
  return key, m, ccpn_name
  

def gui_select_molSystem(ccpProject, dialog=None, gui=None):
  """
  gui: tk parent widget
  """

  args = (ccpProject.molSystems, dialog, gui)
  
  kw = {'title': 'Select molecular system',
        'text' : 'Existing molecular system codes:',
        'urlFile': 'SelectMolSystem'}

  interaction, dialog, gui = setupDialogAndSelect(*args, **kw)

  if interaction.isSelected:
    molsystem = interaction.selection
  else:
    molsystem = None

  return molsystem


def gui_select_shiftList(selList, selDict, dialog=None, gui=None):
          
  dialog, gui = setupDialog(dialog, gui)

  select_list = dialog.SelectionList(gui, selList,
                     selectionDict = selDict,
                     title = 'Select shift list for this peak list',
                     text = 'Valid shift lists:',
                     dismissText = 'No valid shift list')

  if select_list.isSelected:
    shiftList = select_list.selection
  else:
    shiftList = None
    
  return shiftList

def gui_select_peakList(peaklist_labels, peaklist_dict, dialog=None,
             gui=None):

  dialog, gui = setupDialog(dialog, gui)

  select = dialog.SelectionList(gui, peaklist_labels,
                  selectionDict = peaklist_dict,
                  title = 'Select peak list for ARIA',
                  text = 'Remaining valid peak lists:',
                  dismissText = 'Finish selection')

  if select.isSelected:
    peakList = select.selection

  else:
    peakList = None

  return peakList

## BARDIAUX
def gui_select_model(model_labels, model_dict, dialog=None,
                     gui=None):

  dialog, gui = setupDialog(dialog, gui)

  select = dialog.SelectionList(gui, model_labels,
                  selectionDict = model_dict,
                  title = 'Select model for ARIA',
                  text = 'Remaining valid models:',
                  dismissText = 'Finish selection')

  if select.isSelected:
    model = select.selection

  else:
    model = None

  return model

def gui_select_ensemble(ensemble_labels, ensemble_dict, dialog=None,
                     gui=None):

  dialog, gui = setupDialog(dialog, gui)

  select = dialog.SelectionList(gui, ensemble_labels,
                  selectionDict = ensemble_dict,
                  title = 'Select Structure ensemble for ARIA',
                  text = 'Remaining valid Structure ensembles:',
                  dismissText = 'Finish selection')

  if select.isSelected:
    ensemble = select.selection

  else:
    ensemble = None

  return ensemble

## BARDIAUX CCPN restraints
def gui_select_constraintList(constraintlist_labels, constraintlist_dict, dialog=None,
                 gui=None):

  dialog, gui = setupDialog(dialog, gui)

  select = dialog.SelectionList(gui, constraintlist_labels,
                  selectionDict = constraintlist_dict,
                  title = 'Select Constraint list for ARIA',
                  text = 'Remaining valid constraint lists:',
                  dismissText = 'Finish selection')


  if select.isSelected:
    constraintList = select.selection

  else:
    constraintList = None

  return constraintList


def gui_select_constraintStore(labels, dict, dialog=None, gui=None):

  dialog, gui = setupDialog(dialog, gui)

  select = dialog.SelectionList(gui, labels,
                  selectionDict = dict,
                  title = 'Select Constraint Set for ARIA',
                  text = 'Available constraint sets:',
                  dismissText = 'Finish selection')

  if select.isSelected:
    constraintStore = select.selection
  else:
    constraintStore = None

  return constraintStore


def gui_select_structureGeneration(generation_labels, generation_dict, dialog=None,
                  gui=None):

  dialog, gui = setupDialog(dialog, gui)

  select = dialog.SelectionList(gui, generation_labels,
                  selectionDict = generation_dict,
                  title = 'Select Structure Generation for ARIA',
                  text = 'Remaining valid structure generation:',
                  dismissText = 'Finish selection')

  if select.isSelected:
    structureGeneration  = select.selection

  else:
    structureGeneration = None

  return structureGeneration

