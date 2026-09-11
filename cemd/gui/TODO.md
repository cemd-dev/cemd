== Error at opening a datafile with the pypi version (not occuring in the dev version):
qt.accessibility.atspi: AtSpiAdaptor::applicationInterface does not implement "GetApplicationBusAddress" "/org/a11y/atspi/accessible/root"
Traceback (most recent call last):
  File "/home/jerome/miniconda3/envs/cemd/lib/python3.11/site-packages/cemd/gui/main_window.py", line 669, in open_file_clicked
    self.add_structure_tab(system, path)
  File "/home/jerome/miniconda3/envs/cemd/lib/python3.11/site-packages/cemd/gui/main_window.py", line 739, in add_structure_tab
    new_tab = StructureTabWidget(system, self, file_path=path)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/jerome/miniconda3/envs/cemd/lib/python3.11/site-packages/cemd/gui/tabs.py", line 49, in __init__
    self.setup_ui()
  File "/home/jerome/miniconda3/envs/cemd/lib/python3.11/site-packages/cemd/gui/tabs.py", line 57, in setup_ui
    self.plotter = AtomicPlotter(self, config=self.parent_gui.global_config)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/jerome/miniconda3/envs/cemd/lib/python3.11/site-packages/cemd/gui/plotter_widget.py", line 65, in __init__
    self.add_axes()
  File "/home/jerome/miniconda3/envs/cemd/lib/python3.11/site-packages/cemd/gui/plotter_widget.py", line 131, in add_axes
    self.axes_widget.SetInteractor(self.interactor)
TypeError: SetInteractor argument 1: method requires a VTK object

=> Versioning issue => solved