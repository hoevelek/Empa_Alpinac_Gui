"""Simple interface"""


import logging
import sys
import random
from typing import Optional
from PySide6 import QtCore, QtWidgets, QtGui
from PySide6.QtWidgets import QDialog, QLabel, QGroupBox, QGridLayout, QWidget, QVBoxLayout, QPushButton, QHBoxLayout, QFormLayout
from PySide6.QtCore import QSettings, Signal

from matchms import Spectrum
import matplotlib
from matplotlib.figure import Figure
from alpinac_gui.matchms_funcs import AlpinacData, draw_chemical_struct, score_matching
from alpinac_gui.spectrum_widgets import SpectrumDisplay, SpectrumMetadata, SpectrumSelector
from alpinac_gui.plot_widgets import MplCanvas
from alpinac_gui.widgets import FileSelector
from matchms_funcs import plot_mass_spectra
from matchms.importing import load_from_mgf
from pathlib import Path


import numpy as np

import matchms as ms





class SpectrumComparator(QtWidgets.QWidget):
    """Compare a reference spectrum with spectrums from a db.

    The reference spectrum and the db can be selected using the file selectors.
    """
    logger: logging.Logger

    analysed_spectrum_fileselector: FileSelector
    spectrum_db_fileselector: FileSelector

    spectrum_metadata: SpectrumMetadata
    spectrum_selector: SpectrumSelector

    # This is not used in this widget, but is included in the layout
    # because needed by the main alpinac widget
    alpinac_output_selector: SpectrumSelector

    canvas:  MplCanvas

    query_spectrum: ms.Spectrum
    selected_spectrums: list[ms.Spectrum]
    _db_spectrums: list[ms.Spectrum] # All the spectrum from the database

    def __init__(self, parent: QWidget = None) -> None:
        self.logger = logging.getLogger(f"alpinac_gui.{type(self).__name__}")
        super().__init__(parent)

        # Create the several widgets we need
        self.analysed_spectrum_fileselector = FileSelector(self, "analysed_spectrum")
        self.spectrum_db_fileselector = FileSelector(self, "spectrum_db")
        files_layout = QFormLayout()
        files_layout.addRow("Analysed Spectrum", self.analysed_spectrum_fileselector)
        files_layout.addRow("Reference Spectrum Database", self.spectrum_db_fileselector )

        self.spectrum_metadata = SpectrumMetadata(self, title="Analysed Spectrum")
        self.spectrum_selector = SpectrumSelector(self, title="DB search results")
        self.alpinac_output_selector = SpectrumSelector(self, title="Aplinac results")

        # spectrums is a Python list of matchms-type Spectrum objects
        # it store the currently selected spectrums
        self.selected_spectrums = []
        self.query_spectrum = ms.Spectrum(np.array([]), np.array([]))

        self.canvas = MplCanvas(Figure(), self, tool_bar_down=True)

        layout = QtWidgets.QGridLayout(self)
        layout.setColumnStretch(0, 1)
        layout.setColumnStretch(1, 0)



        layout.addLayout(files_layout, 0, 0)
        layout.addWidget(self.canvas, 1, 0, 2, 1)

        layout.addWidget(self.spectrum_metadata, 0, 1,)
        layout.addWidget(self.alpinac_output_selector, 1, 1, 1, 2)
        layout.addWidget(self.spectrum_selector, 2, 1, 1, 2)

        # When we change the selection, update all
        self.spectrum_selector.spectrum_changed.connect(
            self.spectrum_selection_changed
        )
        self.alpinac_output_selector.spectrum_changed.connect(
            self.spectrum_selection_changed
        )
        # When we want to show information on a spectrum, update the metadata
        self.spectrum_selector.show_spectrum_clicked.connect(
            self.open_spectrum_info
        )
        self.alpinac_output_selector.show_spectrum_clicked.connect(
            self.open_spectrum_info
        )

        # Changing the reference spectrum
        self.analysed_spectrum_fileselector.file_changed.connect(
            self.set_query_spectrum
        )
        # Change the database spectra
        self.spectrum_db_fileselector.file_changed.connect(
            self.load_spectrums
        )
        
        self.setLayout(layout)

        # 1) Press Buttom: Run Database search
        # 2) call function
        #    a) spectrums = matchms_funcs.score_matching(query_spectrum: ms.Spectrum, ref_spectra: ms.Spectrum | list[ms.Spectrum])
        #    b) self.spectrum_selector.set_spectra(spectrums)
        #        self.redraw_all()
        # 3) 
    def do_data_base_search(self, query_spectrum: ms.Spectrum, ref_spectra: ms.Spectrum | list[ms.Spectrum]):
        """ Should do the database search and handover spectra, so that they can be shown, ticked, etc. """
        #TODO
        spectrums = score_matching(query_spectrum, ref_spectra)
        self.spectrum_selector.set_spectra(spectrums)
        self.redraw_all()

    def set_query_spectrum(self, spectrum_path: Path):
        """Change the refenrence spectrum."""
        spectrum = AlpinacData.from_alpinac(spectrum_path)
        self.query_spectrum = spectrum
        self.spectrum_metadata.set_spectrum(spectrum)
        self.redraw_all()

    def open_spectrum_info(self, spectrum: ms.Spectrum):
        """Open the info of a spectrum in a new window."""
        dialog = QDialog()
        dialog.setLayout(QVBoxLayout())
        spectrum_metadata = SpectrumMetadata(self)
        spectrum_metadata.set_spectrum(spectrum)
        dialog.layout().addWidget(spectrum_metadata)
        # Add a canvas with the chemical structure
        canvas = MplCanvas(Figure(), self)
        draw_chemical_struct(canvas.axes, spectrum)
        dialog.layout().addWidget(canvas)
        dialog.exec()
    

    def load_spectrums(self, file_path: Path):
        """Load the spectrum form the given path."""
        # TODO: ADD the other possible formats
        # You could also make a standalone function for reading the file type
        # and add some processing in this method
        if file_path.suffix == '.mgf':
            self._db_spectrums = list(load_from_mgf(str(file_path)))
        else:
            raise FileNotFoundError((
                f"Unkown file format {file_path}.\n"
                "Accepted files are [.mgf, ...]" 
            ))
            
        for spectrum in self._db_spectrums:
            print('process ', spectrum)
        
        # TODO: Kaho: 
        # Define which specturms you want to show in the GUI and 
        # pass them to the spectrum_selector
        # Add the comparison at this point and show the best candidates.
        
        self.spectrum_selector.set_spectra(self._db_spectrums)
        self.redraw_all()

    def spectrum_selection_changed(self, spectrum: ms.Spectrum, need_to_show: bool, ):
        """Called when the spectrum selection is changed."""
        print(need_to_show, spectrum)

        if need_to_show:
            self.selected_spectrums.append(spectrum)
        else:
            self.selected_spectrums.remove(spectrum)
        # TODO imporve by replotting or removing only the spectra that was changed
        self.redraw_all()
    
    def redraw_all(self):
        """Redraw all the plots."""
        self.logger.debug(f"Redrawing {self.selected_spectrums}")
        self.canvas.axes.clear()
        plot_mass_spectra(
            self.canvas.axes,
            self.query_spectrum,
            self.selected_spectrums
        )

        # Update what is on the plot
        self.canvas.draw()


if __name__ == "__main__":
    app = QtWidgets.QApplication([])


    widget = SpectrumComparator()

    widget.show()

    sys.exit(app.exec())
