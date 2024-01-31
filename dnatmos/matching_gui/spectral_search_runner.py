# -*- coding: utf-8 -*-
"""
Created on Wed May 25 12:47:51 2022

@author: kaho
"""

import logging
from pathlib import Path
import sys

import numpy as np
from alpinac.mode_identification import make_identification
from PySide6.QtWidgets import QApplication, QFileDialog, QLineEdit, QLabel, QWidget, QVBoxLayout, QPushButton, QHBoxLayout, QFormLayout
from PySide6.QtCore import Signal
from alpinac_gui.matchms_funcs import score_matching
from alpinac_gui.widgets import FileSelector
import matchms as ms
import matplotlib


class DataBaseSearchStarter(QWidget):
    """
    
    """
    number_output_results: QLineEdit
    #file_selector: FileSelector
    start_button: QPushButton

    logger: logging.Logger

    spectral_search_started:Signal = Signal()
    spectral_search_comleted:Signal = Signal(list, list, list)

    query_spectrum: ms.Spectrum
    ref_spectra: list[ms.Spectrum]


    def __init__(self, query_spectrum: ms.Spectrum, ref_spectra: list[ms.Spectrum], parent: QWidget = None) -> None:
        super().__init__(parent)
        self.logger = logging.getLogger('alpinac_gui.DataBaseSearchStarter')

        self.query_spectrum = query_spectrum
        self.ref_spectra = ref_spectra
        
        layout = QVBoxLayout()

        self.number_output_results = QLineEdit('10', self)
        self.delta_mz_tolerance = QLineEdit('0.005', self)
        self.min_matching_peaks = QLineEdit('5', self)

        # TODO add int validator

        layout_edits = QFormLayout()
        layout_edits.addRow('Number of Best Results', self.number_output_results)
        layout_edits.addRow('mz Tolerance', self.delta_mz_tolerance)
        layout_edits.addRow('Min. # Matching Peaks', self.min_matching_peaks)

        layout.addLayout(layout_edits)

        self.start_button = QPushButton("Start spectral search", self)
        self.start_button.clicked.connect(self.start_matching_calculation)
        
        layout.addWidget(self.start_button)
    
        self.setLayout(layout)
    
    @property
    def number_of_returned_database_entries(self) -> int:
        return int(self.number_output_results.text())

    @property
    def delta_mz_tol_level(self) -> float:
        return float(self.delta_mz_tolerance.text())

    @property
    def min_number_of_matching_peaks(self) -> int:
        return int(self.min_matching_peaks.text())

    
    def start_matching_calculation(self):
        """Start spectral matching, with user-defined parameters. Atm: Number N to show N best spectral results
            Might be extended to define delta mz tolerance interval, matching algorithm, etc.
        """ 
        # Get the different variables you need

        # Try to run the function, or except log the error
        spect_list, scores, no_matching_peaks = score_matching(
            self.query_spectrum,
            self.ref_spectra,
            tolerance_level=self.delta_mz_tol_level,
            min_match_peaks=self.min_number_of_matching_peaks,
            no_best_matches=self.number_of_returned_database_entries,
        )

        # emit the signal that the search is done
        self.spectral_search_comleted.emit(spect_list, scores, no_matching_peaks)
        
    def start_mode_identification(self):
        # TODO: when you are done, remove this method
        """Start alpinac with the settings choosen by the user."""
        
        path = self.file_selector.data_path

        self.logger.info("Start identification")
        self.mode_identification_started.emit()
        try:
            make_identification(
                path_file = path,
                fig_path= path.with_suffix(''),
                formulas_path= path.with_suffix(''),
                target_elements= self.target_elements if self.target_elements else None,
                target_formula=self.target_formula if self.target_formula else None,
                show_plt_figures=False
            )
        except Exception as exp:
            self.logger.error("Could not make identification")
            raise exp
        self.mode_identification_stopped.emit()
        self.logger.info("End identification")

if __name__ == "__main__":
    app = QApplication([])

    fake_spectrum = ms.Spectrum(
        np.array([1, 2.1, 2.4]),
        np.array([0.3, 1, 0.5]),
        metadata={"test": "asdf", "test2": "asdf", "test3": None, "test4": 8.5},
    )
    other_spectrum = ms.Spectrum(
        np.array([1, 2.8, 3.4]),
        np.array([0.6, 1, 0.2]),
        metadata={"test": "asdf", "test2": "asdf", "test_lol": "asdfsd", "test4": 8.5},
    )

    widget = DataBaseSearchStarter(fake_spectrum, [fake_spectrum, other_spectrum])
    widget.logger.setLevel(logging.INFO)

    widget.show()

    sys.exit(app.exec())
