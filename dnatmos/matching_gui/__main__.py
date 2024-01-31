"""The main widget for the GUI of alpinac"""

import logging
import sys
from PySide6.QtWidgets import QApplication, QWidget, QDialog, QPushButton, QVBoxLayout, QErrorMessage
from alpinac_gui.alpinac_runner import AlpinacStarter
from alpinac_gui.matchms_funcs import AlpinacData

from alpinac_gui.pyside_alpinaq import SpectrumComparator
from alpinac_gui.spectral_search_runner import DataBaseSearchStarter
import matchms as ms


class MainWidget(SpectrumComparator):
    """The main widget of alpinac"""

    alpinac_starter: AlpinacStarter
    alpinac_starter_dialog: QDialog

    def __init__(self, parent: QWidget = None) -> None:
        super().__init__(parent)

        # Define widgets for alpinac dialog
        self.alpinac_starter = AlpinacStarter()
        self.alpinac_starter_dialog = QDialog(self)
        dialog_layout = QVBoxLayout()
        dialog_layout.addWidget(self.alpinac_starter)
        self.alpinac_starter_dialog.setLayout(dialog_layout)
        # Close the dialog when alpinac finishes
        self.alpinac_starter.mode_identification_stopped.connect(
            self.on_mode_identification_finished
        )
        # TODO: connect the mode_identification_stopped event to other function
        # so that it updates some information in the comparator when the
        # widget is terminated

        button_layout = QVBoxLayout()

        # Add a button to open the dialog
        start_alpinac_button = QPushButton("Start Alpinac", self)
        start_alpinac_button.clicked.connect(self.open_alpinac_dialog)

        # Add a button to open the dialog
        search_button = QPushButton("Start Search", self)
        search_button.clicked.connect(self.open_search_dialog)



        # Simply expend the grid layout from the spectrum comparator to fit it in
        button_layout.addWidget(start_alpinac_button)
        button_layout.addWidget(search_button)

        self.layout().addLayout(button_layout, 0, 2)

    def open_alpinac_dialog(self):
        """Open a window with the AlpinacStarter."""

        self.alpinac_starter_dialog.open()
        # TODO: set additional information on how the starter should be populated here
    

    def open_search_dialog(self):
        """Open a search dialog."""
        if len(self.query_spectrum.peaks.mz) == 0:
            # Empty spectrum
            mess = QErrorMessage(self)
            mess.showMessage(
                "The selected spectrum is empty.\n"
                "Please select a valid one for the search."
            )
            return 
        if not hasattr(self, '_db_spectrums') or not self._db_spectrums:
            # Empty spectrum
            mess = QErrorMessage(self)
            mess.showMessage(
                "The database is not loaded.\n"
                "Please select a database file."
            )
            return 
        search_dialog = QDialog(self)
        dialog_layout = QVBoxLayout()
        search_widget = DataBaseSearchStarter(self.query_spectrum, self._db_spectrums, search_dialog)
        dialog_layout.addWidget(search_widget)
        search_dialog.setLayout(dialog_layout)

        # Connect completed serach to callback
        search_widget.spectral_search_comleted.connect(self.on_spectral_search_finished)

        search_dialog.open()
    
    def on_mode_identification_finished(self):
        """Called when the mode identification is finished"""

        # Get the file
        file = self.alpinac_starter.file_selector.data_path
        self.logger.debug(f"on_mode_identification_finished: Alpinac file is {file}")
        output_dir = file.with_suffix('')

        spectrums = [
            # Find all the alpinac file in this output dir
            AlpinacData.from_alpinac(f)
            for f in output_dir.rglob('*.txt')
        ]

        self.logger.debug(f"on_mode_identification_finished: found {spectrums=}")
            

        # Load the alpinac output
        self.alpinac_output_selector.set_spectra(
            spectrums,
            additional_columns={
                # TODO make the columns you wnat
                'column you want': ['example' for _ in spectrums],
                's': [spec for spec in spectrums],
            }
        )

        self.alpinac_starter_dialog.close()


    def on_spectral_search_finished(self, spectrums: list[ms.Spectrum], scores: list, no_matching_peaks: list):
        """Called when the search is finished."""
        self.logger.setLevel(logging.INFO)
        self.logger.info(
            f"Finisihed search, found\n"
            + '\n'.join([
                f"{spec}, {score}, {nmp}" for spec, score, nmp in zip(spectrums, scores, no_matching_peaks)
            ])
        )
        # update the parts of the widgets which should be updated
        self.spectrum_selector.set_spectra(spectrums, scores=scores, no_matching_peaks=no_matching_peaks)
        self.selected_spectrums = []
        self.redraw_all()

if __name__ == "__main__":
    
    app = QApplication([])

    widget = MainWidget()

    widget.show()

    sys.exit(app.exec())

