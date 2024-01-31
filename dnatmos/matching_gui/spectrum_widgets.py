"""Contains some QWidgets for showing spectra."""

import sys
from typing import AnyStr
from PySide6.QtWidgets import (
    QApplication,
    QLabel,
    QWidget,
    QPushButton,
    QHBoxLayout,
    QFormLayout,
    QVBoxLayout,
    QCheckBox,
)
from PySide6.QtCore import QSettings, Signal

from matplotlib.figure import Figure
import numpy as np

from alpinac_gui.plot_widgets import MplCanvas

import matchms as ms

MAX_METADATA_STR_LEN = 50

class SpectrumSelector(QWidget):
    """A tool for selecting the Specturms shown on the plot.

    Signals:
        spectrum_changed: When the state of a parameter is changed.
            Returns two args.
            1. The name of the parameter (str)
            2. Whether now it is checked or not (bool)

    """

    spectrum_changed: Signal = Signal(ms.Spectrum, bool)
    show_spectrum_clicked: Signal = Signal(ms.Spectrum)

    # Stores the widgets that are currently shown
    _widgets: list[QWidget]

    title: str

    def __init__(self, parent: QWidget = None, title: str = "") -> None:
        super().__init__(parent)

        # The layout of this widget 
        all_layout = QVBoxLayout(self)

        self.setLayout(all_layout)
        self._widgets = []

        self.title = title

    def set_spectra(
        self,
        spectra: list[ms.Spectrum],
        empty: bool = True,
        scores: list[float] = [],
        no_matching_peaks: list[float] = [], 
        additional_columns: dict[str, list[AnyStr]] = {}
    ):
        """Set the spectra that should be shown inside the box.
        
        :arg spectra: The list of spectrum to be displayed.
        :arg empty: Whether the spectra already there should be removed.
        :arg additional_columns: Additional column to show inside the thing.
        """

        # TODO: add test spectra and scores should have same length

        if empty:
            # Delete all the widgets and remove all layout elements
            while item := self.layout().itemAt(0):
                self.layout().removeItem(item)
            
            while self._widgets:
                w = self._widgets.pop(0)
                w.hide()
                w.destroy()
        
        if not self._widgets and self.title:
            # Add the title on top
            layout = QHBoxLayout(self)
            title_label = QLabel(self.title, self)
            layout.addWidget(title_label)
            self.layout().addLayout(layout)
            self._widgets.append(title_label)

        for i, spectrum in enumerate(spectra):
            # Each parameter gets a line with different widgets
            layout = QHBoxLayout(self)
            name = self._name_to_show(spectrum)
            if scores:
                if len(scores) != len(spectra):
                    print("error: spectra list & score list have not the same length")
                else:
                    # add the score to the name
                    name = f"{name} Score: {round(scores[i],4)}"

            if no_matching_peaks:
                if len(scores) != len(no_matching_peaks):
                    print("error: spectra list and #matching peaks list have not the same length")
                else:
                    # add the score to the name
                    name = f"{name} # Matching Peaks: {no_matching_peaks[i]}"

            checkbox = QCheckBox(name, self)
            show_button = QPushButton('Show', self)

            layout.addWidget(checkbox)

            # Finds all the additional column that we need to add
            for name, value_list in additional_columns.items():
                # Create the label to show the value
                label = QLabel(f"{name} : {value_list[i]}", self)
                # Adds it to layout
                layout.addWidget(label)
                self._widgets.append(label)

            layout.addWidget(show_button)

            self.layout().addLayout(layout)

            def open_spectrum(spec=spectrum):
                self.show_spectrum_clicked.emit(spec)

            def state_changed(state: int, spec=spectrum):
                """Called when the button is clicked."""
                self.spectrum_changed.emit(spec, state)

            show_button.pressed.connect(open_spectrum)
            checkbox.stateChanged.connect(state_changed)

            self._widgets.append(checkbox)
            self._widgets.append(show_button)

    
    def _name_to_show(self, spectrum: ms.Spectrum) -> str:
        """Return the name that should be shown in the selector."""

        # Find out which name should be used
        name = spectrum.metadata.get('name', '')
        # If name was not found, look for other names
        name = name or spectrum.metadata.get('compound_name', '')

        # Do the same for formula
        formula = spectrum.metadata.get('formula', '')
        formula = formula or spectrum.metadata.get('precursor_formula', '')


        if name and formula:
            # If both info, show a combined version
            return f"{name} ({formula})"
        elif name:
            return name 
        elif formula:
            return formula
        else:
            # We could not extract the information
            return "Spectra with unkown name/formula."

class SpectrumMetadata(QWidget):
    """Show the metadata of a spectrum.
    
    """
    informations_layout: QFormLayout
    informations_widgets: dict[str, QLabel]
    def __init__(self, parent: QWidget = None, title: str = "Metadata") -> None:
        super().__init__(parent)
        
        self.informations_layout = QFormLayout()
        self.informations_widgets = {}

        self.informations_layout.addRow(QLabel(title))

    
        self.setLayout(self.informations_layout)


    def set_spectrum(self, spectrum: ms.Spectrum) -> None:
        """Add the metadata to the informations layout"""
        # Will store all the metadata which is not updated
        remaining_metadata = list(self.informations_widgets.keys())
        # Access all the keys and values form the metadata
        for key, value in spectrum.metadata.items():
            # Convert the value in the metadata to a string
            value = str(value)
            # Cropp the string if too long
            if len(value) > MAX_METADATA_STR_LEN:
                value = value[:MAX_METADATA_STR_LEN] + ' ... ' + value[-3:]

            if key in self.informations_widgets:
                # Just update the existing label if the key is already there
                self.informations_widgets[key].setText(value)
            else:
                # Create a new label and add it to the widget
                label = QLabel(value, self)
                self.informations_widgets[key] = label
                # Add it to the layout
                self.informations_layout.addRow(key, label)

            if key in remaining_metadata:
                # remove the key as it was updated
                remaining_metadata.remove(key)

        for key in remaining_metadata:
            self.informations_widgets[key].setText("")

class SpectrumDisplay(QWidget):
    """Display the spectrum.

    The plot is shown on the canvas and the metadata in the informations_layout.
    """

    canvas: MplCanvas


    def __init__(self, parent: QWidget = None) -> None:
        super().__init__(parent)

        layout = QHBoxLayout()
        self.canvas = MplCanvas(Figure(), self)

        layout.addWidget(self.canvas)

        self.setLayout(layout)

    def set_spectrum(self, spectrum: ms.Spectrum) -> None:
        """Set the spectrum to be displayed."""
        self.canvas.axes.clear()
        # TODO add the data to the plot
        self.canvas.axes.stem(spectrum.peaks.mz, spectrum.peaks.intensities)
        # Update what is on the plot
        self.canvas.figure.canvas.draw()




if __name__ == "__main__":
    # Small example how the widget looks like and how we can change the spectrum
    app = QApplication([])

    spec_display = SpectrumDisplay()

    spec_display.show()

    metadata = SpectrumMetadata()
    metadata.show()

    button_change = QPushButton("CHANGE THE SPECTURM LIST")
    button_change.show()

    fake_spectrum = ms.Spectrum(
        np.array([1, 2.1, 2.4]),
        np.array([0.3, 1, 0.5]),
        metadata={"name": "fake", "test": "asdf", "test2": "asdf", "test3": None, "test4": 8.5},
    )
    other_spectrum = ms.Spectrum(
        np.array([1, 2.8, 3.4]),
        np.array([0.6, 1, 0.2]),
        metadata={"test": "asdf", "test2": "asdf", "test_lol": "asdfsd", "test4": 8.5},
    )



    selector = SpectrumSelector(title="TESTSF SDFA")
    selector.set_spectra([fake_spectrum, other_spectrum])
    selector.set_spectra(
        [fake_spectrum, other_spectrum],
        empty=False, 
        additional_columns={"example":[1,2], "otherthing":['a', 93.2]}
    )
    selector.show()

    selector.show_spectrum_clicked.connect(metadata.set_spectrum)
    selector.show_spectrum_clicked.connect(spec_display.set_spectrum)

    button_change.clicked.connect(lambda : selector.set_spectra([other_spectrum, fake_spectrum, fake_spectrum]))

    sys.exit(app.exec())
