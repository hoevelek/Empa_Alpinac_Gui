import logging
import sys
from pathlib import Path

from PySide6.QtCore import QSettings, Signal
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QWidget,
)


class FileSelector(QWidget):
    """Choose which file is currently being used.
    
    Contain a LineEdit where one can change the text of the path
    and a button for browsing in the files of the computer.

    Store in the settings which was the last selected file.
    """

    file_changed: Signal = Signal(Path)
    data_path: Path

    def __init__(
        self,
        parent=None,
        unique_save_name: str = ""
    ) -> None:
        """Initialize directory selector Widget.
        
        :arg unique_file_name: Used for saving the file in a specific path.
        """
        self.logger = logging.getLogger("alpinac_gui.FileSelector")
        super().__init__(parent)

        self.settings = QSettings("alpinac", unique_save_name, parent=self)

        layout = QHBoxLayout()
        # Set initial file location
        self.data_path = self.get_initial_path()
        # Create a line edit where we can change the file name
        file_edit = QLineEdit(str(self.data_path), self)
        self.file_edit = file_edit

        browse_button = QPushButton("Browse...", self)

        file_edit.editingFinished.connect(self.on_editing_finished)
        browse_button.clicked.connect(self.select_file_dialog)

        # Adds the widgets to the layout
        layout.addWidget(file_edit)
        layout.addWidget(browse_button)

        self.setLayout(layout)

    def select_file_dialog(self):
        """Open a dialog for the user to select the file."""

        # Note a bug that makes it impossible to view the files
        # https://stackoverflow.com/questions/26078139/qfiledialog-view-files-and-folders-when-choosing-a-directory
        new_path, _ = QFileDialog.getOpenFileName(
            parent=self,
            caption="Select Spectrum File",
            dir=str(self.data_path.parent),
            # options=QFileDialog.Option.DontUseNativeDialog,
        )
        # print(f"{new_path=}")
        if new_path == "":
            # The folder selection was cancelled by user
            return
        new_path = Path(new_path)
        if new_path == self.data_path:
            # The file did not change
            return

        self.file_edit.setText(str(new_path))

        self.save_path(new_path)
        self.emit_filechanged_signal()

    def emit_filechanged_signal(self):
        self.logger.debug("Emitting file_changed signal.")
        self.file_changed.emit(self.data_path)

    def on_editing_finished(self):
        """Called when the editing of the line edit is changed."""
        new_file = self.file_edit.text()
        try:
            new_filepath = Path(new_file)
        except Exception as err:
            self.logger.error(f"{new_file} is not a valid filepath.")
            return

        # The path is valid
        if new_filepath.is_file():
            # Change the file
            self.save_path(new_filepath)
            self.emit_filechanged_signal()
        else:
            self.logger.error(f"No file : '{new_file}'.")

    def save_path(self, path: Path):
        """Save the new path in the settings."""
        self.data_path = path
        self.settings.setValue("SpectrumPath", str(path))

    def get_initial_path(self) -> Path:
        """Return what should be the intial path when launching the app."""
        return Path(self.settings.value("SpectrumPath", "."))


if __name__ == "__main__":
    app = QApplication([])

    widget = FileSelector()
    widget.file_changed.connect(lambda x: print("file changed", x))
    widget.show()

    sys.exit(app.exec())
