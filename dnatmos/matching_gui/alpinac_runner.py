import logging
from pathlib import Path
import sys
from alpinac.mode_identification import make_identification
from PySide6.QtWidgets import QApplication, QFileDialog, QLineEdit, QLabel, QWidget, QVBoxLayout, QPushButton, QHBoxLayout, QFormLayout
from PySide6.QtCore import Signal
from alpinac_gui.widgets import FileSelector



class AlpinacStarter(QWidget):
    """A widget that starts alpinac.
    
    It emits signals when the mode identification is started and stopped.

    Signals:
        mode_identification_started: The mode identification has been
            started by the user.
        mode_identification_stopped: The mode identification has been 
            finished.
    
    TODO: When we split the mode identification in different steps,
        we can show the steps in a progress dialog.
    """
    target_elements_lineedit: QLineEdit
    target_formula_lineedit: QLineEdit
    file_selector: FileSelector
    start_button: QPushButton

    logger: logging.Logger

    mode_identification_started:Signal = Signal()
    mode_identification_stopped:Signal = Signal()


    def __init__(self, parent: QWidget = None, file_selector: FileSelector = None) -> None:
        super().__init__(parent)
        self.logger = logging.getLogger('alpinac_gui.AlpinacStarter')
        
        layout = QVBoxLayout()

        if file_selector is None:
            file_selector = FileSelector(self, "AlpinacLauncher")
            layout.addWidget(file_selector)
        self.file_selector = file_selector

        self.target_elements_lineedit = QLineEdit(self)
        self.target_formula_lineedit = QLineEdit(self)

        layout_edits = QFormLayout()
        layout_edits.addRow('Target Elements', self.target_elements_lineedit )
        layout_edits.addRow('Example: ', QLabel("CHONCrB") )
        layout_edits.addRow('Target Formula', self.target_formula_lineedit )
        layout_edits.addRow('Example: ', QLabel("C6H12O8") )

        layout.addLayout(layout_edits)

        self.start_button = QPushButton("Start alpinac", self)
        self.start_button.clicked.connect(self.start_mode_identification)
        
        layout.addWidget(self.start_button)
    
        self.setLayout(layout)
    
    @property
    def target_elements(self) -> str:
        return self.target_elements_lineedit.text()

    @property
    def target_formula(self) -> str:
        return self.target_formula_lineedit.text()

    def get_input_frag_file_name(self) -> Path:
        #input_frag_file = e1.get()
        filetypes = 'text files (*.txt);; All files (*.*)'

        input_frag_file_name, _filter = QFileDialog.getOpenFileName(
            self,
            caption='Open a Fragment file',
            #initialdir='/',
            dir="C:/Users/mygu/halosearch/data/nontarget_screening/fragments",
            filter=filetypes
        )

        #input_frag_file_name = input_frag_file.name
        self.logger.info(input_frag_file_name)
        return Path(input_frag_file_name)
    
    def start_mode_identification(self):
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

    widget = AlpinacStarter()
    widget.logger.setLevel(logging.INFO)

    widget.show()

    sys.exit(app.exec())
