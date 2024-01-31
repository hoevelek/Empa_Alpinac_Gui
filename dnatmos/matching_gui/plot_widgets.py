import sys
import matplotlib
from matplotlib.figure import Figure

from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

from PySide6.QtWidgets import QWidget, QVBoxLayout, QApplication



class MplCanvas(QWidget):
    """Create a canvas containing a Figure and a toolbar."""
    figure: Figure

    def __init__(self, figure: Figure,  parent: QWidget = None, tool_bar_down: bool = False) -> None:
        super().__init__(parent)
        
        layout = QVBoxLayout(self)
    
        self.figure = figure
        
        dynamic_canvas = FigureCanvas(figure)
        nav_toolbar = NavigationToolbar(dynamic_canvas, self)
        if tool_bar_down:
            layout.addWidget(dynamic_canvas)
            layout.addWidget(nav_toolbar)
        else:
            layout.addWidget(nav_toolbar)
            layout.addWidget(dynamic_canvas)

        self.setLayout(layout)

    @property
    def axes(self):
        """Return the axes of the figure on which the canvas is built."""
        if not hasattr(self, '_axes'):
            self._axes =  self.figure.subplots()
        return self._axes
    
    def draw(self):
        self.figure.canvas.draw()

if __name__ == "__main__":
    app = QApplication([])

    widget = MplCanvas(Figure(figsize=(5, 3)))
    widget.axes.plot([1,2,2,4])
    widget.show()

    sys.exit(app.exec())
