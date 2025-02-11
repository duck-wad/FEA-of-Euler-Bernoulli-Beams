import sys
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                            QHBoxLayout, QPushButton, QGraphicsView, QGraphicsScene,
                            QGraphicsEllipseItem, QGraphicsLineItem)
from PyQt6.QtCore import Qt, QPointF, QLineF, QRectF
from PyQt6.QtGui import QPainter, QPen, QColor

GRID_STEP = 50

class NodeItem(QGraphicsEllipseItem):
    def __init__(self, x, y):
        super().__init__(x-5, y-5, 10, 10)
        self.setBrush(QColor(0, 0, 255))
        self.setPen(QPen(Qt.GlobalColor.black, 1))
        self.pos = QPointF(x, y)

class BeamItem(QGraphicsLineItem):
    def __init__(self, start_node, end_node):
        line = QLineF(start_node.pos, end_node.pos)
        super().__init__(line)
        self.setPen(QPen(Qt.GlobalColor.darkGreen, 3))
        self.start_node = start_node
        self.end_node = end_node

class Canvas(QGraphicsView):
    def __init__(self):
        super().__init__()
        self.scene = QGraphicsScene()
        self.setScene(self.scene)
        self.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.nodes = []
        self.beams = []
        
        self.modes = {
            'add_node': False,
            'draw_beam': False,
            'delete': False
        }
        self.temp_beam_start = None
        self.draw_grid()

    def draw_grid(self):
        pen = QPen(QColor(200, 200, 200), 1, Qt.PenStyle.DotLine)
        rect = QRectF(-5000, -5000, 10000, 10000)
        x = -5000
        while x <= 5000:
            self.scene.addLine(x, rect.top(), x, rect.bottom(), pen)
            x += GRID_STEP
        y = -5000
        while y <= 5000:
            self.scene.addLine(rect.left(), y, rect.right(), y, pen)
            y += GRID_STEP

    def delete_node(self, node):
        # Remove connected beams
        for beam in self.beams.copy():
            if beam.start_node == node or beam.end_node == node:
                self.scene.removeItem(beam)
                self.beams.remove(beam)
        # Remove node
        self.scene.removeItem(node)
        self.nodes.remove(node)

    def delete_beam(self, beam):
        self.scene.removeItem(beam)
        self.beams.remove(beam)

    def snap_to_grid(self, pos):
        x = round(pos.x() / GRID_STEP) * GRID_STEP
        y = round(pos.y() / GRID_STEP) * GRID_STEP
        return QPointF(x, y)

    def node_exists_at(self, pos):
        for node in self.nodes:
            if node.pos == pos:
                return True
        return False

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Beam FEA")
        self.setGeometry(100, 100, 800, 600)
        
        # Create widgets
        self.canvas = Canvas()
        btn_add_node = QPushButton("Add Node")
        btn_draw_beam = QPushButton("Draw Beam")
        btn_delete = QPushButton("Delete Element")
        
        # Setup layout
        control_layout = QHBoxLayout()
        control_layout.addWidget(btn_add_node)
        control_layout.addWidget(btn_draw_beam)
        control_layout.addWidget(btn_delete)
        
        main_layout = QVBoxLayout()
        main_layout.addLayout(control_layout)
        main_layout.addWidget(self.canvas)
        
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)
        
        # Connect signals
        btn_add_node.clicked.connect(self.activate_add_node_mode)
        btn_draw_beam.clicked.connect(self.activate_draw_beam_mode)
        btn_delete.clicked.connect(self.activate_delete_mode)
        self.canvas.mousePressEvent = self.canvas_mouse_press

    def activate_add_node_mode(self):
        self.canvas.modes = {'add_node': True, 'draw_beam': False, 'delete': False}

    def activate_draw_beam_mode(self):
        self.canvas.modes = {'add_node': False, 'draw_beam': True, 'delete': False}
        self.canvas.temp_beam_start = None

    def activate_delete_mode(self):
        self.canvas.modes = {'add_node': False, 'draw_beam': False, 'delete': True}

    def canvas_mouse_press(self, event):
        pos = self.canvas.mapToScene(event.pos())
        
        if self.canvas.modes['add_node']:
            grid_pos = self.canvas.snap_to_grid(pos)
            if not self.canvas.node_exists_at(grid_pos):
                node = NodeItem(grid_pos.x(), grid_pos.y())
                self.canvas.scene.addItem(node)
                self.canvas.nodes.append(node)
                
        elif self.canvas.modes['draw_beam']:
            snapped_node = None
            min_dist = 15
            for node in self.canvas.nodes:
                dist = (pos - node.pos).manhattanLength()
                if dist < min_dist:
                    min_dist = dist
                    snapped_node = node
            
            if snapped_node:
                if not self.canvas.temp_beam_start:
                    self.canvas.temp_beam_start = snapped_node
                else:
                    beam = BeamItem(self.canvas.temp_beam_start, snapped_node)
                    self.canvas.scene.addItem(beam)
                    self.canvas.beams.append(beam)
                    self.canvas.temp_beam_start = None
                    
        elif self.canvas.modes['delete']:
            items = self.canvas.scene.items(pos)
            for item in items:
                if isinstance(item, NodeItem):
                    self.canvas.delete_node(item)
                    self.activate_add_node_mode()  # Exit delete mode
                    break
                elif isinstance(item, BeamItem):
                    self.canvas.delete_beam(item)
                    self.activate_add_node_mode()  # Exit delete mode
                    break
                    
        super(Canvas, self.canvas).mousePressEvent(event)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())