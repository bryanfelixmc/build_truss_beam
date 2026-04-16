from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QMainWindow, QOpenGLWidget


import sys
import os

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)


try:
    from OpenGL import GL, GLU
except Exception as exc:  # pragma: no cover - ayuda en runtime del usuario
    GL = None
    GLU = None
    _OPENGL_IMPORT_ERROR = exc
else:
    _OPENGL_IMPORT_ERROR = None

try:
    import ezdxf
except Exception:  # pragma: no cover - se maneja en exportación
    ezdxf = None


Point3D = Tuple[float, float, float]
Primitive3D = List[Point3D]


def as_point3(p: Sequence[float]) -> Point3D:
    """Convierte listas / np.ndarray / tuplas a una tupla de 3 floats."""
    return (float(p[0]), float(p[1]), float(p[2]))


@dataclass
class PorticoModel:
    """Modelo geométrico independiente del motor de render."""
    primitives: List[Primitive3D] = field(default_factory=list)

    def clear(self) -> None:
        self.primitives.clear()

    def add_segment(self, p1: Sequence[float], p2: Sequence[float]) -> None:
        self.primitives.append([as_point3(p1), as_point3(p2)])

    def add_polyline(self, points: Iterable[Sequence[float]]) -> None:
        pts = [as_point3(p) for p in points]
        if len(pts) >= 2:
            self.primitives.append(pts)

    def all_points(self) -> List[Point3D]:
        pts: List[Point3D] = []
        for prim in self.primitives:
            pts.extend(prim)
        return pts

    def bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        pts = self.all_points()
        if not pts:
            return np.zeros(3), np.ones(3)
        arr = np.array(pts, dtype=float)
        return arr.min(axis=0), arr.max(axis=0)


class DXFExporter:
    """Exporta el modelo geométrico a DXF usando ezdxf."""

    def __init__(self, model: PorticoModel) -> None:
        self.model = model

    def export(self, filename: str) -> None:
        if ezdxf is None:
            raise RuntimeError(
                "No se pudo importar 'ezdxf'. Instala primero: pip install ezdxf"
            )

        doc = ezdxf.new(setup=True)
        msp = doc.modelspace()

        for primitive in self.model.primitives:
            if len(primitive) == 2:
                msp.add_line(primitive[0], primitive[1])
            elif len(primitive) > 2:
                msp.add_polyline3d(primitive)

        doc.saveas(filename)


class Columna:
    def __init__(
        self,
        model: PorticoModel,
        largo_base: float,
        ancho_base: float,
        largo_top: float,
        ancho_top: float,
        altura_tronco: float,
        cruces_tronco: int,
        altura_columna: float,
        cruces_columna: int,
        altura_piramide: float,
        cruces_piramide: int,
        coord_xy: Sequence[float],
    ) -> None:
        self.model = model
        self.largo_base = largo_base
        self.ancho_base = ancho_base
        self.largo_top = largo_top
        self.ancho_top = ancho_top
        self.altura_tronco = altura_tronco
        self.cruces_tronco = cruces_tronco
        self.altura_columna = altura_columna
        self.cruces_columna = cruces_columna
        self.altura_piramide = altura_piramide
        self.cruces_piramide = cruces_piramide
        self.coord_xy = np.array(coord_xy, dtype=float)

    @staticmethod
    def interpolate(p1: Sequence[float], p2: Sequence[float], t: float) -> Point3D:
        p1 = np.array(p1, dtype=float)
        p2 = np.array(p2, dtype=float)
        p = p1 + (p2 - p1) * t
        return as_point3(p)

    def crear_tronco_piramide(self) -> List[dict]:
        altura_nivel = self.altura_tronco / self.cruces_tronco
        aux = []
        for i in range(self.cruces_tronco):
            z0 = i * altura_nivel
            z1 = (i + 1) * altura_nivel

            l0 = self.largo_base - (self.largo_base - self.largo_top) * (i / self.cruces_tronco)
            a0 = self.ancho_base - (self.ancho_base - self.ancho_top) * (i / self.cruces_tronco)
            l1 = self.largo_base - (self.largo_base - self.largo_top) * ((i + 1) / self.cruces_tronco)
            a1 = self.ancho_base - (self.ancho_base - self.ancho_top) * ((i + 1) / self.cruces_tronco)

            cx0 = -l0 / 2 + self.coord_xy[0]
            cy0 = -a0 / 2 + self.coord_xy[1]
            cx1 = -l1 / 2 + self.coord_xy[0]
            cy1 = -a1 / 2 + self.coord_xy[1]

            v = {
                "A": (cx0, cy0, z0),
                "B": (cx0 + l0, cy0, z0),
                "C": (cx0 + l0, cy0 + a0, z0),
                "D": (cx0, cy0 + a0, z0),
                "E": (cx1, cy1, z1),
                "F": (cx1 + l1, cy1, z1),
                "G": (cx1 + l1, cy1 + a1, z1),
                "H": (cx1, cy1 + a1, z1),
            }
            aux.append(v)
        return aux

    def crear_columna_con_cruces(self) -> List[dict]:
        altura_nivel = self.altura_columna / self.cruces_columna
        aux = []
        for i in range(self.cruces_columna):
            z0 = self.altura_tronco + i * altura_nivel
            z1 = self.altura_tronco + (i + 1) * altura_nivel

            cx = -self.largo_top / 2 + self.coord_xy[0]
            cy = -self.ancho_top / 2 + self.coord_xy[1]

            v = {
                "A": (cx, cy, z0),
                "B": (cx + self.largo_top, cy, z0),
                "C": (cx + self.largo_top, cy + self.ancho_top, z0),
                "D": (cx, cy + self.ancho_top, z0),
                "E": (cx, cy, z1),
                "F": (cx + self.largo_top, cy, z1),
                "G": (cx + self.largo_top, cy + self.ancho_top, z1),
                "H": (cx, cy + self.ancho_top, z1),
            }
            aux.append(v)
        return aux

    def crear_piramide_mastil(self):
        half_length = self.largo_top / 2
        half_width = self.ancho_top / 2

        base_vertices = [
            (self.coord_xy[0] - half_length, self.coord_xy[1] - half_width, self.altura_tronco + self.altura_columna),
            (self.coord_xy[0] + half_length, self.coord_xy[1] - half_width, self.altura_tronco + self.altura_columna),
            (self.coord_xy[0] + half_length, self.coord_xy[1] + half_width, self.altura_tronco + self.altura_columna),
            (self.coord_xy[0] - half_length, self.coord_xy[1] + half_width, self.altura_tronco + self.altura_columna),
        ]
        apex = (
            self.coord_xy[0],
            self.coord_xy[1],
            self.altura_tronco + self.altura_columna + self.altura_piramide,
        )
        vertices = base_vertices + [apex]
        return [base_vertices, apex, vertices]

    def generate_geometry(self) -> None:
        # Tronco de pirámide
        for v in self.crear_tronco_piramide():
            for p1, p2 in [("A", "E"), ("B", "F"), ("C", "G"), ("D", "H")]:
                self.model.add_segment(v[p1], v[p2])

            for p1, p2, p3, p4 in [
                ("A", "F", "B", "E"),
                ("B", "G", "C", "F"),
                ("C", "H", "D", "G"),
                ("D", "E", "A", "H"),
            ]:
                self.model.add_segment(v[p1], v[p2])
                self.model.add_segment(v[p3], v[p4])

        # Columna intermedia
        for v in self.crear_columna_con_cruces():
            for p1, p2 in [("A", "E"), ("B", "F"), ("C", "G"), ("D", "H")]:
                self.model.add_segment(v[p1], v[p2])

            for p1, p2, p3, p4 in [
                ("A", "F", "B", "E"),
                ("B", "G", "C", "F"),
                ("C", "H", "D", "G"),
                ("D", "E", "A", "H"),
            ]:
                self.model.add_segment(v[p1], v[p2])
                self.model.add_segment(v[p3], v[p4])

        # Castillete
        base_vertices, apex, vertices = self.crear_piramide_mastil()

        for i, j in [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)]:
            self.model.add_segment(vertices[i], vertices[j])

        cara_indices = [(0, 1), (1, 2), (2, 3), (3, 0)]
        for i1, i2 in cara_indices:
            base1 = base_vertices[i1]
            base2 = base_vertices[i2]
            for i in range(1, self.cruces_piramide + 1):
                t1 = (i - 1) / (self.cruces_piramide + 1)
                t2 = i / (self.cruces_piramide + 1)

                p1_start = self.interpolate(base1, apex, t1)
                p1_end = self.interpolate(base2, apex, t2)
                self.model.add_segment(p1_start, p1_end)

                p2_start = self.interpolate(base2, apex, t1)
                p2_end = self.interpolate(base1, apex, t2)
                self.model.add_segment(p2_start, p2_end)


class Viga:
    def __init__(
        self,
        model: PorticoModel,
        length: float,
        width_viga: float,
        height_viga: float,
        num_crosses: int,
        nivel_de_conexion: float,
        A: Sequence[float],
        B: Sequence[float],
    ) -> None:
        self.model = model
        self.length = float(length)
        self.width_viga = float(width_viga)
        self.height_viga = float(height_viga)
        self.num_crosses = int(num_crosses)
        self.nivel_de_conexion = float(nivel_de_conexion)
        self.A = np.append(np.array(A, dtype=float), self.nivel_de_conexion)
        self.B = np.append(np.array(B, dtype=float), self.nivel_de_conexion)

    def parametrizacion(self):
        AB = self.B - self.A
        if np.allclose(AB, np.zeros(3)):
            return [[self.A.tolist()], []]

        u_raw = np.cross(AB, np.array([1.0, 0.0, 0.0]))
        if np.allclose(u_raw, np.zeros(3)):
            u_raw = np.cross(AB, np.array([0.0, 1.0, 0.0]))
        if np.allclose(u_raw, np.zeros(3)):
            u_raw = np.cross(AB, np.array([0.0, 0.0, 1.0]))

        u_hat = u_raw / np.linalg.norm(u_raw)
        u = u_hat * (self.width_viga * 0.5)

        v_raw = np.cross(AB, u_raw)
        v_hat = v_raw / np.linalg.norm(v_raw)
        v = v_hat * (self.height_viga * 0.5)

        vertices = [
            (self.A - u - v).tolist(),
            (self.A + AB - u - v).tolist(),
            (self.A + u - v).tolist(),
            (self.A + v - u).tolist(),
            (self.A + AB + u - v).tolist(),
            (self.A + AB + v - u).tolist(),
            (self.A + u + v).tolist(),
            (self.A + AB + u + v).tolist(),
        ]

        edges = [
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 4),
            (1, 5),
            (2, 4),
            (2, 6),
            (3, 5),
            (3, 6),
            (4, 7),
            (5, 7),
            (6, 7),
        ]
        return [vertices, edges]

    def generate_geometry(self) -> None:
        vertices, edges = self.parametrizacion()
        if len(vertices) == 1:
            return

        for i, j in edges:
            self.model.add_segment(vertices[i], vertices[j])

        self.construir_beam_truss()

    def construir_beam_truss(self) -> None:
        vertices, _ = self.parametrizacion()

        faces = [[(0, 1), (3, 5)], [(3, 5), (6, 7)], [(6, 7), (2, 4)], [(2, 4), (0, 1)]]
        for face in faces:
            p1 = np.array(vertices[face[0][0]], dtype=float)
            p2 = np.array(vertices[face[0][1]], dtype=float)
            q1 = np.array(vertices[face[1][0]], dtype=float)
            q2 = np.array(vertices[face[1][1]], dtype=float)

            vv = (p2 - p1) / self.num_crosses
            m = q1 - p1 + vv
            n = vv - q1 + p1

            i1 = 0
            i2 = 0
            aux = True
            lista_3dpoly_truss = []

            while i1 + i2 <= self.num_crosses:
                if aux:
                    pf = p1 + i1 * m + i2 * n
                    i1 += 1
                else:
                    pf = p1 + i1 * m + i2 * n
                    i2 += 1
                lista_3dpoly_truss.append(as_point3(pf))
                aux = not aux

            self.model.add_polyline(lista_3dpoly_truss)


class PorticoBuilder:
    """Genera el modelo completo del pórtico usando las mismas fórmulas geométricas."""

    def __init__(self) -> None:
        self.model = PorticoModel()

    def generar(self, params: dict | None = None) -> PorticoModel:
        self.model.clear()

        defaults = {
            "ancho_de_campo": 12.0,
            "altura_de_nivel_de_conexion": 5.0,
            "largo_base": 2.0,
            "ancho_base": 1.0,
            "largo_top": 0.85,
            "ancho_top": 1.0,
            "altura_tronco": 3.0,
            "cruces_tronco": 3,
            "altura_columna": 4.0,
            "cruces_columna": 4,
            "altura_piramide": 2.0,
            "cruces_piramide": 2,
            "num_crosses_box": 12,
        }

        if params is None:
            params = defaults
        else:
            params = {**defaults, **params}

        #--------Parametros Generales Pórtico------
        ancho_de_campo = float(params["ancho_de_campo"])
        altura_de_nivel_de_conexion = float(params["altura_de_nivel_de_conexion"])
        punto_inicio_portico = np.array([0.0, 0.0])
        punto_fin_portico = np.array([0.0, float(ancho_de_campo)])

        #-----------------COLUMNA------------------
        # Tronco de piramide
        largo_base = float(params["largo_base"])  # Largo de la base del tronco
        ancho_base = float(params["ancho_base"])  # Ancho de la base del tronco
        largo_top = float(params["largo_top"])  # Largo de la cara superior del tronco (base columna)
        ancho_top = float(params["ancho_top"])  # Ancho de la cara superior del tronco (base columna)
        altura_tronco = float(params["altura_tronco"])  # Altura del tronco
        cruces_tronco = int(params["cruces_tronco"])  # Numero de secciones en  tronco

        # Columna Intermedia
        altura_columna = float(params["altura_columna"])  # Altura de  columna
        cruces_columna = int(params["cruces_columna"])  # Numero de secciones en  columna

        # Castillete
        altura_piramide = float(params["altura_piramide"])  # Altura del castillete
        cruces_piramide = int(params["cruces_piramide"])  # Numero de secciones por cara en castillete

        # Construccion de columnas de pórtico
        columna1 = Columna(
            self.model,
            largo_base, ancho_base, largo_top, ancho_top,
            altura_tronco, cruces_tronco,
            altura_columna, cruces_columna,
            altura_piramide, cruces_piramide,
            punto_inicio_portico,
        )
        columna1.generate_geometry()

        columna2 = Columna(
            self.model,
            largo_base, ancho_base, largo_top, ancho_top,
            altura_tronco, cruces_tronco,
            altura_columna, cruces_columna,
            altura_piramide, cruces_piramide,
            punto_fin_portico,
        )
        columna2.generate_geometry()

        #-------------------VIGA---------------------
        # Ingrese dimensiones de la viga
        length = ancho_de_campo  # Largo
        width_viga = float(params["ancho_top"])  # Ancho (vinculado a ancho_top)
        height_viga = float(params["ancho_top"])  # Alto (vinculado a ancho_top)

        num_crosses_box = int(params["num_crosses_box"])  # Cantidad de cruces para la viga
        aa = punto_inicio_portico
        bb = punto_fin_portico
        punto_inicio_viga = aa + ((bb - aa) / np.linalg.norm(bb - aa)) * largo_top * 0.5
        punto_fin_viga = bb - ((bb - aa) / np.linalg.norm(bb - aa)) * largo_top * 0.5

        # Construccion de viga de pórtico
        viga1 = Viga(
            self.model,
            length, width_viga, height_viga,
            num_crosses_box,
            altura_de_nivel_de_conexion,
            punto_inicio_viga,
            punto_fin_viga,
        )
        viga1.generate_geometry()

        return self.model

class PorticosGLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.model = PorticoModel()
        self.origin = np.zeros(3, dtype=float)
        self.center = np.zeros(3, dtype=float)
        self.camera_focus = np.zeros(3, dtype=float)
        self.radius = 10.0
        self.yaw = 45.0
        self.pitch = 20.0
        self.distance = 30.0
        self._last_pos = None
        self.setMinimumSize(300, 300)

        self.pan_active = False

    def set_model(self, model: PorticoModel) -> None:
        self.model = model
        mn, mx = self.model.bounds()
        self.center = (mn + mx) / 2.0
        self.origin = np.zeros(3, dtype=float)
        self.camera_focus = self.center.copy()
        span = np.linalg.norm(mx - mn)
        self.radius = max(float(span), 1.0) * 0.55
        self.distance = max(self.radius * 2.5, 10.0)
        self.update()

    def set_camera_focus_to_origin(self) -> None:
        self.camera_focus = self.origin.copy()
        self.update()

    def set_camera_focus_to_center(self) -> None:
        self.camera_focus = self.center.copy()
        self.update()

    def initializeGL(self) -> None:
        if GL is None:
            raise RuntimeError(
                "No se pudo importar PyOpenGL. Instala 'PyOpenGL' y 'PyOpenGL_accelerate' si aplica."
            ) from _OPENGL_IMPORT_ERROR

        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_LINE_SMOOTH)
        GL.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST)
        GL.glClearColor(0.10, 0.10, 0.12, 1.0)
        GL.glLineWidth(1)

    def resizeGL(self, w: int, h: int) -> None:
        if GL is None:
            return
        h = max(h, 1)
        GL.glViewport(0, 0, w, h)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        aspect = w / float(h)
        GLU.gluPerspective(45.0, aspect, 0.1, 1000.0)
        GL.glMatrixMode(GL.GL_MODELVIEW)

    def _camera_position(self) -> Tuple[float, float, float]:
        yaw = math.radians(self.yaw)
        pitch = math.radians(self.pitch)
        focus = self.camera_focus
        x = focus[0] + self.distance * math.cos(pitch) * math.sin(yaw)
        y = focus[1] + self.distance * math.cos(pitch) * math.cos(yaw)
        z = focus[2] + self.distance * math.sin(pitch)
        return x, y, z

    def _pan_camera(self, dx: float, dy: float) -> None:
        if dx == 0 and dy == 0:
            return

        focus = self.camera_focus
        eye = np.array(self._camera_position(), dtype=float)
        direction = focus - eye
        if np.linalg.norm(direction) < 1e-9:
            return
        direction /= np.linalg.norm(direction)

        world_up = np.array([0.0, 0.0, 1.0], dtype=float)
        right = np.cross(direction, world_up)
        if np.linalg.norm(right) < 1e-6:
            right = np.array([1.0, 0.0, 0.0], dtype=float)
        else:
            right /= np.linalg.norm(right)

        up = np.cross(right, direction)
        up /= np.linalg.norm(up)

        factor = self.radius * 0.0020
        self.camera_focus += (-dx * right + dy * up) * factor

    def paintGL(self) -> None:
        if GL is None:
            return

        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glLoadIdentity()

        eye = self._camera_position()
        focus = self.camera_focus
        GLU.gluLookAt(
            eye[0], eye[1], eye[2],
            focus[0], focus[1], focus[2],
            0.0, 0.0, 1.0,
        )

        # Ejes de referencia
        self._draw_axes()

        GL.glColor3f(0.90, 0.90, 0.95)
        for primitive in self.model.primitives:
            if len(primitive) < 2:
                continue
            GL.glBegin(GL.GL_LINE_STRIP)
            for x, y, z in primitive:
                GL.glVertex3f(float(x), float(y), float(z))
            GL.glEnd()

    def _draw_axes(self) -> None:
        if GL is None:
            return
        axis_len = max(self.radius * 1.2, 5.0)
        origin = self.origin
        GL.glBegin(GL.GL_LINES)

        GL.glColor3f(1.0, 0.3, 0.3)  # X
        GL.glVertex3f(origin[0], origin[1], origin[2])
        GL.glVertex3f(origin[0] + axis_len, origin[1], origin[2])

        GL.glColor3f(0.3, 1.0, 0.3)  # Y
        GL.glVertex3f(origin[0], origin[1], origin[2])
        GL.glVertex3f(origin[0], origin[1] + axis_len, origin[2])

        GL.glColor3f(0.3, 0.6, 1.0)  # Z
        GL.glVertex3f(origin[0], origin[1], origin[2])
        GL.glVertex3f(origin[0], origin[1], origin[2] + axis_len)

        GL.glEnd()
        GL.glColor3f(0.90, 0.90, 0.95)

    def mousePressEvent(self, event) -> None:
        self._last_pos = event.pos()
        if event.button() == QtCore.Qt.MiddleButton:
            self.pan_active = True

    def mouseReleaseEvent(self, event) -> None:
        if event.button() == QtCore.Qt.MiddleButton:
            self.pan_active = False
        self._last_pos = None

    def mouseMoveEvent(self, event) -> None:
        if self._last_pos is None:
            self._last_pos = event.pos()
            return

        dx = event.x() - self._last_pos.x()
        dy = event.y() - self._last_pos.y()

        if event.buttons() & QtCore.Qt.LeftButton:
            self.yaw += dx * 0.5
            self.pitch -= dy * 0.5
            self.pitch = max(-89.0, min(89.0, self.pitch))
            self.update()

        if self.pan_active:
            self._pan_camera(dx, dy)
            self.update()

        self._last_pos = event.pos()

    def wheelEvent(self, event) -> None:
        delta = event.angleDelta().y()
        factor = 0.85 if delta > 0 else 1.15
        self.distance = max(0.5, min(1000.0, self.distance * factor))
        self.update()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        ui_file = resource_path("main_ui_porticos_opengl_version05.ui")

        if not os.path.exists(ui_file):
            raise FileNotFoundError(f"No se encontró la interfaz: {ui_file}")
        uic.loadUi(ui_file, self)

        placeholder = self.findChild(QtWidgets.QOpenGLWidget, "openGLWidget")
        if placeholder is None:
            raise RuntimeError("La UI no contiene un widget llamado 'openGLWidget'.")

        geo = placeholder.geometry()
        parent = placeholder.parentWidget()

        placeholder.setParent(None)
        placeholder.deleteLater()

        self.gl_widget = PorticosGLWidget(parent)
        self.gl_widget.setGeometry(geo)
        self.gl_widget.show()

        self.pushButton.clicked.connect(self.actualizar)

        export_btn = self.findChild(QtWidgets.QPushButton, "pushButton_2")
        if export_btn is not None:
            export_btn.clicked.connect(self.exportar_dxf)

        self.builder = PorticoBuilder()
        self.model: PorticoModel | None = None

        self._configurar_campos()
        self._cargar_valores_iniciales()
        self._sincronizar_viga_con_ancho_top()

        if self.statusBar() is not None:
            self.statusBar().showMessage("Listo para generar el pórtico")

    def _configurar_campos(self) -> None:
        """Configura validadores y vínculos de la UI."""
        double_fields = [
            self.lineEdit_ancho_de_campo,
            self.lineEdit_altura_de_nivel_de_conexion,
            self.lineEdit_largo_base,
            self.lineEdit_ancho_base,
            self.lineEdit_largo_top,
            self.lineEdit_ancho_top,
            self.lineEdit_altura_tronco,
            self.lineEdit_altura_columna,
        ]
        int_fields = [
            self.lineEdit_cruces_tronco,
            self.lineEdit_cruces_columna,
            self.lineEdit_num_crosses_box,
        ]

        double_validator = QtGui.QDoubleValidator(0.0, 1e9, 6, self)
        double_validator.setNotation(QtGui.QDoubleValidator.StandardNotation)

        int_validator = QtGui.QIntValidator(1, 100000, self)

        for widget in double_fields:
            widget.setValidator(double_validator)
        for widget in int_fields:
            widget.setValidator(int_validator)

        self.lineEdit_width_viga.setReadOnly(True)
        self.lineEdit_height_viga.setReadOnly(True)

        self.lineEdit_ancho_top.textChanged.connect(self._sincronizar_viga_con_ancho_top)

    def _cargar_valores_iniciales(self) -> None:
        """Carga los valores de muestra en los QLineEdit."""
        defaults = {
            "lineEdit_ancho_de_campo": "12",
            "lineEdit_altura_de_nivel_de_conexion": "5",
            "lineEdit_largo_base": "2",
            "lineEdit_ancho_base": "1",
            "lineEdit_largo_top": "0.85",
            "lineEdit_ancho_top": "1",
            "lineEdit_altura_tronco": "3",
            "lineEdit_cruces_tronco": "3",
            "lineEdit_altura_columna": "4",
            "lineEdit_cruces_columna": "4",
            "lineEdit_num_crosses_box": "12",
        }

        for name, value in defaults.items():
            widget = getattr(self, name, None)
            if widget is not None:
                widget.setText(value)

    def _sincronizar_viga_con_ancho_top(self) -> None:
        """Vincula width_viga y height_viga con el valor de ancho_top."""
        ancho_top = self.lineEdit_ancho_top.text().strip()
        self.lineEdit_width_viga.setText(ancho_top)
        self.lineEdit_height_viga.setText(ancho_top)

    def _leer_float(self, widget_name: str, label: str) -> float:
        widget = getattr(self, widget_name)
        text = widget.text().strip()
        if not text:
            raise ValueError(f"El campo '{label}' está vacío.")
        return float(text)

    def _leer_int(self, widget_name: str, label: str) -> int:
        widget = getattr(self, widget_name)
        text = widget.text().strip()
        if not text:
            raise ValueError(f"El campo '{label}' está vacío.")
        return int(text)

    def leer_parametros(self) -> dict | None:
        """Lee los parámetros desde la UI y los devuelve en un diccionario."""
        try:
            params = {
                "ancho_de_campo": self._leer_float("lineEdit_ancho_de_campo", "ancho_de_campo"),
                "altura_de_nivel_de_conexion": self._leer_float(
                    "lineEdit_altura_de_nivel_de_conexion", "altura_de_nivel_de_conexion"
                ),
                "largo_base": self._leer_float("lineEdit_largo_base", "largo_base"),
                "ancho_base": self._leer_float("lineEdit_ancho_base", "ancho_base"),
                "largo_top": self._leer_float("lineEdit_largo_top", "largo_top"),
                "ancho_top": self._leer_float("lineEdit_ancho_top", "ancho_top"),
                "altura_tronco": self._leer_float("lineEdit_altura_tronco", "altura_tronco"),
                "cruces_tronco": self._leer_int("lineEdit_cruces_tronco", "cruces_tronco"),
                "altura_columna": self._leer_float("lineEdit_altura_columna", "altura_columna"),
                "cruces_columna": self._leer_int("lineEdit_cruces_columna", "cruces_columna"),
                "num_crosses_box": self._leer_int("lineEdit_num_crosses_box", "num_crosses_box"),
                "altura_piramide": 2.0,
                "cruces_piramide": 2,
            }
            return params
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, "Dato inválido", str(exc))
            return None

    def actualizar(self) -> None:
        params = self.leer_parametros()
        if params is None:
            return

        self.model = self.builder.generar(params)
        self.gl_widget.set_model(self.model)
        if self.statusBar() is not None:
            self.statusBar().showMessage(
                f"Modelo generado: {len(self.model.primitives)} primitivas geométricas"
            )

    def exportar_dxf(self) -> None:
        # Si no se ha generado aún, generamos antes de exportar.
        if self.model is None or not self.model.primitives:
            self.actualizar()

        if self.model is None or not self.model.primitives:
            if self.statusBar() is not None:
                self.statusBar().showMessage("No hay geometría para exportar")
            return

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Guardar DXF",
            "portico.dxf",
            "DXF Files (*.dxf)",
        )
        if not filename:
            return
        if not filename.lower().endswith(".dxf"):
            filename += ".dxf"

        try:
            DXFExporter(self.model).export(filename)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Error al exportar DXF", str(exc))
            if self.statusBar() is not None:
                self.statusBar().showMessage("Error al exportar DXF")
            return

        if self.statusBar() is not None:
            self.statusBar().showMessage(f"DXF exportado correctamente: {filename}")

def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
