from pyautocad import Autocad, aDouble
class columna:
    def __init__(self,acad,largo_base, ancho_base,largo_top,ancho_top,altura_tronco,cruces_tronco,altura_columna,cruces_columna,altura_piramide,cruces_piramide,coord_xy):
        self.acad=acad
        self.largo_base=largo_base
        self.ancho_base=ancho_base
        self.largo_top=largo_top
        self.ancho_top=ancho_top
        self.altura_tronco=altura_tronco
        self.cruces_tronco=cruces_tronco
        self.altura_columna=altura_columna
        self.cruces_columna=cruces_columna
        self.altura_piramide=altura_piramide
        self.cruces_piramide=cruces_piramide
        self.coord_xy=coord_xy

    def interpolate(self,p1, p2, t):
        return (
            p1[0] + (p2[0] - p1[0]) * t,
            p1[1] + (p2[1] - p1[1]) * t,
            p1[2] + (p2[2] - p1[2]) * t
        )

    def crear_tronco_piramide(self):
        ''' Parametros
            largo_base, ancho_base: Dimensiones de la cara inferior.
            largo_top, ancho_top: Dimensiones de cara superior.
            altura_tronco: Altura del tronco (Altura total del tronco).
            num_secciones: Numero de secciones de corte
            
            Variables
            altura_nivel: Altura de cada seccion de corte.
            z0, z1: Las z-coordinates para la actual y la siguiente seccion de corte.
            l0, a0, l1, a1: Dimensiones de la actual y la siguiente seccion de corte.
            cx0, cy0, cx1, cy1: Cordenadas del centro de la actual y la siguiente seccion de corte .
            v: A dictionary containing the vertices (A, B, C, D, E, F, G, H) of the current cross-section, where each vertex is represented as a tuple of (x, y, z).
            ''' 
        altura_nivel = self.altura_tronco / self.cruces_tronco
        aux=[]
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

    def crear_columna_con_cruces(self):

        ''' Parameters:
            largo, ancho: Dimensions of the column.
            altura_columna: Height of the column.
            num_secciones: Number of cross-sections.
            altura_base: The base height from which the column starts.
            
            Variables:
            Similar to crear_tronco_piramide, it calculates the z-coordinates and dimensions for each cross-section and stores the vertices in the dictionary v.'''
                
        altura_nivel = self.altura_columna / self.cruces_columna
        aux=[]
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
        '''
        Parameters:
        base_length, base_width: Dimensions of the base of the pyramid.
        height: Height of the pyramid.
        num_crosses: Number of cross-sections on the pyramid.
        altura_base: The base height from which the pyramid starts.

        Variables:
        half_length, half_width: Half dimensions of the base.
        base_vertices: List of vertices for the base of the pyramid.
        apex: The apex (top point) of the pyramid.
        The function draws lines between the base vertices and the apex, as well as cross-sections.
        '''
        
        half_length = self.largo_top / 2
        half_width = self.ancho_top / 2

        base_vertices = [
            (self.coord_xy[0]-half_length, self.coord_xy[1]-half_width, self.altura_tronco + self.altura_columna),  # V0
            (self.coord_xy[0]+half_length, self.coord_xy[1]-half_width, self.altura_tronco + self.altura_columna),   # V1
            (self.coord_xy[0]+half_length, self.coord_xy[1]+half_width, self.altura_tronco + self.altura_columna),    # V2
            (self.coord_xy[0]-half_length, self.coord_xy[1]+half_width, self.altura_tronco + self.altura_columna)    # V3
        ]
        apex = (0 + self.coord_xy[0], 0 + self.coord_xy[1], self.altura_tronco + self.altura_columna + self.altura_piramide)

        vertices = base_vertices + [apex]
        

        return [base_vertices,apex,vertices]

    def dibujar(self):
        #dibujar tronco piramide
        aux=self.crear_tronco_piramide()
        for v in aux:
            for p1, p2 in [("A", "E"), ("B", "F"), ("C", "G"), ("D", "H")]:
                self.acad.model.Add3Dpoly(aDouble(*v[p1], *v[p2]))

            for p1, p2, p3, p4 in [("A", "F", "B", "E"), ("B", "G", "C", "F"), ("C", "H", "D", "G"), ("D", "E", "A", "H")]:
                self.acad.model.Add3Dpoly(aDouble(*v[p1], *v[p2]))
                self.acad.model.Add3Dpoly(aDouble(*v[p3], *v[p4]))
        
        #dibujar  columna intermedia
        aux=self.crear_columna_con_cruces()
        for v in aux:
            for p1, p2 in [("A", "E"), ("B", "F"), ("C", "G"), ("D", "H")]:
                self.acad.model.Add3Dpoly(aDouble(*v[p1], *v[p2]))

            for p1, p2, p3, p4 in [("A", "F", "B", "E"), ("B", "G", "C", "F"), ("C", "H", "D", "G"), ("D", "E", "A", "H")]:
                self.acad.model.Add3Dpoly(aDouble(*v[p1],*v[p2]))
                self.acad.model.Add3Dpoly(aDouble(*v[p3], *v[p4]))

        #dibujar  castillete
        base_vertices=self.crear_piramide_mastil()[0]
        apex=self.crear_piramide_mastil()[1]
        vertices=self.crear_piramide_mastil()[2]

        for i, j in [(0,1),(1,2),(2,3),(3,0),(0,4),(1,4),(2,4),(3,4)]: # Dibujar aristas principales
            self.acad.model.Add3Dpoly(aDouble(*vertices[i],*vertices[j]))
 
        cara_indices = [(0,1), (1,2), (2,3), (3,0)] # Dibujar num_secciones en caras laterales
        for i1, i2 in cara_indices:
            base1 = base_vertices[i1]
            base2 = base_vertices[i2]

            for i in range(1, self.cruces_piramide + 1):
                t1 = (i - 1) / (self.cruces_piramide + 1)
                t2 = i / (self.cruces_piramide + 1)

                p1_start = self.interpolate(base1, apex, t1)
                p1_end = self.interpolate(base2, apex, t2)
                self.acad.model.Add3Dpoly(aDouble(*p1_start,*p1_end))

                p2_start = self.interpolate(base2, apex, t1)
                p2_end = self.interpolate(base1, apex, t2)
                self.acad.model.Add3Dpoly(aDouble(*p2_start,*p2_end))

class viga:
    def __init__(self, acad, length, width, height, num_crosses,nivel_de_conexion):
        self.acad = acad
        self.length = length
        self.width = width
        self.height = height
        self.num_crosses = num_crosses
        self.nivel_de_conexion=nivel_de_conexion

    def dibujar(self):
        V = [
            (0, 0, 0+self.nivel_de_conexion),
            (self.length, 0, 0+self.nivel_de_conexion),
            (self.length, self.width, 0+self.nivel_de_conexion),
            (0, self.width, 0+self.nivel_de_conexion),
            (0, 0, self.height+self.nivel_de_conexion),
            (self.length, 0, self.height+self.nivel_de_conexion),
            (self.length, self.width, self.height+self.nivel_de_conexion),
            (0, self.width, self.height+self.nivel_de_conexion)
        ]

        edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
        for i, j in edges:
            self.acad.model.Add3Dpoly(aDouble(*V[i],*V[j]))

        self._agregar_cruces(0, self.length, 0, 'XZ')
        self._agregar_cruces(0, self.length, self.width, 'XZ')
        self._agregar_cruces(0, self.length, 0, 'XY')
        self._agregar_cruces(0, self.length, self.height, 'XY')

    def _agregar_cruces(self, start, end, fixed, axis):
        step = (end - start) / self.num_crosses
        for i in range(self.num_crosses):
            a = start + i * step
            b = start + (i + 1) * step

            if axis == 'XZ':
                y = fixed
                z0 = 0+self.nivel_de_conexion
                z1 = self.height+self.nivel_de_conexion
                self.acad.model.Add3Dpoly(aDouble(a, y, z0,b, y, z1))
                self.acad.model.Add3Dpoly(aDouble(a, y, z1,b, y, z0))

            elif axis == 'XY':
                z = fixed+self.nivel_de_conexion
                self.acad.model.Add3Dpoly(aDouble(a, 0, z, b, self.width, z))
                self.acad.model.Add3Dpoly(aDouble(a, self.width, z, b, 0, z))


def main():
    acad = Autocad(create_if_not_exists=True)
    #--------Parametros Generales Pórtico------
    ancho_de_campo=12
    altura_de_nivel_de_conexion=4
    punto_inicio_portico=[0,0]
    punto_fin_portico=[ancho_de_campo,0]


    #-----------------COLUMNA------------------
    # Tronco de piramide
    largo_base = 2 #Largo de la base del tronco
    ancho_base = 1 #Ancho de la base del tronco
    largo_top = 1 #Largo de la cara superior del tronco (base columna)
    ancho_top = 0.5 #Ancho de la cara superior del tronco (base columna) #--$
    altura_tronco = 3 #Altura del tronco
    cruces_tronco = 3 #Numero de secciones en  tronco

    # Columna Intermedia
    altura_columna =4 # Altura de  columna
    cruces_columna = 4 # Numero de secciones en  columna

    # Castillete
    altura_piramide =2 #Altura del castillete
    cruces_piramide = 2 #Numero de secciones por cara en castillete

    columna1 = columna(acad, largo_base, ancho_base,largo_top,ancho_top,altura_tronco,cruces_tronco,altura_columna,cruces_columna,altura_piramide,cruces_piramide,punto_inicio_portico)
    columna1.dibujar()
    columna2 = columna(acad, largo_base, ancho_base,largo_top,ancho_top,altura_tronco,cruces_tronco,altura_columna,cruces_columna,altura_piramide,cruces_piramide,punto_fin_portico)
    columna2.dibujar()
    print("Columnas Generada en AutoCAD")

    #-------------------VIGA---------------------
    #Ingrese dimensiones de la viga
    length = ancho_de_campo #Largo
    width = ancho_top #Ancho #--$
    height = ancho_top #Alto (modificable)
    altura_de_nivel_de_conexion=5 #Nivel de conexion
    num_crosses_box = ancho_de_campo #Cantidad de cruces para la viga (modificable)

    viga1 = viga(acad, length, width, height, num_crosses_box,altura_de_nivel_de_conexion)
    viga1.dibujar()

    print("Viga generada en AutoCAD.")

if __name__ == "__main__":
    main()
