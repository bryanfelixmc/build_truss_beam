from pyautocad import Autocad, aDouble
import numpy as np
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
    def __init__(self, acad, length, width, height, num_crosses,nivel_de_conexion,A,B):
        self.acad = acad
        self.length = length
        self.width = width
        self.height = height
        self.num_crosses = num_crosses
        self.nivel_de_conexion=nivel_de_conexion
        self.A=np.append(A, nivel_de_conexion) #vector inicio de viga
        self.B=np.append(B, nivel_de_conexion) #vector fin de viga

    def parametrizacion(self):
        AB=self.B-self.A
        if np.allclose(AB, np.zeros(3)):
            return [self.A.tolist()]
        u_raw = np.cross(AB, np.array([1, 0, 0]))

        if np.allclose(u_raw, np.zeros(3)):
            u_raw = np.cross(AB, np.array([0, 1, 0]))
        
        if np.allclose(u_raw, np.zeros(3)):
            u_raw = np.cross(AB, np.array([0, 0, 1]))
        
        u_hat = u_raw / np.linalg.norm(u_raw)
        u = u_hat * (self.width *0.5)
        
        v_raw = np.cross(AB, u_raw)
        v_hat = v_raw / np.linalg.norm(v_raw)
        v = v_hat * (self.height*0.5)


        vertices = [
            (self.A - u - v).tolist(),  # Vertex 0
            (self.A + AB - u - v).tolist(),  # Vertex 1
            (self.A + u - v).tolist(),  # Vertex 2
            (self.A + v - u).tolist(),  # Vertex 3
            (self.A + AB + u - v).tolist(),  # Vertex 4
            (self.A + AB + v - u).tolist(),  # Vertex 5
            (self.A + u + v).tolist(),  # Vertex 6
            (self.A + AB + u + v).tolist()  # Vertex 7
        ]
        if len(vertices) == 1:
            print("The parallelepiped collapses to a single point.")

        # Define the edges of the parallelepiped
        edges = [
            (0, 1),  # A to A + AB
            (0, 2),  # A to A + u
            (0, 3),  # A to A + v
            (1, 4),  # A + AB to A + AB + u
            (1, 5),  # A + AB to A + AB + v
            (2, 4),  # A + u to A + AB + u
            (2, 6),  # A + u to A + u + v
            (3, 5),  # A + v to A + AB + v
            (3, 6),  # 
            (4, 7),  # A + AB + u to A + u + v
            (5, 7),  # 
            (6, 7)   # A + u + v to A + AB + u + v
        ]
        return [vertices,edges]

    def dibujar(self):
        V = self.parametrizacion()[0] 
        print(V)
        edges = self.parametrizacion()[1]
        print(edges)
        for i, j in edges:
            self.acad.model.Add3Dpoly(aDouble(*V[i],*V[j]))
        self.construir_beam_truss()
    
    def construir_beam_truss(self):
        V = self.parametrizacion()[0]
        edges = self.parametrizacion()[1]

        # Identify edges that are parallel to AB (length of the beam)
        ab_edges = [(0, 1), (2, 4), (3, 5), (6, 7)]
        
        # Identify faces that are parallel to AB (length of the beam)
        #faces= [[(0, 1), (2, 4)], [(2, 4), (6, 7)], [(3, 5), (6, 7)], [(0, 1), (3, 5)]] #face1,face2,face3,face4
        faces= [[(0, 1), (3, 5)], [(3, 5), (6, 7)], [(6, 7), (2, 4)], [(2, 4), (0, 1)]] #face1,face2,face3,face4
        for face in faces:
            p1=np.array(V[face[0][0]])
            p2=np.array(V[face[0][1]])
            q1=np.array(V[face[1][0]])
            q2=np.array(V[face[1][1]])
            vv=(p2-p1)/self.num_crosses #vector que va longitudinalmente, y de tamaño igual a las particiones 
            m = q1-p1+vv #vv / np.linalg.norm(vv)  :  vector en direccion hacia arriba, en la misma cara
            n=vv-q1+p1   # : vector en direccion hacia abajo, en la misma cara
            #Debe formarse la siguiente secuencia de ordenamiento para los nuevos puntos del truss
            #CASO self.num_crosses =5
            #A=p1+(0)*m+(0)*n
            #B=p1+(1)*m+(0)*n
            #C=p1+(1)*m+(1)*n
            #D=p1+(2)*m+(1)*n
            #E=p1+(2)*m+(2)*n
            #F=p1+(3)*m+(2)*n

            i1=0
            i2=0
            aux=True
            lista_3dpoly_truss=[]
            while i1+i2<=self.num_crosses:
                if aux:
                    pf=p1+i1*m+i2*n
                    i1=i1+1
                    #lista_3dpoly_truss.append(pf)
                    lista_3dpoly_truss.extend(pf.tolist()) 
                    aux=not aux
                else:
                    pf=p1+i1*m+i2*n
                    i2=i2+1
                    #lista_3dpoly_truss.append(pf)
                    lista_3dpoly_truss.extend(pf.tolist()) 
                    aux=not aux
            self.acad.model.Add3Dpoly(aDouble(*lista_3dpoly_truss))

def main():
    acad = Autocad(create_if_not_exists=True)
    #--------Parametros Generales Pórtico------
    ancho_de_campo=12
    altura_de_nivel_de_conexion=4
    punto_inicio_portico=np.array([0,0])
    punto_fin_portico=np.array([0,ancho_de_campo])


    #-----------------COLUMNA------------------
    # Tronco de piramide
    largo_base = 2 #Largo de la base del tronco
    ancho_base = 1 #Ancho de la base del tronco
    largo_top = 0.85 #Largo de la cara superior del tronco (base columna)
    ancho_top = 1 #Ancho de la cara superior del tronco (base columna) #--$
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
    aa=punto_inicio_portico
    bb=punto_fin_portico
    punto_inicio_viga=aa+((bb-aa)/np.linalg.norm(bb-aa))*largo_top*0.5
    punto_fin_viga=bb-((bb-aa)/np.linalg.norm(bb-aa))*largo_top*0.5

    viga1 = viga(acad, length, width, height, num_crosses_box,altura_de_nivel_de_conexion,punto_inicio_viga,punto_fin_viga)
    viga1.dibujar()

    print("Viga generada en AutoCAD.")

if __name__ == "__main__":
    main()
#Falta parametrizar clase columna
#clase viga ya esta parametrizada