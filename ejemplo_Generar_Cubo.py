from pyautocad import Autocad, APoint, aDouble

# Initialize AutoCAD application
acad = Autocad(create_if_not_exists=True)

# Define the coordinates for the cube's vertices
vertices = [
    (0, 0, 0),  # Vertex 0
    (1, 0, 0),  # Vertex 1
    (1, 1, 0),  # Vertex 2
    (0, 1, 0),  # Vertex 3
    (0, 0, 1),  # Vertex 4
    (1, 0, 1),  # Vertex 5
    (1, 1, 1),  # Vertex 6
    (0, 1, 1)   # Vertex 7
]

# Define each edge of the cube as a polyline
edges = [
    (0, 1),  # Edge from Vertex 0 to Vertex 1
    (1, 2),  # Edge from Vertex 1 to Vertex 2
    (2, 3),  # Edge from Vertex 2 to Vertex 3
    (3, 0),  # Edge from Vertex 3 to Vertex 0
    (4, 5),  # Edge from Vertex 4 to Vertex 5
    (5, 6),  # Edge from Vertex 5 to Vertex 6
    (6, 7),  # Edge from Vertex 6 to Vertex 7
    (7, 4),  # Edge from Vertex 7 to Vertex 4
    (0, 4),  # Edge from Vertex 0 to Vertex 4
    (1, 5),  # Edge from Vertex 1 to Vertex 5
    (2, 6),  # Edge from Vertex 2 to Vertex 6
    (3, 7)   # Edge from Vertex 3 to Vertex 7
]

# Create polylines for each edge
for edge in edges:
    start_vertex = vertices[edge[0]]
    end_vertex = vertices[edge[1]]
    polilinea = aDouble(start_vertex[0], start_vertex[1], start_vertex[2],
                        end_vertex[0], end_vertex[1], end_vertex[2])
    acad.model.Add3Dpoly(polilinea)

print("Cube created with each edge as a separate polyline.")