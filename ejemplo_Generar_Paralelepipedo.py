from pyautocad import Autocad, aDouble

def create_parallelepiped(length, width, height):
    # Initialize AutoCAD application
    acad = Autocad(create_if_not_exists=True)

    # Define the coordinates for the parallelepiped's vertices
    vertices = [
        (0, 0, 0),          # Vertex 0
        (length, 0, 0),     # Vertex 1
        (length, width, 0), # Vertex 2
        (0, width, 0),      # Vertex 3
        (0, 0, height),     # Vertex 4 (height)
        (length, 0, height), # Vertex 5 (height)
        (length, width, height), # Vertex 6 (height)
        (0, width, height)   # Vertex 7 (height)
    ]

    # Define each edge of the parallelepiped as a polyline
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

    print("Parallelepiped created with each edge as a separate polyline.")

# Example usage
create_parallelepiped(length=5, width=1, height=1)  # You can change the dimensions as needed