from pyautocad import Autocad, aDouble

def create_pyramid(base_length, base_width, height):
    # Initialize AutoCAD application
    acad = Autocad(create_if_not_exists=True)

    # Define the coordinates for the pyramid's vertices
    half_length = base_length / 2
    half_width = base_width / 2

    vertices = [
        (-half_length, -half_width, 0),  # Vertex 0 (Bottom left)
        (half_length, -half_width, 0),   # Vertex 1 (Bottom right)
        (half_length, half_width, 0),    # Vertex 2 (Top right)
        (-half_length, half_width, 0),   # Vertex 3 (Top left)
        (0, 0, height)                    # Vertex 4 (Apex)
    ]

    # Define each edge of the pyramid as a polyline
    edges = [
        (0, 1),  # Edge from Vertex 0 to Vertex 1
        (1, 2),  # Edge from Vertex 1 to Vertex 2
        (2, 3),  # Edge from Vertex 2 to Vertex 3
        (3, 0),  # Edge from Vertex 3 to Vertex 0
        (0, 4),  # Edge from Vertex 0 to Apex
        (1, 4),  # Edge from Vertex 1 to Apex
        (2, 4),  # Edge from Vertex 2 to Apex
        (3, 4)   # Edge from Vertex 3 to Apex
    ]

    # Create polylines for each edge
    for edge in edges:
        start_vertex = vertices[edge[0]]
        end_vertex = vertices[edge[1]]
        polilinea = aDouble(start_vertex[0], start_vertex[1], start_vertex[2],
                            end_vertex[0], end_vertex[1], end_vertex[2])
        acad.model.Add3Dpoly(polilinea)

    print("Pyramid created with quadrangular base.")

# Example usage
create_pyramid(base_length=2, base_width=2, height=3)  # You can change the dimensions as needed