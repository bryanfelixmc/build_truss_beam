from pyautocad import Autocad, aDouble

def create_truss_beam(base_length, base_width, num_triangles_length, num_triangles_width):
    # Initialize AutoCAD application
    acad = Autocad(create_if_not_exists=True)

    # Define the coordinates for the quadrangular face vertices
    half_length = base_length / 2
    half_width = base_width / 2

    vertices = [
        (-half_length, -half_width, 0),  # Vertex 0 (Bottom left)
        (half_length, -half_width, 0),   # Vertex 1 (Bottom right)
        (half_length, half_width, 0),    # Vertex 2 (Top right)
        (-half_length, half_width, 0)    # Vertex 3 (Top left)
    ]

    # Calculate the spacing for the triangles
    delta_length = base_length / num_triangles_length
    delta_width = base_width / num_triangles_width

    # Create triangles
    for i in range(num_triangles_length):
        for j in range(num_triangles_width):
            # Calculate the vertices of the triangle
            p1 = (vertices[0][0] + i * delta_length, vertices[0][1] + j * delta_width, 0)  # Bottom left
            p2 = (vertices[0][0] + (i + 1) * delta_length, vertices[0][1] + j * delta_width, 0)  # Bottom right
            p3 = (vertices[0][0] + i * delta_length, vertices[0][1] + (j + 1) * delta_width, 0)  # Top left
            p4 = (vertices[0][0] + (i + 1) * delta_length, vertices[0][1] + (j + 1) * delta_width, 0)  # Top right

            # Create two triangles for each rectangle
            triangles = [
                (p1, p2, p3),  # Triangle 1
                (p2, p4, p3)   # Triangle 2
            ]

            # Create polylines for each triangle
            for triangle in triangles:
                for k in range(3):
                    start_vertex = triangle[k]
                    end_vertex = triangle[(k + 1) % 3]
                    polilinea = aDouble(start_vertex[0], start_vertex[1], start_vertex[2],
                                        end_vertex[0], end_vertex[1], end_vertex[2])
                    acad.model.Add3Dpoly(polilinea)

    print("Truss beam created for the quadrangular face.")

# Example usage
create_truss_beam(base_length=8, base_width=2, num_triangles_length=3, num_triangles_width=1)