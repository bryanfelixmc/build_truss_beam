import numpy as np
from pyautocad import Autocad, aDouble

def parallelepiped_vertices(A, B, width, height):
    """
    Generates the vertices of a parallelepiped given vectors A and B,
    with predefined width and height. The parallelepiped is oriented
    in the direction of vector B - A.

    Args:
        A (list or tuple): Starting point.
        B (list or tuple): Ending point.
        width (float): Width of the parallelepiped.
        height (float): Height of the parallelepiped.

    Returns:
        list: List of 8 vertices of the parallelepiped.
    """
    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)
    
    AB = B - A
    
    if np.allclose(AB, np.zeros(3)):
        return [A.tolist()]
    
    # Find a vector perpendicular to AB
    # First, try to find a perpendicular vector using cross product with x-axis
    u_raw = np.cross(AB, np.array([1, 0, 0]))
    
    # If the result is zero, try with y-axis
    if np.allclose(u_raw, np.zeros(3)):
        u_raw = np.cross(AB, np.array([0, 1, 0]))
    
    # If still zero, use z-axis (though this case should not occur if AB is non-zero)
    if np.allclose(u_raw, np.zeros(3)):
        u_raw = np.cross(AB, np.array([0, 0, 1]))
    
    u_hat = u_raw / np.linalg.norm(u_raw)
    u = u_hat * width
    
    v_raw = np.cross(AB, u_raw)
    v_hat = v_raw / np.linalg.norm(v_raw)
    v = v_hat * height
    delta = -0.5 * u - 0.5 * v
    
    vertices = [
        (A + delta).tolist(),
        (A + AB + delta).tolist(),
        (A + u + delta).tolist(),
        (A + v + delta).tolist(),
        (A + AB + u + delta).tolist(),
        (A + AB + v + delta).tolist(),
        (A + u + v + delta).tolist(),
        (A + AB + u + v + delta).tolist()
    ]
    
    return vertices

def create_parallelepiped(A, B, width, height):
    """
    Creates a parallelepiped in AutoCAD based on the given points A and B,
    with specified width and height. The parallelepiped is oriented in the
    direction of vector B - A.

    Args:
        A (list or tuple): Starting point.
        B (list or tuple): Ending point.
        width (float): Width of the parallelepiped.
        height (float): Height of the parallelepiped.
    """
    # Initialize AutoCAD application
    acad = Autocad(create_if_not_exists=True)

    # Calculate the vertices of the parallelepiped
    vertices = parallelepiped_vertices(A, B, width, height)

    if len(vertices) == 1:
        print("The parallelepiped collapses to a single point.")
        return

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
        (3, 7),  # A + v to A + u + v
        (4, 5),  # A + AB + u to A + AB + v
        (4, 7),  # A + AB + u to A + u + v
        (5, 6),  # A + AB + v to A + AB + u + v
        (6, 7)   # A + u + v to A + AB + u + v
    ]

    # Create polylines for each edge
    for edge in edges:
        start_vertex = vertices[edge[0]]
        end_vertex = vertices[edge[1]]
        polyline = aDouble(start_vertex[0], start_vertex[1], start_vertex[2],
                           end_vertex[0], end_vertex[1], end_vertex[2])
        acad.model.Add3Dpoly(polyline)

    print("Parallelepiped created with each edge as a separate polyline.")

# Example usage:
A = [0, 0, 0]
B = [5, 0, 0]
width = 2
height = 1

create_parallelepiped(A, B, width, height)
