from pyautocad import Autocad, APoint, aDouble

# Initialize AutoCAD application
acad = Autocad(create_if_not_exists=True)

polilinea = aDouble(0, 0, 0,  # Bottom face
                    1, 0, 0,
                    1, 1, 0,
                    0, 1, 0,
                    0, 0, 1,  # Top face
                    1, 0, 1,
                    1, 1, 1,
                    0, 1, 1)

acad.model.Add3Dpoly(polilinea)

