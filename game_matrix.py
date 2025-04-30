def create_game_matrix(num_vertical_axis, num_horizontal_axis, values):
    # Initialize the main matrix with zeros
    main_matrix = [[0 for _ in range(num_horizontal_axis)] for _ in range(num_vertical_axis)]
    
    # Fill the matrix with the provided values
    for (row, col, value) in values:
        main_matrix[row][col] = value

    # Generate game_columns and game_beams
    game_columns = []
    game_beams = []

    # Collecting game_columns
    for col in range(num_horizontal_axis):
        for row in range(num_vertical_axis):
            if main_matrix[row][col] == 1:
                game_columns.append([row + 1, col + 1])  # 1-indexed

    # Collecting game_beams
    # Check horizontal neighbors (left to right)
    for row in range(num_vertical_axis):
        for col in range(num_horizontal_axis - 1):
            if main_matrix[row][col] == 1 and main_matrix[row][col + 1] == 1:
                beam = [[row + 1, col + 1], [row + 1, col + 2]]
                game_beams.append(beam)

    # Check vertical neighbors (top to bottom)
    for col in range(num_horizontal_axis):
        for row in range(num_vertical_axis - 1):
            if main_matrix[row][col] == 1 and main_matrix[row + 1][col] == 1:
                beam = [[row + 1, col + 1], [row + 2, col + 1]]
                game_beams.append(beam)

    return game_columns, game_beams, main_matrix

# Example usage
if __name__ == "__main__":
    num_vertical_axis = 4
    num_horizontal_axis = 3
    values = [
        (0, 0, 1),  
        (0, 1, 1),  
        (0, 2, 1),  
        (1, 1, 1),  
        (2, 1, 1),  
        (3, 1, 1),  
    ]


    game_columns, game_beams, main_matrix = create_game_matrix(num_vertical_axis, num_horizontal_axis, values)
    
    print("main_matrix:", main_matrix)
    print("Game Columns:", game_columns)
    print("Game Beams:", game_beams)