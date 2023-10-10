class Matrix:
    # Z-row is row 0
    # First are main variables, then slack variables
    # n-1: Ratio, n-2: Right-hand
    def __init__(self, rows, cols):
        self.rows = rows
        self.cols = cols
        self.data = [[0 for _ in range(cols)] for __ in range(rows)]

    def set(self, row, col, value):
        if 0 <= row < self.rows and 0 <= col < self.cols:
            self.data[row][col] = value
        else:
            raise ValueError("Invalid indices")

    def get(self, row, col):
        if 0 <= row < self.rows and 0 <= col < self.cols:
            return self.data[row][col]
        else:
            raise ValueError("Invalid indices")

    def set_data_from_list(self, input_list):
        if len(input_list) != self.rows * self.cols:
            raise ValueError("Input list size does not match matrix dimensions")

        for i in range(self.rows):
            for j in range(self.cols):
                self.set(j, i, input_list[i * self.cols + j])

    def get_rows(self):
        return self.rows

    def get_cols(self):
        return self.cols

    # Get the minimum from a specific row, use this to find the new entering variable
    def min_from_row(self, row):
        min_value = self.get(row, 0)
        min_position = 0
        for col in range(self.cols):
            if self.get(row, col) < min_value:
                min_value = self.get(row, col)
                min_position = col
        return min_value, min_position

    # Get the minimum from a specific column.
    def min_from_col(self, col, negative_flag=True):
        min_value = 100000
        min_position = 1
        for row in range(self.rows - 1):
            if self.get(row+1, col) <= 0 and not negative_flag:
                continue
            if self.get(row+1, col) < min_value:
                min_value = self.get(row+1, col)
                min_position = row+1
        return min_value, min_position

    def subtracting_rows(self, first_row, second_row, coefficient=1, last_col=0):
        for col in range(self.cols + last_col):
            new_value = self.get(first_row, col) - self.get(second_row, col) * coefficient
            self.set(first_row, col, new_value)

    def diagonal_null_check(self):
        for i in range(self.rows-1):
            for j in range(self.rows-1):
                if self.get(i, j) != 0 and i != j:
                    return False
        return True

    def z_row_zero(self):
        for i in range(self.cols-1):
            if self.get(0, i) < 0:
                return True
        return False
    
    def show_matrix(self):
        for i in range(self.rows):
            for j in range(self.cols):
                print(self.get(i, j), end='\t')
            print(end = '\n')
        print(end = '\n')
        print(end = '\n')
    def __str__(self):
        return str(self.data)


def lpp_matrix_generator(coefficient_c, coefficient_a, coefficient_right_hand, nb_variables, nb_equations):
    matrix = Matrix(nb_equations + 1, nb_variables + nb_equations + 2)

    # Заполняю коэффициентами целевой функции
    for i in range(nb_variables):
        matrix.set(0, i, -1 * coefficient_c[i])

    # Заполняю коэффициентами огрничений
    for row in range(nb_equations):
        for col in range(nb_variables):
            matrix.set(row+1, col, coefficient_a.get(row, col))

    
    # Заполняю добавочными коэффициентами
    for row in range(nb_equations):
            matrix.set(row+1, row+nb_variables, 1)

    # Заполняю коэффициентами right-hand
    for row in range(nb_equations):
        matrix.set(row+1, nb_equations+nb_variables, coefficient_right_hand[row])
    return matrix


def optimization_step(matrix, nb_variables, nb_equations, solution):
    resolution_column = matrix.min_from_row(0)[1]
    #Entering variable
    
    for row in range(nb_equations):
        #Check whether ratio is okay to calculate
        if (matrix.get(row+1,resolution_column) * matrix.get(row+1,matrix.get_cols()-2) >= 0):
            matrix.set( row+1, matrix.get_cols()-1, matrix.get(row+1, matrix.get_cols()-2) / matrix.get(row+1, resolution_column))

    #Exiting variable
    resolution_row = matrix.min_from_col(matrix.get_cols()-1, False)[1]
    matrix = get_null_col(matrix, resolution_row, resolution_column, nb_variables, nb_equations)
    solution[resolution_row-1] = resolution_column
    return matrix


def get_null_col(matrix, resolution_row, resolution_column, nb_variables, nb_equations):
    print("Exiting variable: ", resolution_row)
    print("Entering variable: ", resolution_column)
    resolution_value = matrix.get(resolution_row, resolution_column)
    for col in range(matrix.get_cols()-1):
        matrix.set(resolution_row,col,matrix.get(resolution_row,col)/resolution_value)
    resolution_value = 1
    for row in range(matrix.get_rows()):
        if row == resolution_row:
            continue
        coefficient = matrix.get(row, resolution_column) / resolution_value
        matrix.subtracting_rows(row, resolution_row, coefficient, -1)
    
    return matrix


def main():
    minimize = 2;
    while (minimize != 1 and minimize != 0):
        minimize = int(input("Please specify whether this is a minimization or maximization problem.\n1 - Minimzation\n0 - Maximization\n"))
    
    nb_variables = 0;
    while (nb_variables <= 0):
        nb_variables = int(input("Please specify the number of variables\n"))

    print("Please submit the coefficients for the z expression")
    coefficient_c = [0]*nb_variables
    for i in range(nb_variables):
        coefficient_c[i] = int(input())
        
    # Make it a maximization problem
    if (minimize):
        for i in range(nb_variables):
            coefficient_c[i] *= -1
    nb_equations = 0
    while (nb_equations <= 0):
        nb_equations = int(input("Please specify the number of inequations (don't count sign constraints)\n"))

    # Reading the inequations
    coefficient_a = Matrix(nb_equations, nb_variables)
    for i in range(nb_equations):
        print("Equation ", i)
        for j in range(nb_variables):
            x = int(input())
            coefficient_a.set(i,j,x)
        print(coefficient_a)
    coefficient_right_hand = [0]*nb_equations
    print("Please submit the right hand values")
    for i in range(nb_equations):
        coefficient_right_hand[i] = int(input())

    main_matrix = lpp_matrix_generator(coefficient_c, coefficient_a, coefficient_right_hand, nb_variables, nb_equations)
    solution = [0]*nb_equations
    for i in range (nb_equations):
        solution[i] = i+nb_variables

    iteration = 1
    while main_matrix.z_row_zero():
        print("Iteration #" + str(iteration))
        main_matrix = optimization_step(main_matrix, nb_variables, nb_equations, solution)
        main_matrix.show_matrix()
        iteration+=1
        
    variables = [0]*(nb_equations+nb_variables)
    for i in range(nb_equations):
        variables[solution[i]] = main_matrix.get(i+1, main_matrix.get_cols()-2)
    print("The solution for this LLP is the following:")
    if (not(minimize)):
        print("The maximum value for Z is " + str(main_matrix.get(0,main_matrix.get_cols()-2)) + " such that:")
    else:
        print("The minimum value for Z is " + str(-main_matrix.get(0,main_matrix.get_cols()-2)) + " such that:")

    for i in range(nb_variables):
        print("x" + str(i+1) + " = " + str(variables[i]), end = ", ")


if __name__ == "__main__":
    main()


