# Gauss Elimination using pivoting and scaling
TOLERANCE = 5
significance = (10 ** (-1 * TOLERANCE))


# This function is for rounding a float using tolerance
def my_round(some_float):
    p = float(10 ** TOLERANCE)
    return int(some_float * p + 0.5) / p


# This function prints the matrix and the vector
def print_matrix(arr):
    Result = ""
    n = len(arr)
    for i in range(n):
        for j in range(n + 1):
            if j == 0:
                Result += "||  "
            x = "{:0." + str(TOLERANCE) + "f}"
            Result += x.format(arr[i][j]) + "  "
            if j == n - 1:
                Result += " | "
            if j == n:
                Result += "||\n"
    return Result


# This function logs the scaling procedure steps
def log_scaling(arr, row):
    Result = "\nThe biggest coefficient after scaling is the row :" + str(row + 1) + "\nwhich is :" + str(arr[row])
    return Result


''' Scaling procedure
This Function takes the whole array and scale each row to make the largest coefficient equals to 1
And it returns an integer with the biggest value of the column i so we can make this row the first one'''


def Scaling(arr, column, steps):
    n = len(arr)
    column_values = list()
    for i in range(column, n):
        if abs(max(arr[i][column + 1:-1], key=abs)) < significance:
            steps.append(log_scaling(arr, i))
            return i
        column_values.append((arr[i][column] / abs(max(arr[i][column + 1:-1], key=abs))))
    result = column_values.index(max(column_values, key=abs)) + column
    steps.append(log_scaling(arr, result))
    return result


def log_pivoting(arr):
    Result = "\nThe Matrix After being pivoted is :\n"
    Result += print_matrix(arr)
    return Result


''' Partial Pivoting procedure
This Function Takes the Array and two rows and it swaps them after choosing the biggest coefficient after being scaled'''


def Pivoting(arr, i, row, steps):
    temp = arr[i]
    arr[i] = arr[row]
    arr[row] = temp
    steps.append(log_pivoting(arr))


def log_elimination(arr, i):
    Result = "Elimination steps :\n"
    for x in range(i + 1, len(arr)):
        Result += "Row " + str(x + 1) + " = Row " + str(x + 1) + " - m" + str(x + 1) + str(i + 1) + " Row" + str(
            i + 1) + "\n"
    Result += print_matrix(arr)
    return Result


# Elimination procedure
'''This Function takes the array and it calls scaling and pivoting procedures and it uses the eliminations steps
and finally it process the backward substitution procedure and returns the result as a vector'''


def Gauss_Elimination(matrix, vector):
    n = len(matrix)

    arr = list()  # Augmented matrix

    for x in range(n):
        arr.append(list())
        for y in range(n):
            arr[x].append(matrix[x][y])
            if y == n - 1:
                arr[x].append(vector[x])

    steps = list()  # Logs
    steps.append("The initial matrix is :\n\n" + print_matrix(arr))

    # Checking for zero columns
    for y in range(n):
        for x in range(n):
            if arr[x][y] > significance:
                break
            if x == n - 1:
                steps.append("The matrix is singular")
                return steps

    # Elimination steps
    for i in range(n - 1):

        # Checking for zero rows
        for x in range(i, n):
            for y in range(i, n):
                if abs(arr[x][y]) > significance:
                    break
                if y == n - 1:
                    steps.append("The matrix is singular")
                    return steps

        # Scaling
        row = Scaling(arr, i, steps)  # Scaling
        # Pivoting
        if row != i:
            Pivoting(arr, i, row, steps)  # Pivoting
        else:
            steps.append("\nPivoting is not required\n")

        # Elimination Steps
        for x in range(i + 1, n):
            multiplier = arr[x][i] / arr[i][i]
            for y in range(i, n + 1):
                arr[x][y] = arr[x][y] - multiplier * arr[i][y]
                if abs(arr[x][y]) < significance:
                    arr[x][y] = 0

        steps.append(log_elimination(arr, i))

    # Checking if the matrix is singular
    if abs(arr[n - 1][n - 1]) < significance:
        if abs(arr[n - 1][n]) < significance:
            steps.append("The matrix is singular")
            return steps

    # Backward substitution
    Result = [0 for _ in range(n)]
    Result[n - 1] = arr[n - 1][n] / arr[n - 1][n - 1]
    Result[n - 1] = my_round(Result[n - 1])
    for i in range(n - 2, -1, -1):
        Result[i] = arr[i][n]
        for j in range(i + 1, n):
            Result[i] -= (arr[i][j] * Result[j])
        Result[i] /= arr[i][i]
        if abs(Result[i]) < significance:
            Result[i] = 0
        else:
            if Result[i] < 0:
                Result[i] = -1 * Result[i]
                Result[i] = my_round(Result[i]) * -1
            else:
                Result[i] = my_round(Result[i])
    steps.append("The Result is :\n" + str(Result))
    return steps


n = int(input("Enter matrix dimensions : "))
print("Enter the matrix and the vector next it like this :\n1 2 3\n4 5 6\nWhere v = [3,6] and the matrix is 2x2\n")
matrix = list()
vector = list()
for i in range(n):
    line = input()
    matrix.append([int(x) for x in line.split()[:-1]])
    vector.append(int(line.split()[-1]))


for x in Gauss_Elimination(matrix, vector):
    print(x)
