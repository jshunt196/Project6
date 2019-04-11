import copy

class subProblem:

    def __init__(self, myMatrix, myCost, myCity):
        # reduce it's matrix and expand it.
        self.cost = myCost
        self.lowerBound = 0
        self.path = []
        self.matrix = myMatrix
        self.cityFrom = myCity # index of city

    def __lt__(self, other):
        if self.cost < other.cost: # return lower cost
            return self
        elif self.cost > other.cost:
            return other
        else: # if tie return smaller lower bound
            if self.lowerBound < other.lowerBound:
                return self
            elif self.lowerBound > other.lowerBound:
                return other
            else: # if tie return bigger path 
                if len(self.path) >= other.lowerBound:
                    return self
                else:
                    return other

    def reduce(self):
        self.getLowerBound(self.path[-2], self.path[-1])

    #  TIME COMPLEXITY: O(n^2), two for loops with size n, calling function with time n
    # SPACE COMPLEXITY: O(n^2), uses the matrix of the cities, n x n
    def getLowerBound(self, fromR, toC):
        matrixCost = 0

        # reduce all rows in the matrix, skips city coming from
        for i in range(len(self.matrix)):
            if i is not fromR:
                myVar = self.checkRow(i)
                if myVar != float("inf"):
                    matrixCost += myVar

        # reduce all columns in the matrix, skips city coming from
        for i in range(len(self.matrix)):
            if i is not toC:
                myVar = self.checkColumn(i)
                if myVar != float("inf"):
                    matrixCost += myVar

        self.lowerBound += matrixCost

    #  TIME COMPLEXITY: O(n), technically n+n because of two for loops size n
    # SPACE COMPLEXITY: O(n^2), uses the matrix of the cities, n x n
    def checkColumn(self, colI):
        cost = 0

        # get smallest num
        smallest = float("inf")
        for i in range(len(self.matrix)):
            if self.matrix[i][colI] < smallest:
                smallest = self.matrix[i][colI]

        # check if whole column is inf
        if smallest == float("inf"):
            return 0
        
        # subtract smallest from each value in column
        for i in range(len(self.matrix)):
            if self.matrix[i][colI] != float("inf"):
                self.matrix[i][colI] -= smallest
        
        cost += smallest
        return cost

    #  TIME COMPLEXITY: O(n), technically n+n because of two for loops size n
    # SPACE COMPLEXITY: O(n^2), uses the matrix of the cities, n x n
    def checkRow(self, rowI):
        cost = 0

        # get smallest num
        smallest = float('inf')
        for i in range(len(self.matrix)):
            if self.matrix[rowI][i] < smallest:
                smallest = self.matrix[rowI][i]

        # check if whole row is inf
        if smallest == float('inf'):
            return 0

        # subtract smallest from each value in row
        for i in range(len(self.matrix)):
            if self.matrix[rowI][i] != float('inf'):
                self.matrix[rowI][i] -= smallest

        cost += smallest
        return cost

    # will initialize but not reduce matrix from 'A' to city 'B'
    #  TIME COMPLEXITY: O(n), for loop is length of num cities
    # SPACE COMPLEXITY: O(n^2), uses the matrix of the cities, n x n
    def initMatrix(self, problemMatrix, parent, toC):
        fromR = self.path[-1]
        self.cost += problemMatrix[fromR][toC]
        # lowerBound = parent lowerbound + cost from A to B + reduction of current (done later)
        self.lowerBound += problemMatrix[fromR][toC] + parent.lowerBound

        tempMatrix = problemMatrix
        
        # infinity out row 'A' and column 'B'
        for i in range(len(problemMatrix)):
            tempMatrix[fromR][i] = float("inf")
            tempMatrix[i][toC]   = float("inf")

        # infinity out 'B' to 'A'
        tempMatrix[toC][fromR] = float("inf")

        return tempMatrix

    def __repr__(self):
        return self.matrix.__repr__()