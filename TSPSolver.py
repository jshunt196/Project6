#!/usr/bin/python3

import copy
import itertools
import heapq
from SubProb import *
from TSPClasses import *
import numpy as np
import time
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None
        self.totalStates = 1
        self.prunedStates = 0

    def setupWithScenario(self, scenario):
        self._scenario = scenario

    ''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution,
		time spent to find solution, number of permutations tried during search, the
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time()-start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    # b, number of subproblems a subproblem can generate when expanded
    # n, number of cities
    #  TIME COMPLEXITY: O(n^2*b^2), b^2 comes from two loops, a while loop and a for
    #		loop. Each are both size b, representing the subproblems that could be
    # 		produced from subproblems in the priority queue. The n^2 is generated from
    #		the function calls that are made. Calling expand and reduce are both in
    #		n^2 time.
    # SPACE COMPLEXITY: O(n^2*b^2), uses the matrix of the cities, n x n. b is used
    #		as a length of measurement for each problem in the queue. Another b could
    #		be argued to be used to measure the length of the priorityQ itself.
    def branchAndBound(self, time_allowance=60.0):
        start_time = time.time()
        cities = self._scenario.getCities()
        matrix = self.makeMatrix(cities)

        # make first problem
        firstProblem = subProblem(matrix, 0, 0)
        firstProblem.path.append(0)
        firstProblem.cityFrom = 0

        # initialize priorityQ
        priorityQ = [firstProblem]
        heapq.heapify(priorityQ)

        # initialize other variables
        randomResults = self.defaultRandomTour()
        bestSoFar = randomResults['cost']
        bestPath = []
        counter = 1
        queueSizeMax = 0
        numBSFChanges = 0

        # b*b*n*n
        # b*b?*n^2
        while len(priorityQ) is not 0 and time.time() - start_time < time_allowance:
            if len(priorityQ) > queueSizeMax:  # get max size
                queueSizeMax = len(priorityQ)
            problem = heapq.heappop(priorityQ)
            listProblems = self.expand(problem)
            # b
            for i in listProblems:  # examine each expanded problem
                if i.cost < bestSoFar:  # skips reduce(), which is expensive
                    i.reduce()
                    if len(i.path) is len(matrix):  # if 'i' is a tour
                        tour = []
                        tourPath = []
                        for j in i.path:
                            tour.append(cities[j])
                            tourPath.append(cities[j]._index)

                        if TSPSolution(tour).cost < bestSoFar:
                            bestSoFar = TSPSolution(tour).cost
                            bestPath = tourPath
                            numBSFChanges += 1
                    elif i.lowerBound < bestSoFar:  # chance of finding a smaller solution
                        heapq.heappush(priorityQ, i)
                        # priorityQ.append(i)
                    else:  # if problem lower bound is bigger than BSSF
                        self.prunedStates += 1
                else:  # if problem cost is bigger than BSSF
                    self.prunedStates += 1

        finalPath = []
        for i in bestPath:  # get path of cities
            finalPath.append(cities[i])

        end_time = time.time()
        # print("Best: " + str(bestSoFar))
        results = {}
        results['cost'] = TSPSolution(finalPath).cost
        results['time'] = end_time - start_time
        results['count'] = numBSFChanges
        results['soln'] = TSPSolution(finalPath)
        results['max'] = queueSizeMax
        results['total'] = self.totalStates
        results['pruned'] = self.prunedStates
        return results

    # expand given problem into more problems
    #  TIME COMPLEXITY: O(n^2), n represents number of cities in the graph. Given a problem's
    #		matrix, expand will create new subproblems, going through each city to see if it
    #		has been added to the path; if not, then it creates a problem with that city.
    #		for loop is size n. newProblem.initMatrix() is also size n time.
    # SPACE COMPLEXITY: O(b), a new subproblem is created. This means that the problem has a
    #		set of more problems that can be expanded. b is simplified, also including the
    #		n x n matrix it contains.
    def expand(self, problem):
        listProblems = []
        currentPath = problem.path
        problemMatrix = copy.deepcopy(problem.matrix)
        for i in range(len(problem.matrix)):
            if i in problem.path:
                continue
            else:
                # create a new subProblem
                newProblem = subProblem(copy.deepcopy(
                    problemMatrix), copy.deepcopy(problem.cost), i)
                newProblem.path = currentPath.copy()
                newProblem.matrix = newProblem.initMatrix(
                    copy.deepcopy(problemMatrix), problem, i)
                newProblem.cityFrom = i
                newProblem.path.append(i)
                self.totalStates += 1

                listProblems.append(newProblem)
        return listProblems

    ''' <summary>
		This is the entry point for the greedy solver, which you must implement for
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

    def greedy(self, time_allowance=60.0):
        results = {}
        start_time = time.time()
        cities = self._scenario.getCities()
        matrix = self.makeMatrix(cities)
        tours = []

        for i in range(len(cities)):
            path = []
            path.append(cities[i])
            current = cities[i]
            while len(path) != len(cities):  # iterate through each node
                smallest = float("inf")
                newCity = None
                for j in cities:  # find smallest path
                    if j not in path and matrix[current._index][j._index] < smallest:
                        smallest = matrix[current._index][j._index]
                        newCity = j
                if newCity is None:  # if there is no possible connection
                    break  # specific tour should return inf cost
                path.append(newCity)
                current = newCity
            tours.append(path)

        finalPath = None
        smallestTour = float("inf")
        for i in range(len(tours)):
            if len(tours[i]) == len(cities) and TSPSolution(tours[i]).cost < smallestTour:
                smallestTour = TSPSolution(tours[i]).cost
                finalPath = tours[i]

        if finalPath is None:
            return None

        end_time = time.time()
        results['cost'] = TSPSolution(finalPath).cost
        results['time'] = end_time - start_time
        results['count'] = len(finalPath)
        results['soln'] = TSPSolution(finalPath)
        results['max'] = 0
        results['total'] = 0
        results['pruned'] = 0
        return results

    def makeMatrix(self, cities):  # make the first matrix for the set of cities
        matrix = []
        for i in range(len(cities)):
            matrix.append([])
            for j in range(len(cities)):
                matrix[i].append(cities[i].costTo(cities[j]))

        return matrix

    ''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

    def fancy(self, time_allowance=60.0):
        start_time = time.time()
        cities = self._scenario.getCities()
        self.matrix = np.zeros((len(cities), len(cities)))
        for i in range(len(cities)):
            for j in range(len(cities)):
                self.matrix[i, j] = cities[i].costTo(cities[j])
        self.inf_matrix = self.matrix.copy()
        array1 = []
        self.map = {}
        for i in range(len(cities)):
            path = [cities[i]]
            array1.append(path)
            self.map[i] = path
        while len(array1) > 1:
            array1 = self.combine(array1)
        end_time = time.time()

        results = {}
        results['cost'] = TSPSolution(array1[0]).cost
        results['time'] = end_time - start_time
        results['count'] = 0
        results['soln'] = TSPSolution(array1[0])
        results['max'] = 0
        results['total'] = 0
        results['pruned'] = 0
        return results

    def combine(self, array1):
        array2 = []
        temp_matrix = self.inf_matrix.copy()
        while len(array1) > 1:
            city1, city2 = self.minEdgeInMatrix(temp_matrix)
            path1, path2 = self.findPaths(city1, city2)
            if not(path1 in array1 and path2 in array1) or path1 is path2:
                temp_matrix[city1][city2] = np.inf
                continue
            bestCost = (np.inf, [])
            len1 = len(path1)
            len2 = len(path2)
            for i in range(len(path1)):
                nexti = (i+1) % len1
                for j in range(len(path2)):
                    nextj = (j+1) % len2
                    cost = (self.matrix[path1[i]._index][path2[nextj]._index] + self.matrix[path1[nexti]._index][path2[j]._index] - self.matrix[path1[i]._index][path1[nexti]._index] - self.matrix[path2[j]._index][path2[nextj]._index])
                    if cost < bestCost[0]:
                        bestCost = (cost, (i, nexti, j, nextj))
            if bestCost[0] is np.inf:
                temp_matrix[city1][city2] = np.inf
                continue
            new_path = path1[:bestCost[1][0]]
            new_path += path2[bestCost[1][3]:]
            new_path += path2[:bestCost[1][3]]
            new_path += path1[bestCost[1][0]:]
            array2.append(new_path)
            array1.remove(path1)
            array1.remove(path2)
            for city in new_path:
                self.map[city._index] = new_path
            for city1 in new_path:
                for city2 in new_path:
                    self.inf_matrix[city1._index][city2._index] = np.inf
        if len(array1) is 1:
            array2 += array1
            return array2
        return array2

    def minEdgeInMatrix(self, matrix):
        ind = np.unravel_index(
            np.argmin(matrix), matrix.shape)
        city1 = ind[0]
        city2 = ind[1]
        return city1, city2
        #	https://docs.scipy.org/doc/numpy/reference/generated/numpy.argmin.html

    # take two cities,
    def findPaths(self, city1, city2):
        return self.map[city1], self.map[city2]

    ''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''
