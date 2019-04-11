from SubProb import *

if __name__ == "__main__":
    inf = float('inf')
    A = [[inf, inf, inf, inf],
        [3, inf, 6, 14],
        [5, inf, inf, 6],
        [9, inf, 5, inf]]



    print(A, flush=True)

    a_prob = subProblem(A, 0, 0)
    a_prob.getLowerBound(0, 1)

    print(a_prob, flush=True)
    print(a_prob)
    print(cost)