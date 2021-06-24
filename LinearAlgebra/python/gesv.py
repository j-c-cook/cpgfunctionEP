# Jack C. Cook
# July 30, 2020

import numpy as np


def main():
    A = np.array([[5, 9, 10], [2, 7, 3], [8, 2, 4]])
    b = np.array([22, 13, 17])

    x = np.linalg.solve(A, b)

    print(x)


if __name__ == '__main__':
    main()
