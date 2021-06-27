import numpy as np


def main():
    A = np.array([
        [1., 2., 3., 4., 5.],
        [6., 7., 8., 9., 10.]
    ])
    x = np.array([4., 8., 10., 6., 8.])
    x = x[:, np.newaxis]

    y = A @ x

    print(y)


if __name__ == '__main__':
    main()
