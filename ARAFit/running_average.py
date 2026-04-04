#!/usr/bin/env python3

import sys
import numpy as np

def running_average(data, window_size):
    kernel = np.ones(window_size) / window_size
    return np.convolve(data, kernel, mode='valid')

def main():
    if len(sys.argv) < 3:
        print("Usage: python running_avg.py <input_file.csv> <window_size>")
        sys.exit(1)

    filename = sys.argv[1]
    window_size = int(sys.argv[2])

    # Load two-column CSV data
    data = np.loadtxt(filename, delimiter=',')

    x = data[:, 0]
    y = data[:, 1]

    # Apply running average to y
    y_avg = running_average(y, window_size)

    # Adjust x to match the shorter y_avg length
    x_avg = x[window_size - 1:]

    # Output result
    for xi, yi in zip(x_avg, y_avg):
        print(f"{xi},{yi}")

if __name__ == "__main__":
    main()
