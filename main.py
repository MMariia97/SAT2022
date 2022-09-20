import itertools
import argparse

def sat(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16):
    return not x1 and not x2 and not x3 and not x4 and not x5 and x6 and not x7 and x8 and x9 and not x10 and x11 and x12 and x13 and x14 and x15 and not x16

def main():
    for (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16) in list(itertools.product([False, True], repeat=16)):
        if (sat(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16) == True):
            print(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, '-->', True)

if __name__ == '__main__':
    main()
