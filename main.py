import itertools
def sat(x1, x2):
    return (x1 or not x2) and (not x1 or x2)
for (x1, x2) in list(itertools.product([False, True], repeat=2)):
    print(x1, x2, '-->', sat(x1, x2))
