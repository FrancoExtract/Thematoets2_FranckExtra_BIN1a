# Thematoets 2
# Franck Extra

try:
    import re
    import matplotlib.pyplot as plt
    import tkinter as tk
except ModuleNotFoundError:
    print("Matplotlib module not found")  # It can't find the module for some reason...


def read_file():
    gff_file = open(gff, "r")
    # print(gff.readable())

    gff_file.close()


# def graph(x, y):
# plt.plot(x, y, 'g-')
# plt.title("Total counts of kinase")
# plt.xlabel("Chromosome")
# plt.ylabel("Kinase count")
# plt.show()


if __name__ == '__main__':
    gff = "GCF_000013425.1_ASM1342v1_genomic.gff"

# Prototype graph
# graph(([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), ([2, 4, 6, 8, 10, 12, 14, 16, 18, 20]))
