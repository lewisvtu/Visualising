from utils import plot_from_file
import matplotlib.pyplot as plt
fname = "Paths\weave_2.txt"

def run():
    plt = plot_from_file(fname)
    plt.show()
if __name__ == "__main__":
    f = plot_from_file(fname)
    plt.show()