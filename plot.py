from utils import plot_from_file
import matplotlib.pyplot as plt
fname = "galaxy_tour2_.txt"

def run():
    plt = plot_from_file(fname)
    plt.show()
if __name__ == "__main__":
    f = plot_from_file(fname)
    plt.show()