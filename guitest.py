import Tkinter as tk
import utils
import matplotlib
matplotlib.use("TkAgg", warn=False)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler


class App(object):
    def __init__(self, master):
        print "creating obj"
        self.master = master
        self.inputs = {}
        self.tab_frame = tk.Frame(master)
        self.tab_frame.pack(side=tk.LEFT)
        self.but_frame = tk.Frame(master)
        self.but_frame.pack(side=tk.RIGHT)
        create = tk.Button(self.but_frame, text="New Galaxy", command=self.create_entry_row)
        create.grid(column=0, row=0)
        quit_button = tk.Button(self.but_frame, text="quit", command=self._quit)
        quit_button.grid(column=0,row=6)
        read = tk.Button(self.but_frame, text="Read Gals", command=self.read_rows)
        read.grid(column=0, row=1)
        gb = tk.Button(self.but_frame, text="Draw Graph", command=self.plot_graph)
        self.fne = tk.Entry(self.but_frame)
        self.fne.grid(column=0,row=3)
        gb.grid(column=1,row=3)

    def _quit(self):
        self.master.quit()
        self.master.destroy()

    def create_entry_row(self):
        if len(self.inputs) > 0:
            max_row = max(self.inputs.keys())
        else:
            max_row = 0
        self.create_entries(self.tab_frame, max_row + 1, 14)

    def create_entries(self, root, row, cols):
        lab = tk.Label(root, text="Galaxy_no: %s" % row)
        lab.grid(column=0, row=row)
        entries = []
        self.inputs[row] = entries
        for col in range(cols):
            entry = tk.Entry(root)
            entry.grid(column=col + 1, row=row)
            entries.append(entry)
        
    def read_rows(self):
        for key in self.inputs.keys():
            entries = self.inputs[key]
            args = [float(entry.get()) for entry in entries if entry.get() != ""]
            print args

    def plot_graph(self):
        f_name = self.fne.get()
        f = utils.plot_from_file(f_name)
        canvas = FigureCanvasTkAgg(f, self.master)
        ax = f.gca()
        ax.mouse_init()
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH)



rt = tk.Tk()
app = App(rt)
rt.mainloop()
