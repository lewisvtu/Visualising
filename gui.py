from Tkinter import *
import utils
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg", warn=False)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler

class MainWindow(object):
    '''Main window object for aplication'''
    def __init__(self, master):
        self.master = master
        self.data_store = []
        self.fig = plt.figure()
        self.fig.add_subplot(111, projection="3d")
        self.draw()

    def draw(self):
        self.main_frame = Frame(self.master)
        quit_b = Button(self.main_frame, text="Quit >:(", command=self._quit)
        quit_b.grid(row=0, column=0)
        data_b = Button(self.main_frame, text="Open Data Set", command=self.open_data)
        data_b.grid(row=0, column=1)
        graph_b = Button(self.main_frame, text="Open Graph Viewer", command=self.open_graph)
        graph_b.grid(row=0, column=2)
        self.main_frame.pack()

    def open_graph(self):
        window = Toplevel(self.master)
        self.graph_window = GraphWindow(window, self.fig)

    def open_data(self):
        window = Toplevel(self.master)
        self.data_window = DataWindow(window, self.data_store)
    
    def _quit(self):
        self.master.quit()
        self.master.destroy()

class GraphWindow(object):
    '''Window for drawing graphs'''
    def __init__(self, master, fig):
        self.master = master
        self.draw()
        self.set_graph(self, fig, des=False)

    def draw(self):
        close_b = Button(self.master, text="Close Window",
                         command=self.close_window).grid(row=0, column=0)
        draw_b = Button(self.master, text="Draw Graph",
                        command=self.draw_graph().grid(row=0, column=3)
        draw_b
        Label(self.master, text="Flight File Name: ").grid(row=0,column=1)
        self.fname_e = Entry(self.master).grid(row=0, column=2)

    def close_window(self):
        self.master.destroy()

    def set_graph(self, fig, des=True):
        self.canvas = FigureCanvasTkAgg(fig)
        if des:
            self.canv_widget.destroy()
        return canvas.get_tk_widget().grid(row=3)
        

class DataWindow(object):
    '''Data input table window'''
    def __init__(self, master, data_store):
        self.master = master
        self.data_store = data_store
        self.data_entries = []
        self.frame = Frame(self.master)
        close_b = Button(self.frame, text="Close Window", command=self.close_window)
        close_b.grid(row=0, column=0)
        read_b = Button(self.frame, text="Read Data", command=self.read_entry_boxes)
        read_b.grid(row=0,column=1)
        add_b = Button(self.frame, text="New Gal", command=self.add_row)
        add_b.grid(row=0, column=2)
        clear_b = Button(self.frame, text="Clear All", command=self.clear_entries)
        clear_b.grid(row=0, column=3)
        self.frame.grid()
        self.draw_entry_boxes()
        
    def clear_entries(self):
        self.data_store = []
        self.draw_entry_boxes()

    def add_row(self):
        self.read_entry_boxes()
        self.data_store.append([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
        self.draw_entry_boxes()

    def draw_entry_boxes(self):
        #Destroy all old widgets
        for entry_row in self.data_entries:
            for entry in entry_row:
                entry.destroy()
        #Remove all old widgets from list
        self.data_entries = []
        for row_no, data_row in enumerate(self.data_store):
            lab = Label(self.master, text="Galaxy no: %2u"% row_no)
            lab.grid(row=row_no + 1, column = 0)
            entry_row = [lab]
            for col_no, data_piece in enumerate(data_row):
                entry = Entry(master=self.master, width=5)
                entry.grid(row=row_no + 1, column=col_no + 1)
                entry_row.append(entry)
                entry.insert(0, data_piece)
            self.data_entries.append(entry_row)

    def read_entry_boxes(self):
        #Store data from each entry box
        for row_no, row in enumerate(self.data_entries):
            for col_no, data_piece in enumerate(row[1:]): #row[0] is a label
                self.data_store[row_no][col_no] = float(data_piece.get())



    def close_window(self):
        self.master.destroy()

rt = Tk()
w = MainWindow(rt)
rt.mainloop()
