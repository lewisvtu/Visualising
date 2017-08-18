from Tkinter import *
import utils
import matplotlib
matplotlib.use("TkAgg", warn=False)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler

class MainWindow(object):
    '''Main window object for aplication'''
    def __init__(self, master):
        self.master = master
        self.main_frame = Frame(self.master)
        self.q_button = Button(self.main_frame, text="Quit >:(", command=self._quit)
        self.q_button.grid(row=3, column=0)
        self.d_button = Button(self.main_frame, text="Open Data Set", command=self.open_data)
        self.d_button.grid(row=1, column=0)
        self.main_frame.pack()
        self.data_store = [
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [1,2,3,4,5,6,7,8,9,0,1,2,3,4]
        ]


    def open_data(self):
        window = Toplevel(self.master)
        self.data_window = DataWindow(window, self.data_store)
    
    def _quit(self):
        self.master.quit()
        self.master.destroy()

class DataWindow(object):
    '''Data input table window'''
    def __init__(self, master, data_store):
        self.master = master
        self.data_store = data_store
        self.data_entries = []
        self.frame = Frame(self.master)
        self.c_button = Button(self.frame, text="Close Window", command=self.close_window)
        self.c_button.grid()
        self.frame.grid()
        self.draw_entry_boxes()
        
    def draw_entry_boxes(self):
        for row_no, data_row in enumerate(self.data_store):
            lab = Label(self.master, text="Galaxy no: %2u"% row_no)
            lab.grid(row=row_no + 1, column = 0)
            entry_row = []
            for col_no, data_piece in enumerate(data_row):
                entry = Entry(master=self.master)
                entry.grid(row=row_no + 1, column=col_no + 1)
                entry_row.append(entry)
                entry.insert(0, data_piece)
            self.data_entries.append(entry_row)



  

    def close_window(self):
        self.master.destroy()

rt = Tk()
w = MainWindow(rt)
rt.mainloop()
