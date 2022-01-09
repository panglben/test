import numpy as np
import tkinter as tk
import tkinter.ttk as ttk

class GUIWindow:
    def __init__(self, windowname):
        self.win = tk.Tk()
        self.win.title(windowname)
        self.win.size = np.array([800,800])

        self.canvas = tk.Canvas(self.win, bg='white')

    def circle(self, color, selected, *bbox):
        if selected:
            outline = '#004500'
            width = 3
        else:
            outline = 'black'
            width = 1
        self.canvas.create_oval(*tuple(int(x) for x in bbox), fill=color,
                                outline=outline, width=width)


gui = GUIWindow('test')
gui.circle("#ffffff", True, 10, 20, 10, 20)

gui.canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

gui.win.mainloop()