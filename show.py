import numpy as np
import tkinter as tk
import tkinter.ttk as ttk

win = tk.Tk()
win.size = np.array([450,450])

canvas = tk.Canvas(win, width=win.size[0], height=win.size[1], highlightthickness=0)

def circle(self, color, selected, *bbox):
    if selected:
        outline = '#004500'
        width = 3
    else:
        outline = 'black'
        width = 1
    self.canvas.create_oval(*tuple(int(x) for x in bbox), fill=color,
                            outline=outline, width=width)

from ase.cli.main import old
old()
