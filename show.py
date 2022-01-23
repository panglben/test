import numpy as np
import tkinter as tk
import tkinter.ttk as ttk

class GUIWindow:
    def __init__(self, windowname):
        self.win = tk.Tk()
        self.win.title(windowname)
        self.win.size = np.array([800,800])

        self.canvas = tk.Canvas(self.win, bg='white')
        
        self.canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

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
r = 20
A1 = [120, 120]
A2 = [125, 125]
gui.circle("#ffffff", False, A2[0], A2[1], A2[0]+r, A2[1]+r)
gui.circle("#ffffff", False, A1[0], A1[1], A1[0]+r, A1[1]+r)


gui.win.mainloop()