import tkinter as tk
import tkinter.ttk as ttk
from tktooltip import ToolTip

app = tk.Tk()
b = ttk.Button(app, text="Button")
b.pack()
ToolTip(b, msg="Hover info")
app.mainloop()