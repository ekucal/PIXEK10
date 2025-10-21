import tkinter as tk
from tkinter import ttk

# Tworzenie okna głównego
root = tk.Tk()
root.title("Multi-Select Treeview")

# Tworzenie Treeview
tree = ttk.Treeview(root, columns=("ID", "Value"), show="headings")
tree.heading("ID", text="ID")
tree.heading("Value", text="Value")

# Dodawanie danych
data = [(1, 10), (2, 20), (3, 15), (4, 30), (5, 25)]
for item in data:
    tree.insert("", "end", values=item)

tree.pack()

# Ręczne zaznaczanie wielu wierszy
selected_items = []

def on_item_click(event):
    global selected_items
    item = tree.selection()
    
    # Jeśli wiersz już był zaznaczony, usuń go; jeśli nie, dodaj
    if item in selected_items:
        selected_items.remove(item)
    else:
        selected_items.append(item)
    
    print("Zaznaczone:", [tree.item(i)["values"] for i in selected_items])

tree.bind("<ButtonRelease-1>", on_item_click)

root.mainloop()