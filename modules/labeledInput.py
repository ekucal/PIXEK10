

from tkinter import ttk, Tk
import tkinter as tk
#import tkinter.font as tkFont


"""  definition/implementation from: Python GUI Programming, A Complete Reference-Guide """

class LabeledInput(tk.Frame):
    tag=''
            
    """A widget merging a label and edit together."""
    """ page :  9% """
    """  
        mandatory :  parent,label='',input_class=ttk.Entry,
        optional  : input_var=None, input_args=None, label_args=None,
        tk.Frame options:  kwargs
        
        input_var: assign input
        
        
* parent: This argument is a reference to the parent widget; all widgets we create will take this as the first argument.
* label: This the text for the label part of the widget.
* input_class: This is the class of the widget we want to create. It should be an actual callable class object, not a string. If left blank, ttk.Entry will be used.
* input_var: This is a Tkinter variable to assign to the input. It's optional, since some widgets don't use variables.
* input_args: This is an optional dictionary of any additional arguments for the input constructor.
* label_args: This is an optional dictionary of any additional arguments for the label constructor.
* **kwargs: Finally, we catch any additional keyword arguments in **kwargs. These will be passed to the Frame constructor.
    """
    
    def __init__(self,parent,label='',input_class=ttk.Entry,
         input_var=None, input_args=None, label_args=None,
         **kwargs):
        
        super().__init__(parent,**kwargs)
        input_args=input_args or {}
        label_args=label_args or {}
        
        self.variable=input_var
        
        if input_class in (ttk.Checkbutton, ttk.Button,ttk.Radiobutton):
            input_args["text"] = label
            input_args["variable"] = input_var
        else:
            self.label = ttk.Label(self, text=label, **label_args)
            self.label.grid(row=0, column=0, sticky=(tk.W + tk.E))
            input_args["textvariable"] = input_var
            
        self.input=input_class(self,**input_args)
        self.input.grid(row=0,column=1,sticky=(tk.W+tk.E))                
        self.columnconfigure(0, weight=1)
        
        
    #def grid(self,sticky=(tk.E+tk.W),**kwargs):
    #    super().grid(sticky=sticky,**kwargs)
        
        
    #def getInput():
    #    return 1
        
        
    def get(self):
        
        try:
            if self.variable:
                return self.variable.get()
            elif type(self.input) == tk.Text:
                return self.input.get('1.0', tk.END)
            else:
                return self.input.get()
        except (TypeError, tk.TclError):
            # happens when numeric fields are empty.
            return ''
            
            
    def set(self,value,*args,**kwargs):
        
        if type(self.variable) == tk.BooleanVar:            
            self.variable.set(bool(value))

        elif self.variable:
                self.variable.set(value, *args, **kwargs)
        elif type(self.input) in (ttk.Checkbutton, 
        ttk.Radiobutton):
            if value:
                self.input.select()
            else:
                self.input.deselect()
        elif type(self.input) == tk.Text:
            self.input.delete('1.0', tk.END)
            self.input.insert('1.0', value)
        else: # input must be an Entry-type widget with no variable
            self.input.delete(0, tk.END)
            self.input.insert(0, value)            



'''
 * If we have a variable of class BooleanVar, cast value to bool and set it. BooleanVar.set() will only take a bool, not other falsy or truthy values. This ensures our variable only gets an actual boolean value.
 * If we have any other kind of variable, just pass value to its .set() method.
 * If we have no variable, and a button-style class, we use the .select() and .deselect() methods to select and deselect the button based on the truthy value of the variable.
 * If it's a tk.Text class, we can use its .delete and .insert methods.
 * Otherwise, we use the .delete and .insert methods of input, which work on the Entry, Spinbox, and Combobox classes. We have to do this separately from the tk.Text inputs, because the indexing values work differently.
'''
