import tkinter as tk
from tkinter import Label, Entry, Tk, Button
from tkinter import messagebox
import numpy as np
from matplotlib.figure import Figure


import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

def plot():
    x = ['Col A', 'Col B', 'Col C']
    
    y = [50, 20, 80]
    
    fig = plt.figure(figsize=(4, 5))
    plt.bar(x=x, height=y)
    
    # You can make your x axis labels vertical using the rotation
    plt.xticks(x, rotation=90)
    
    # specify the window as master
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().grid(row=2, column=3, ipadx=40, ipady=20)
    
    # navigation toolbar
    #toolbarFrame = tk.Frame(master=window)
    #toolbarFrame.grid(row=2,column=2)
    #oolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
    

window = Tk()
window.title("Visual PSV")
#window.geometry('600x300')
plot()
lon_lb = Label(
   window,
   text="Longitude ",
)
lon_lb.grid(row=1, column=1)

lon_tf = Entry(
    window
)
lon_tf.grid(row = 2, column = 1)

lat_lb = Label(
   window,
   text="Latitude ",
)
lat_lb.grid(row=3, column=1)

lat_tf = Entry(
    window
)
lat_tf.grid(row = 4, column = 1)


cal_bt = Button(
   window, #Заготовка с настроенными отступами.
   text='Calcualte', #Надпись на кнопке.
   command=plot
   
)
cal_bt.grid(row = 5, column = 1)



#this piece of shit must be in the end
window.mainloop()