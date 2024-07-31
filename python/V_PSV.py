import tkinter as tk
from tkinter import Label, Entry, Tk, Button
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,  
NavigationToolbar2Tk) 
  
# plot function is created for  
# plotting the graph in  
# tkinter window 
def plot(): 
  
    # the figure that will contain the plot 
    fig = Figure(figsize = (5, 5), 
                 dpi = 100) 
  
    # list of squares 
    y = [i**2 for i in range(101)] 
  
    # adding the subplot 
    ax1 = fig.add_subplot(211) 
  
    # plotting the graph 
    ax1.plot(y) 
    
    y = [i**(0.5) for i in range(101)] 
    ax2 = fig.add_subplot(212) 
    ax2.plot(y) 
  
    # creating the Tkinter canvas 
    # containing the Matplotlib figure 
    canvas = FigureCanvasTkAgg(fig, 
                               master = window)   
    canvas.draw() 
  
    # placing the canvas on the Tkinter window 
    canvas.get_tk_widget().pack() 
  
    # creating the Matplotlib toolbar 
    toolbar = NavigationToolbar2Tk(canvas, 
                                   window) 
    toolbar.update() 
  
    # placing the toolbar on the Tkinter window 
    canvas.get_tk_widget().pack(expand=True, fill=tk.BOTH) 
  
# the main Tkinter window 
window = Tk() 
  
# setting the title  
window.title('Plotting in Tkinter') 
  
# dimensions of the main window 
window.geometry("500x500") 
  
lon_lb = Label(
   window,
   text="Longitude ",
)
lon_lb.pack()

lat_lb = Label(
   window,
   text="Longitude ",
)
lat_lb.pack()

# button that displays the plot 
plot_button = Button(master = window,  
                     command = plot, 
                     height = 2,  
                     width = 10, 
                     text = "Calculate") 
  
# place the button  
# in main window 
plot_button.pack(side=tk.LEFT) 
  
# run the gui 
window.mainloop() 