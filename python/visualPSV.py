import tkinter as tk
from tkinter import Label, Entry, Tk
from tkinter import messagebox

window = Tk()
window.title("Visual PSV")
window.geometry('1200x900')

lon_lb = Label(
   window,
   text="Longitude ",
)
lon_lb.grid(row=4, column=1)

lon_tf = Entry(
    window
)
lon_tf.grid(row = 5, column = 1)

lat_lb = Label(
   window,
   text="Latitude ",
)
lat_lb.grid(row=6, column=1)

lat_tf = Entry(
    window
)
lat_tf.grid(row = 7, column = 1)

#this piece of shit must be in the end
window.mainloop()