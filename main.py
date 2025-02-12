from tkinter import *
from tkinter import ttk

import numpy as np

import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.ticker import LinearLocator

def u(X, Y):
    return np.exp(X) * np.cos(Y)

def plot_real_function():
    step = float(entry_step.get())
    X = np.arange(0, 1, step)
    Y = np.arange(0, 1, step)
    X, Y = np.meshgrid(X, Y)
    Z = u(X, Y)

    plt.imshow(Z, origin='lower', cmap=cm.coolwarm, extent=(0,1,0,1))
    
    plt.title("Точное решение")

    plt.show()

root = Tk()
root.geometry('600x400')
frm = ttk.Frame(root, padding=10)
frm.grid()

label_step = ttk.Label(frm, text='Шаг сетки:')
label_step.grid(column = 0, row = 0)

entry_step = ttk.Entry(frm)
entry_step.insert(0, '0.01')
entry_step.grid(column = 1, row = 0)

label_error = ttk.Label(frm, text='Ошибка:')
label_error.grid(column = 0, row = 1)

entry_error = ttk.Entry(frm)
entry_error.insert(0, '0.01')
entry_error.grid(column = 1, row = 1)


label_min_eigen = ttk.Label(frm, text='Мин. собственное число:')
label_min_eigen.grid(column = 0, row = 2)

label_min_eigen_value = ttk.Label(frm, text='0')
label_min_eigen_value.grid(column = 1, row = 2)


label_max_eigen = ttk.Label(frm, text='Макс. собственное число:')
label_max_eigen.grid(column = 0, row = 3)

label_max_eigen_value = ttk.Label(frm, text='0')
label_max_eigen_value.grid(column = 1, row = 3)

ttk.Button(frm, text='Показать точное решение', command=plot_real_function).grid(column = 2, row = 5)
root.mainloop()
