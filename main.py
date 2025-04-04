from tkinter import *
from tkinter import ttk

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

a = 1
b = 1.2

def real_u(X, Y):
    return np.exp(X) * np.cos(Y)

def f(x, y):
    return 0.2 * np.exp(x) * np.cos(y)

def find_eigen_values():
    global a, b
    h = float(entry_step.get())
    tau = float(entry_step.get())

    n = (int)(1 // h) + 1

    eps = float(entry_error.get())

    u = np.ones((n, n))
    for i in range(n):
        u[i][0] = 0
        u[i][n-1] = 0
        u[0][i] = 0
        u[n-1][i] = 0

    u = u / np.linalg.norm(u)

    v = np.empty((n, n))

    prev_eigen = float('inf')
    max_eigen_value = 0

    while abs(prev_eigen - max_eigen_value) >= eps:
        for i in range(1, n - 1):
            for j in range(1, n - 1):
                v[i][j] = -a*(u[i+1][j] - 2 * u[i][j] + u[i-1][j])/(h*h) - b*(u[i][j+1] - 2 * u[i][j] + u[i][j-1])/(tau*tau)
        
        prev_eigen = max_eigen_value
        max_eigen_value = 0
        for i in range(n):
            for j in range(n):
                max_eigen_value += u[i][j] * v[i][j]
        
        u = v / np.linalg.norm(v)

    B_max_eigen_value = 0

    prev_eigen = float('inf')
    B_max_eigen_value = 0

    while abs(prev_eigen - B_max_eigen_value) >= eps:
        for i in range(1, n - 1):
            for j in range(1, n - 1):
                v[i][j] = (max_eigen_value + 1) * u[i][j] - (-a*(u[i+1][j] - 2 * u[i][j] + u[i-1][j])/(h*h) - b*(u[i][j+1] - 2 * u[i][j] + u[i][j-1])/(tau*tau))
        
        prev_eigen = B_max_eigen_value
        B_max_eigen_value = 0
        for i in range(n):
            for j in range(n):
                B_max_eigen_value += u[i][j] * v[i][j]

        u = v / np.linalg.norm(v)

    min_eigen_value = max_eigen_value - B_max_eigen_value

    label_max_eigen_value['text'] = str(max_eigen_value)
    label_min_eigen_value['text'] = str(min_eigen_value)

ani = []

def plot_fast_desc_method():
    global ani, a, b
    h = float(entry_step.get())
    tau = float(entry_step.get())

    n = (int)(1 // h) + 1
    eps = float(entry_error.get())

    u = np.zeros((n, n))
    next_u = np.zeros((n, n))
    for i in range(n):
        next_u[i][0] = real_u(i*h, 0)
        next_u[i][n-1] = real_u(i*h, 1)
        next_u[0][i] = real_u(0, i*tau)
        next_u[n-1][i] = real_u(1, i*tau)
    
    artists = []

    fig = plt.figure(2)
    plt.title("Метод скорейшего спуска")

    r = np.zeros((n, n))
    norm_r = float('+inf')
    while (norm_r > eps):
        for i in range(1, n - 1):
            for j in range(1, n - 1):
                r[i][j] = a * (u[i+1][j] - 2*u[i][j] + u[i-1][j]) / h + b * (u[i][j+1] - 2*u[i][j] + u[i][j-1]) / tau + f(i * h, j * tau)

        numer = 0
        denom = 0
        for i in range(1, n - 1):
            for j in range(1, n - 1):
                numer += r[i][j] * r[i][j]
                denom += r[i][j] * (a * (r[i+1][j] - 2*r[i][j] + r[i-1][j]) / h + b * (r[i][j+1] - 2*r[i][j] + r[i][j-1]) / tau)

        alpha_k = numer/denom

        for i in range(1, n - 1):
            for j in range(1, n - 1):
                x_i = h * i
                t_j = tau * j

                next_u[i][j] = u[i][j] - alpha_k * r[i][j]
        u = next_u
        artists.append([plt.imshow(u.T, origin='lower', cmap='coolwarm', extent=(0,1,0,1))])
        norm_r = np.linalg.norm(r)

    
    ani = animation.ArtistAnimation(fig, artists=artists, interval=100)

    plt.show()


def plot_jacobi_method(): 
    global ani, a, b
    h = float(entry_step.get())
    tau = float(entry_step.get())

    n = (int)(1 // h) + 1


    u = np.zeros((n, n))
    next_u = np.zeros((n, n))
    for i in range(n):
        next_u[i][0] = real_u(i*h, 0)
        next_u[i][n-1] = real_u(i*h, 1)
        next_u[0][i] = real_u(0, i*tau)
        next_u[n-1][i] = real_u(1, i*tau)
    
    artists = []

    fig = plt.figure(2)
    plt.title("Метод Якоби")

    #TODO: Вместо 1000 итераций, условие что ||u_k - u(x,y)|| < eps
    #                                        текущее  точное
    for k in range(1000):
        for i in range(1, n - 1):
            for j in range(1, n - 1):
                x_i = h * i
                t_j = tau * j
                next_u[i][j] = (f(x_i, t_j)*h*h*tau*tau + a*tau*tau*(u[i+1][j] + u[i-1][j]) + b*h*h*(u[i][j+1] + u[i][j-1])) / (2*(a*tau*tau + b*h*h))
        u = next_u
        artists.append([plt.imshow(u.T, origin='lower', cmap='coolwarm', extent=(0,1,0,1))])

    #TODO: Отрисовывать in time
    #TODO: На каждой итерации выводить ее номер, невязку, ошибку, погрешность
    ani = animation.ArtistAnimation(fig, artists=artists, interval=100)

    plt.show()

def plot_real_function():
    step = float(entry_step.get())
    X = np.arange(0, 1, step)
    Y = np.arange(0, 1, step)
    X, Y = np.meshgrid(X, Y)
    Z = real_u(X, Y)

    plt.figure(1)
    plt.imshow(Z, origin='lower', cmap='coolwarm', extent=(0,1,0,1))
    
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

#TODO: Сделать переключатель метода остановки (ошибка, погрешность, количество итераций)
label_error = ttk.Label(frm, text='Эпсилон:')
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

ttk.Button(frm, text='Показать точное решение', command=plot_real_function).grid(column = 1, row = 5)
ttk.Button(frm, text='Метод Якоби', command=plot_jacobi_method).grid(column = 2, row = 5)
ttk.Button(frm, text='Метод Скорейшего спуска', command=plot_fast_desc_method).grid(column = 3, row = 5)
ttk.Button(frm, text='Найти собственные значения', command=find_eigen_values).grid(column = 1, row = 6)
root.mainloop()
