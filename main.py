from tkinter import *
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm

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

def update_iteration_info(method, iteration, diff_norm, solution_error, alpha_k=None):
    if method == 'Якоби':
        label_method_jacobi['text'] = f"Метод: {method}"
        label_iteration_jacobi['text'] = f"Итерация: {iteration}"
        label_diff_norm_jacobi['text'] = f"Норма разности: {diff_norm:.16f}"
        label_solution_error_jacobi['text'] = f"Ошибка решения: {solution_error:.16f}"
    else:
        label_method_desc['text'] = f"Метод: {method}"
        label_iteration_desc['text'] = f"Итерация: {iteration}"
        label_diff_norm_desc['text'] = f"Норма разности: {diff_norm:.16f}"
        label_solution_error_desc['text'] = f"Ошибка решения: {solution_error:.16f}"
        if alpha_k is not None:
            label_alpha_desc['text'] = f"Alpha_k: {alpha_k:.16f}"
        else:
            label_alpha_desc['text'] = ""

def plot_fast_desc_method():
    global ani, a, b
    h = float(entry_step.get())
    tau = float(entry_step.get())

    n = (int)(1 // h) + 1
    eps = float(entry_error.get())

    u = np.zeros((n, n))
    next_u = np.zeros((n, n))
    prev_u = np.zeros((n, n))
    
    for i in range(n):
        next_u[i][0] = real_u(i*h, 0)
        next_u[i][n-1] = real_u(i*h, 1)
        next_u[0][i] = real_u(0, i*tau)
        next_u[n-1][i] = real_u(1, i*tau)
    
    fig = plt.figure(3)
    ax = fig.add_subplot(111)
    plt.title("Метод скорейшего спуска")
    
    img = ax.imshow(u.T, origin='lower', cmap='coolwarm', extent=(0, 1, 0, 1))
    
    r = np.zeros((n, n))
    diff_norm = float('+inf')
    solution_error = float('+inf')
    iteration = 0
    
    def update(frame):
        nonlocal u, next_u, prev_u, r, diff_norm, solution_error, iteration
        
        if diff_norm <= eps:
            return img,
        
        prev_u = u.copy()
        
        for i in range(1, n - 1):
            for j in range(1, n - 1):
                r[i][j] = (a * (u[i+1][j] - 2*u[i][j] + u[i-1][j]) / h**2 + 
                          b * (u[i][j+1] - 2*u[i][j] + u[i][j-1]) / tau**2 + 
                          f(i*h, j*tau))
        
        numer = np.sum(r[1:-1, 1:-1]**2)
        Ar = np.zeros_like(r)
        for i in range(1, n-1):
            for j in range(1, n-1):
                Ar[i,j] = (a * (r[i+1][j] - 2*r[i][j] + r[i-1][j]) / h**2 + b * (r[i][j+1] - 2*r[i][j] + r[i][j-1]) / tau**2)
        denom = np.sum(r[1:-1, 1:-1] * Ar[1:-1, 1:-1])
        alpha_k = numer / denom if denom != 0 else 0
        
        next_u[1:-1, 1:-1] = u[1:-1, 1:-1] - alpha_k * r[1:-1, 1:-1]
        u = next_u.copy()
        
        diff_norm = np.linalg.norm(u - prev_u)
        
        # Вычисляем ошибку по сравнению с точным решением
        exact_solution = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                exact_solution[i][j] = real_u(i*h, j*tau)
        #solution_error = np.linalg.norm(u - exact_solution)
        solution_error = np.max(np.abs(u - exact_solution))
        
        iteration += 1
        
        update_iteration_info("Скорейшего спуска", iteration, diff_norm, solution_error, alpha_k)
        
        skip_value = int(entry_skip.get())
        if iteration % skip_value != 0:
            return img,

        img.set_array(u.T)
        img.set_clim(np.min(u), np.max(u)) 
        
        
        return img,
    
    ani = FuncAnimation(
        fig, 
        update, 
        frames=range(1000),
        interval=100,      
        blit=True,         
        repeat=False       
    )
    
    plt.tight_layout()
    plt.show()

def plot_jacobi_method():
    global ani, a, b
    h = float(entry_step.get())
    tau = float(entry_step.get())

    n = (int)(1 // h) + 1
    eps = float(entry_error.get())

    u = np.zeros((n, n))
    next_u = np.zeros((n, n))
    prev_u = np.zeros((n, n))
    
    for i in range(n):
        next_u[i][0] = real_u(i * h, 0)
        next_u[i][n - 1] = real_u(i * h, 1)
        next_u[0][i] = real_u(0, i * tau)
        next_u[n - 1][i] = real_u(1, i * tau)
    
    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    plt.title("Метод Якоби")
    
    img = ax.imshow(u.T, origin='lower', cmap='coolwarm', extent=(0, 1, 0, 1))
    
    diff_norm = float('+inf')
    solution_error = float('+inf')
    iteration = 0
    
    def update(frame):
        nonlocal u, next_u, prev_u, diff_norm, solution_error, iteration
        
        if diff_norm <= eps:
            return img,
        
        prev_u = u.copy()
        
        for i in range(1, n - 1):
            for j in range(1, n - 1):
                x_i = h * i
                t_j = tau * j
                next_u[i][j] = (f(x_i, t_j) * h * h * tau * tau + 
                               a * tau * tau * (prev_u[i + 1][j] + prev_u[i - 1][j]) + 
                               b * h * h * (prev_u[i][j + 1] + prev_u[i][j - 1])) / (2 * (a * tau * tau + b * h * h))
        
        u = next_u.copy()
        
        diff_norm = np.linalg.norm(u - prev_u)
        
        # Вычисляем ошибку по сравнению с точным решением
        exact_solution = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                exact_solution[i][j] = real_u(i*h, j*tau)
        #solution_error = np.linalg.norm(u - exact_solution)
        solution_error = np.max(np.abs(u - exact_solution))
        
        iteration += 1

        update_iteration_info("Якоби", iteration, diff_norm, solution_error)

        skip_value = int(entry_skip.get())
        if iteration % skip_value != 0:
            return img,
        
        img.set_array(u.T)
        img.set_clim(np.min(u), np.max(u))
        
        return img,
    
    ani = FuncAnimation(
        fig, 
        update, 
        frames=range(1000),
        interval=100,      
        blit=True,         
        repeat=False       
    )
    
    plt.tight_layout()
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

# Create main window
root = Tk()
root.geometry('900x500')
frm = ttk.Frame(root, padding=10)
frm.grid()

root.title("Вычислительный практикум решение задачи Дирихле")

# Parameters frame
params_frame = ttk.LabelFrame(frm, text="Параметры", padding=10)
params_frame.grid(column=0, row=0, padx=5, pady=5, sticky="nsew")

label_step = ttk.Label(params_frame, text='Шаг сетки:')
label_step.grid(column=0, row=0, sticky="w")

entry_step = ttk.Entry(params_frame)
entry_step.insert(0, '0.01')
entry_step.grid(column=1, row=0, padx=5)

label_error = ttk.Label(params_frame, text='Эпсилон:')
label_error.grid(column=0, row=1, sticky="w")

entry_error = ttk.Entry(params_frame)
entry_error.insert(0, '0.01')
entry_error.grid(column=1, row=1, padx=5)

label_skip = ttk.Label(params_frame, text='Отрисовывать каждую итерацию кратную:')
label_skip.grid(column=0, row=2, sticky="w")

entry_skip = ttk.Entry(params_frame)
entry_skip.insert(0, '10')
entry_skip.grid(column=1, row=2, padx=5)

# Eigenvalues frame
eigen_frame = ttk.LabelFrame(frm, text="Собственные значения", padding=10)
eigen_frame.grid(column=0, row=1, padx=5, pady=5, sticky="nsew")

label_min_eigen = ttk.Label(eigen_frame, text='Мин. собственное число:')
label_min_eigen.grid(column=0, row=0, sticky="w")

label_min_eigen_value = ttk.Label(eigen_frame, text='0')
label_min_eigen_value.grid(column=1, row=0)

label_max_eigen = ttk.Label(eigen_frame, text='Макс. собственное число:')
label_max_eigen.grid(column=0, row=1, sticky="w")

label_max_eigen_value = ttk.Label(eigen_frame, text='0')
label_max_eigen_value.grid(column=1, row=1)

# Iteration info frame
iter_frame_jacobi = ttk.LabelFrame(frm, text="Информация о решении Якоби", padding=10)
iter_frame_jacobi.grid(column=0, row=2, padx=5, pady=5, sticky="nsew")

label_method_jacobi = ttk.Label(iter_frame_jacobi, text="Метод: ")
label_method_jacobi.grid(column=0, row=0, sticky="w")

label_iteration_jacobi = ttk.Label(iter_frame_jacobi, text="Итерация: ")
label_iteration_jacobi.grid(column=0, row=1, sticky="w")

label_diff_norm_jacobi = ttk.Label(iter_frame_jacobi, text="Норма разности: ")
label_diff_norm_jacobi.grid(column=0, row=2, sticky="w")

label_solution_error_jacobi = ttk.Label(iter_frame_jacobi, text="Ошибка решения: ")
label_solution_error_jacobi.grid(column=0, row=3, sticky="w")

iter_frame_desc = ttk.LabelFrame(frm, text="Информация о решении Скорейшего спуска", padding=10)
iter_frame_desc.grid(column=1, row=2, padx=5, pady=5, sticky="nsew")

label_method_desc = ttk.Label(iter_frame_desc, text="Метод: ")
label_method_desc.grid(column=0, row=0, sticky="w")

label_iteration_desc = ttk.Label(iter_frame_desc, text="Итерация: ")
label_iteration_desc.grid(column=0, row=1, sticky="w")

label_diff_norm_desc = ttk.Label(iter_frame_desc, text="Норма разности: ")
label_diff_norm_desc.grid(column=0, row=2, sticky="w")

label_solution_error_desc = ttk.Label(iter_frame_desc, text="Ошибка решения: ")
label_solution_error_desc.grid(column=0, row=3, sticky="w")

label_alpha_desc = ttk.Label(iter_frame_desc, text="")
label_alpha_desc.grid(column=0, row=4, sticky="w")

# Buttons frame
buttons_frame = ttk.Frame(frm, padding=10)
buttons_frame.grid(column=0, row=3, padx=5, pady=5, sticky="nsew")

ttk.Button(buttons_frame, text='Точное решение', command=plot_real_function).grid(column=0, row=0, padx=5)
ttk.Button(buttons_frame, text='Метод Якоби', command=plot_jacobi_method).grid(column=1, row=0, padx=5)
ttk.Button(buttons_frame, text='Метод Скорейшего спуска', command=plot_fast_desc_method).grid(column=2, row=0, padx=5)
ttk.Button(buttons_frame, text='Найти собственные значения', command=find_eigen_values).grid(column=3, row=0, padx=5)

root.mainloop()
