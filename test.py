import numpy as np
from matplotlib import pyplot as plt
import time

n = 10
dt = 0.0001
diff = 0.01
visc = 1
size = (n + 2)**2
water = np.zeros((n + 2, n + 2), dtype = np.float64)
past = water.copy()
source = water.copy()
source[::, 5] = 1
image = water[1::n, 1::n]
u = water.copy()
v = water.copy()
u0 = water.copy()
v0 = water.copy()

def swap(a, b):
	tmp = a
	a = b
	b = tmp


def add_source(water, source, dt):
	water += source * dt



def diffuse_v1(b, water, past, diff, dt, n):
	a = dt * diff * n * n
	for k in range(20):
		for i in range(1, n + 1):
			for j in range(1, n + 1):
				water[i, j] = past[i, j] + a * (water[i - 1, j] + water[i + 1, j] + water[i, j - 1] + water[i, j + 1]) / (4 * a)
		set_bnd(n, b, water)



def advect(n, b, water, past, u, v, dt):
	dt0 = dt * n
	for i in range(1, n + 1):
		for j in range(1, n + 1):
			x = i - dt0 * u[i, j]
			y = j - dt0 * v[i, j]
			if x < 0.5:
				x = 0.5
			else if x > n + 0.5:
				x = n + 0.5
			i0 = int(x)
			i1 = i0 + 1;
			if y < 0.5:
				y = 0.5
			else if y > n + 0.5:
				y = n + 0.5
			j0 = int(y)
			j1 = j0 + 1;
			s1 = x - i0
			s0 = 1 - s1
			t1 = y - j0
			t0 = 1 - t1
			water[i, j] = s0 * (t0 * past[i0, j0] + t1 * past[i0, j1]) + s1 * (t0 * past[i1, j0] + t1 * past[i1, j1])
	set_bnd(n, b, water)

def vel_step(n, u, v, u0, v0, visc, dt):
	add_source(u, u0, dt)
	add_source(v, v0, dt)
	swap(u0, u)
	diffuse(1, u, u0, visc, dt)
	swap(v0, v)
	diffuse(2, v, v0, visc, dt)
	project(n, u, v, u0, v0)
	swap(u0, u)
	swap(v0, v)
	advect(n, u, v, u0, v0)


def projet(n, u, v, p, div):
	h = 1. / n
	for i in range(1, n + 1):
		for j in range(1, n + 1):
			div[i, j] = -0.5 * h * (u[i + 1, j] - u[i - 1, j] + v[i, j + 1] - v[i, j - 1])
			p[i, j] = 0
	set_bnd(n, 0, div)
	set_bnd(n, 0, p)

	for k in range(20):
		for i in range(1, n + 1):
			for j in range(1, n + 1):
				p[i, j] = (div[i, j] + p[i - 1, j] + p[i + 1, j] + p[i, j - 1] + p[i, j + 1]) / 4
		set_bnd(n, 0, p)

	for i in range(1, n + 1):
		for j in range(1, n + 1):
			u[i, j] -= 0.5 * (p[i + 1, j] - p[i - 1, j]) / h
			v[i, j] -= 0.5 * (p[i + 1, j] - p[i - 1, j]) / h
	set_bnd(n, 1, u)
	set_bnd(n, 2, v)

def set_bnd(n, b, x):
	for i in range(1, n + 1):
		x[0, i] = (-x[1, i] if b == 1 else x[1, i])
		x[n + 1, i] = (-x[n, i] if b == 1 else x[n, i])
		x[i, 0] = (-x[i, 1] if b == 2 else x[i, 1])
		x[i, n + 1] = (-x[i, n] if b == 2 else x[i, n])

	x[0, 0] = 0.5 * (x[1, 0] + x[0, 1])
	x[0, n + 1] = 0.5 * (x[1, n + 1] + x[0, n])
	x[n + 1, 0] = 0.5 * (x[n, 0] + x[n + 1, 1])
	x[n + 1, n + 1] = 0.5 * (x[n, n + 1] + x[n + 1, n])



for i in range(100):

	print water, "\n"
	plt.imshow(water)
	plt.colorbar()
	plt.show()

plt.imshow(water)
plt.colorbar()
plt.show()

