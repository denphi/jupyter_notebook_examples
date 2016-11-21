""" Various helper functions and resources for blog:
    http://www.soton.ac.uk/~ol1g13
    -----------------------------------------------
    Copyright (C) 2014  Oliver W. Laslett
    <O.Laslett@soton.ac.uk>

    See LICENSE file in project home folder ./LICENSE for
    the full copyright notice.
"""

from __future__ import division
from math import pi, exp, sin, cos, acos, asin
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt


def SWenergy(tm, h, th):
    """ Calculates the free energy for a Stoner-Wohlfarth particle
    with uniaxial ansiotropy. Given the angle between the positive
    z-axis and the magnetic moment, it returns this energy. """
    return 0.5*(np.sin(tm) ** 2) - h*np.cos(tm-th)


def minMaxSW(h, th):
    """ Returns the easy and hard axis for an SW particle with
    uniaxial anisotropy.
    Result is four tuples of the form (angle, energy)
    in the order of min(min) max(min) min(max) max(max)."""

    # Use a brent search with bounds starting around the anisotropy axis
    min1 = op.brent(SWenergy, brack=(0.0, 0.001),
                    args = (h, th))
    min2 = op.brent(SWenergy, brack=(pi, pi+0.001),
                    args = (h, th))

    min1 = [min1, SWenergy(min1, h, th)]
    min2 = [min2, SWenergy(min2, h, th)]

    # Check whether the minima are the same solution
    if min1[1] == min2[1]:
        max1 = op.brent(lambda x: SWenergy(x, h, th)*(-1),
                        brack=(min1[0], min1[0]+0.01))

        max1 = (max1, SWenergy(max1, h, th))
        max2 = max1
    else:
        try:
            max1 = op.brent(lambda x: SWenergy(x, h, th)*(-1),
                            brack=(min(min1[0], min2[0]), (min1[0]+min2[0])/2,
                                   max(min1[0], min2[0])))
        except ValueError:
            max1 = op.minimize(lambda x: SWenergy(x,h,th)*(-1), (min1[0]+min2[0])/2.,
                               method='TNC',
                               bounds=[(min(min1[0],min2[0]), max(min1[0],min2[0]))])
            if max1.success:
                max1 = max1.x
            else:
                raise RuntimeError
        if max1<0: max1+=2*pi
        max1 = (max1, SWenergy(max1, h, th))
        
        try:
            max2 = op.brent(lambda x: SWenergy(x, h, th)*(-1),
                            brack=(min(min1[0], min2[0])+2*pi, (min1[0]+min2[0]+2*pi)/2,
                                   max(min1[0], min2[0])))
        except ValueError:
            max2 = op.minimize(lambda x: SWenergy(x, h, th)*(-1), (min1[0]+min2[0]+2*pi)/2.,
                               method='TNC',
                               bounds=[(max(min1[0],min2[0]), min(min1[0],min2[0])+2*pi)])
            if max2.success:
                max2 = max2.x
            else:
                max2 = max1[0]
        if max2<0: max2+=2*pi
        max2 = (max2, SWenergy(max2, h, th))

    if min1[0]<0: min1[0]+=2*pi
    if min2[0]<0: min2[0]+=2*pi

    tuple_energy = lambda x: x[1]
    minmin = min([min1, min2], key=tuple_energy)
    maxmin = max([min1, min2], key=tuple_energy)
    minmax = min([max1, max2], key=tuple_energy)
    maxmax = max([max1, max2], key=tuple_energy)
    
    return minmin, maxmin, minmax, maxmax

def analytic_minMaxSW(h, orientation):
    """ Uses the analytic result to compute the energy barriers
    for the aligned and perpendicular applied field cases. """
    sols = []
    mins = []
    maxs = []
    if orientation=='align':
        sols.append(0)
        sols.append(pi)
        if -1<h<1:
            sols.append(acos(-h))
            sols.append(2*pi - acos(-h))

        swe = lambda x: SWenergy(x, h, 0)
        d2e = lambda x: cos(2*x)+h*cos(x)

    elif orientation=='perp':
        sols.append(pi/2.)
        sols.append(3*pi/2.)
        if -1<h<1:
            sols.append(asin(h))
            sols.append(pi-asin(h))

        swe = lambda x: SWenergy(x, h, pi/2)
        d2e = lambda x: cos(2*x)+h*sin(x)

    else:
        error("orientation must be 'align' or 'perp'")

    for s in sols:
        if d2e(s) > 0:
            mins.append((s, swe(s)))
        else:
            maxs.append((s, swe(s)))

    tuple_energy = lambda x: x[1]
    min1 = min(mins, key = tuple_energy)
    min2 = max(mins, key = tuple_energy)
    max1 = min(maxs, key = tuple_energy)
    max2 = max(maxs, key = tuple_energy)

    return min1, min2, max1, max2


def drawSWsketch():
    """ Used in IPython notebook to draw a Stoner-Wohlfarth
        particle sketch. """
    fg = plt.figure(figsize=(5,5))
    ax = fg.add_subplot(111, polar=True)
    th = np.linspace(0,2*pi,10000)
    ax.plot([pi/2, 3*pi/2], [1.1, 1.1], 'k-', label='x-axis')
    ax.plot(th, 0.9*np.ones(10000), 'b-', label='Particle')
    ax.plot(np.linspace(0, pi/3, 1000), 0.3*np.ones(1000), 'm-')
    ax.plot(np.linspace(0, pi/6, 500), 0.5*np.ones(500), 'b-')
    ax.arrow(pi/3, 0, 0, 1.1, label='Applied field', fc='m', ec='m', head_width=0.1)
    ax.arrow(0, 0, 0, 1.1, fc='r', ec='r', head_width=0.1)
    ax.arrow(pi, 0, 0, 1.1, fc='r', ec='r', head_width=0)
    ax.arrow(pi/6, 0, 0, 1.1, fc='b', ec='b', head_width=0.1, width=0.02)
    ax.text(pi/3, 1.3, 'Applied field', fontsize=14)
    ax.text(-0.5, 1.5, 'Anisotropy axis (z-axis)', fontsize=14)
    ax.text(pi/6, 1.3, 'Magnetic moment vector', fontsize=14)
    ax.text(pi/5, 0.35, '$\\theta_H$', fontsize=14)
    ax.text(pi/12, 0.55, '$\\theta$', fontsize=14)
    lg = ax.legend(loc=(0.8, 0.1), fontsize=16)
    ax.set_theta_zero_location(loc='N')
    ax.set_theta_direction(-1)
    ax.set_rmax(1.3)
    ax.set_axis_off()

def eswidget(h, th):
    # Compute the free energy of a Stoner-Wohlfarth particle at various magnetic moment angles.
    tm = np.linspace(0,2*pi,10000)
    em = SWenergy(tm, h, th)
    
    # Use the energyscape search to find the easy and hard axes
    minaxes = minMaxSW(h, th)
    
    # Plot the results
    fg = plt.figure(figsize=(12,6))
    ax = fg.add_subplot(111)
    ax.plot(tm, em, 'b-', lw=0.5, label='Free energy')
    ax.plot([x[0] for x in minaxes], [x[1] for x in minaxes], 'r.', ms=7)
    ybounds = (-0.8,1.5)
    ax.plot([th, th], ybounds, 'g-', lw=1.2, label='Applied field angle')
    ax.plot(5*[th], np.linspace(ybounds[0], ybounds[1], 5), 'gv', ms=5)
    ax.plot(2*[pi], ybounds, 'm-', lw=0.8, label='Anisotropy axis')
    ax.plot(5*[pi], np.linspace(ybounds[0], ybounds[1], 5), 'mv', ms=5)
    ax.plot(2*[0], ybounds, 'm-', lw=0.8)
    ax.plot(5*[0], np.linspace(ybounds[0], ybounds[1], 5), 'mv', ms=5)
    ax.plot(2*[2*pi], ybounds, 'm-', lw=0.8)
    ax.plot(5*[2*pi], np.linspace(ybounds[0], ybounds[1], 5), 'mv', ms=5)
    ax.set_ylim((-0.8,1.5))
    ax.set_xlim((0,2*pi))
    lg = ax.legend(loc='upper right')
